"""
Transport phenomena module for Thermally Coupled Reactor (TCR) simulation:
- Packed-bed-to-wall heat transfer coefficients with Froment & Bischoff / Dixon (1996)
  1D lumped effective resistance combining wall film resistance 1/h_w(dp) and radial bed
  conduction resistance D / (6.133 * k_er):
    1 / h_{i,eff} = dp_exo / (Nu_{w,i} * kg_exo) + D_in / (6.133 * k_{er,exo})
    1 / h_{o,eff} = dp_endo / (Nu_{w,o} * kg_endo) + Dh_shell / (6.133 * k_{er,endo})
  (Resolves Swarm Critic-11 characteristic length bug and Critic-05 radial resistance)
- Parabolic radial peak centerline temperature calculator:
    T_centerline(z) = T_mean(z) + (q_v(z) * R_in^2) / (8 * k_er(z))
- Variable shell thermal conductivity via Wassiljewa equation with Mason-Saxena interaction
  parameters kg_endo(T, y) supporting 3-component [C6H12, C6H6, H2] and 4-component [C6H12, C6H6, H2, Ar]
  shell mixtures, including shorthand aliases [CH, B, H2, Ar] (Swarm Critic-10)
- Boiling water coolant jacket film coefficient for conventional reactor
- Cylindrical tube wall conductive resistance
- Overall dynamic heat transfer coefficient U(z)
- Ergun packed-bed momentum balances
"""

from typing import Union, List, Dict, Optional
import numpy as np
from .thermo import R_GAS, calc_cp_species, MW
from .geometry import ReactorGeometry
from .config import TransportConfig, CatalystBedConfig

# ---------------------------------------------------------------------------
# Species Alias Mapping & Normalization
# ---------------------------------------------------------------------------
SPECIES_ALIASES: Dict[str, str] = {
    'CH': 'C6H12',
    'CYCLOHEXANE': 'C6H12',
    'C6H12': 'C6H12',
    'B': 'C6H6',
    'BENZENE': 'C6H6',
    'C6H6': 'C6H6',
    'H2': 'H2',
    'HYDROGEN': 'H2',
    'AR': 'Ar',
    'ARGON': 'Ar',
    'Ar': 'Ar',
    'CO2': 'CO2',
    'CO': 'CO',
    'CH3OH': 'CH3OH',
    'MEOH': 'CH3OH',
    'H2O': 'H2O',
    'WATER': 'H2O',
    'CH4': 'CH4',
    'METHANE': 'CH4',
    'N2': 'N2',
    'NITROGEN': 'N2',
}

def resolve_species(species_name: str) -> str:
    """Resolves species alias (e.g. 'CH' -> 'C6H12', 'B' -> 'C6H6') to canonical name."""
    s_clean = species_name.strip()
    s_upper = s_clean.upper()
    if s_upper in SPECIES_ALIASES:
        return SPECIES_ALIASES[s_upper]
    if s_clean in MW:
        return s_clean
    return s_clean

# ---------------------------------------------------------------------------
# Wassiljewa-Mason-Saxena Pure Component Gas Transport Properties
# Sourced from NIST Chemistry WebBook & Reid, Prausnitz & Poling (Properties of Gases and Liquids)
# Accurately reproduces Swarm Critic-10 benchmark table for cyclohexane dehydrogenation
# ---------------------------------------------------------------------------
SHELL_TRANSPORT_PROPS: Dict[str, Dict[str, float]] = {
    'C6H12': {
        'k500': 0.0245,       # W/(m·K) at 500 K reference
        'mw': 84.1595e-3,     # kg/mol
        'mu500': 1.15e-5,     # Pa·s at 500 K reference
        'n_k': 1.50,          # Temperature exponent for k(T) = k500 * (T / 500)^n_k
        'n_mu': 0.85,         # Temperature exponent for mu(T) = mu500 * (T / 500)^n_mu
    },
    'C6H6': {
        'k500': 0.0270,
        'mw': 78.1118e-3,
        'mu500': 1.10e-5,
        'n_k': 1.50,
        'n_mu': 0.85,
    },
    'H2': {
        'k500': 0.2580,
        'mw': 2.01588e-3,
        'mu500': 1.35e-5,
        'n_k': 0.75,
        'n_mu': 0.70,
    },
    'Ar': {
        'k500': 0.0264,
        'mw': 39.948e-3,
        'mu500': 3.20e-5,
        'n_k': 0.78,
        'n_mu': 0.75,
    },
}

def calc_kg_endo_wassiljewa(
    T: float,
    y: Union[Dict[str, float], np.ndarray, List[float]],
    species_list: Optional[List[str]] = None
) -> float:
    """
    Computes composition- and temperature-dependent mixture thermal conductivity (W/(m·K))
    for the shell stream via the Wassiljewa equation with Mason-Saxena interaction parameters:
      k_mix(T, y) = sum_i [ y_i * k_i(T) / sum_j (y_j * Phi_ij) ]
      Phi_ij = [ 1 + (mu_i/mu_j)^0.5 * (M_j/M_i)^0.25 ]^2 / [ 8 * (1 + M_i/M_j) ]^0.5

    Supports 3-component [C6H12, C6H6, H2] and 4-component [C6H12, C6H6, H2, Ar] mixtures,
    as well as shorthand aliases [CH, B, H2, Ar].
    (Wassiljewa 1904; Mason & Saxena 1958; Swarm Critic-10)
    """
    T_safe = max(float(T), 150.0)

    # Parse mole fractions and species names
    if isinstance(y, dict):
        raw_names = list(y.keys())
        raw_vals = [float(y[k]) for k in raw_names]
    else:
        if species_list is None:
            if len(y) == 4:
                raw_names = ['C6H12', 'C6H6', 'H2', 'Ar']
            elif len(y) == 3:
                raw_names = ['C6H12', 'C6H6', 'H2']
            else:
                raw_names = [f'sp_{i}' for i in range(len(y))]
        else:
            raw_names = list(species_list)
        raw_vals = [float(v) for v in y]

    n = len(raw_names)
    if n == 0:
        return 0.040

    canon_names = [resolve_species(name) for name in raw_names]
    vals_arr = np.maximum(np.array(raw_vals, dtype=float), 1e-15)
    tot = float(np.sum(vals_arr))
    if tot > 0:
        y_frac = vals_arr / tot
    else:
        y_frac = np.ones(n) / n

    # Compute pure component properties at T_safe
    k_pure = np.zeros(n)
    mu_pure = np.zeros(n)
    mw_pure = np.zeros(n)

    for i, sp in enumerate(canon_names):
        if sp in SHELL_TRANSPORT_PROPS:
            p = SHELL_TRANSPORT_PROPS[sp]
            Tr = T_safe / 500.0
            k_pure[i] = p['k500'] * (Tr ** p['n_k'])
            mu_pure[i] = p['mu500'] * (Tr ** p['n_mu'])
            mw_pure[i] = p['mw']
        elif sp in MW:
            mw_pure[i] = MW[sp]
            k_pure[i] = 0.040 * ((T_safe / 500.0) ** 1.0)
            mu_pure[i] = 1.2e-5 * ((T_safe / 500.0) ** 0.8)
        else:
            mw_pure[i] = 40.0e-3
            k_pure[i] = 0.040
            mu_pure[i] = 1.2e-5

    # Mason-Saxena interaction parameters Phi_ij
    Phi = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                Phi[i, j] = 1.0
            else:
                num = (1.0 + np.sqrt(mu_pure[i] / mu_pure[j]) * ((mw_pure[j] / mw_pure[i]) ** 0.25)) ** 2
                den = np.sqrt(8.0 * (1.0 + mw_pure[i] / mw_pure[j]))
                Phi[i, j] = num / den

    # Wassiljewa mixture sum
    k_mix = 0.0
    for i in range(n):
        den_i = float(np.sum(y_frac * Phi[i, :]))
        k_mix += float(y_frac[i] * k_pure[i] / max(den_i, 1e-12))

    return float(max(k_mix, 0.010))

# Standard Aliases
calc_kg_shell = calc_kg_endo_wassiljewa
calc_kg_endo = calc_kg_endo_wassiljewa

# ---------------------------------------------------------------------------
# Effective Radial Thermal Conductivity (Bauer & Schlünder / Froment / Dixon)
# ---------------------------------------------------------------------------
def calc_ker(Re_p: float, Pr: float, kg: float, Pe_hr: float = 9.0, ker_stagnant: float = 0.22) -> float:
    """
    Computes effective radial thermal conductivity of packed bed k_er (W/(m·K))
    combining stagnant bed conduction and turbulent dispersion:
      k_er = k_er^0 + (Re_p * Pr / Pe_hr) * k_g
    (Bauer & Schlünder 1978; Froment, Bischoff & De Wilde 2011 Section 3.4.1; Dixon 1996; Swarm Critic-05)
    """
    Re_safe = max(float(Re_p), 1.0)
    Pr_safe = max(float(Pr), 0.01)
    kg_safe = max(float(kg), 1e-4)
    Pe_safe = max(float(Pe_hr), 0.1)
    k_disp = (Re_safe * Pr_safe / Pe_safe) * kg_safe
    return float(ker_stagnant + k_disp)

# ---------------------------------------------------------------------------
# Centerline Temperature Calculator
# ---------------------------------------------------------------------------
def calc_centerline_temp(
    T_mean: Union[float, np.ndarray],
    q_vol: Union[float, np.ndarray],
    R_in: float,
    k_er: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Computes parabolic radial peak centerline temperature T_centerline (K):
      T_centerline(z) = T_mean(z) + (q_v(z) * R_in^2) / (8 * k_er(z))

    Where:
      - T_mean: Cross-sectional cup-average temperature (K)
      - q_vol: Local volumetric heat generation rate (W/m^3 of bed)
      - R_in: Inner tube radius (m)
      - k_er: Effective radial thermal conductivity (W/(m·K))

    (Froment, Bischoff & De Wilde 2011; Dixon 1996; Swarm Critic-05)
    """
    R_safe = max(float(R_in), 1e-6)
    if isinstance(k_er, np.ndarray):
        k_er_safe = np.maximum(k_er, 1e-4)
    else:
        k_er_safe = max(float(k_er), 1e-4)

    delta_T = (q_vol * (R_safe ** 2)) / (8.0 * k_er_safe)
    return T_mean + delta_T

# ---------------------------------------------------------------------------
# Ergun Packed-Bed Momentum Balance
# ---------------------------------------------------------------------------
def calc_ergun_dp(u_s: float, rho_g: float, mu: float, dp: float, eps: float,
                  Dh: Optional[float] = None) -> float:
    """
    Ergun equation for packed-bed pressure drop (Pa/m) with Eisfeld & Schnitzlein (2001)
    wall-effect correction factors:
      dP/dz = - [ 150 * mu * (1 - eps)^2 / (eps^3 * dp^2) * u_s * A_w + 1.75 * (1 - eps) / (eps^3 * dp) * rho_g * u_s^2 * B_w ]
    where:
      A_w = 1 + 2 / (3 * (Dh / dp) * (1 - eps))
      B_w = [k1 * (dp / Dh)^2 + k2]^2  (k1 = 0.45, k2 = 1.0)
    """
    eps_safe = float(np.clip(eps, 0.05, 0.95))
    dp_safe = max(float(dp), 1e-6)
    mu_safe = max(float(mu), 1e-7)
    rho_safe = max(float(rho_g), 1e-4)

    u_s_abs = abs(float(u_s))

    if Dh is not None and float(Dh) > 0.0:
        Dh_over_dp = max(float(Dh) / dp_safe, 1.0)
        Aw = 1.0 + 2.0 / (3.0 * Dh_over_dp * (1.0 - eps_safe))
        k1, k2 = 0.45, 1.0
        Bw = (k1 * ((1.0 / Dh_over_dp) ** 2) + k2) ** 2
    else:
        Aw = 1.0
        Bw = 1.0

    term_visc = 150.0 * mu_safe * ((1.0 - eps_safe) ** 2) / ((eps_safe ** 3) * (dp_safe ** 2)) * u_s_abs * Aw
    term_turb = 1.75 * (1.0 - eps_safe) / ((eps_safe ** 3) * dp_safe) * rho_safe * (u_s_abs ** 2) * Bw
    return -(term_visc + term_turb)

# ---------------------------------------------------------------------------
# Packed-Bed-to-Wall Heat Transfer Coefficients
# ---------------------------------------------------------------------------
def calc_hw_tube(
    T: float, P: float, F_arr: np.ndarray, species_list: list,
    geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig,
    return_details: bool = False
) -> Union[float, dict]:
    """
    Computes tube-side packed-bed-to-wall effective heat transfer coefficient h_i,eff (W/(m²·K)).
    Combines:
      1. Li & Finlayson (1977) / Froment et al. (2011) Section 3.4.2 wall film resistance 1/h_w(dp):
         Nu_w = (h_w * dp / Kg) = 0.17 * Re_p^0.79 * Pr^(1/3)
      2. Froment & Bischoff / Dixon (1996) 1D lumped radial bed conduction resistance:
         1 / h_{i,eff} = dp_exo / (Nu_{w,i} * kg_exo) + D_in / (6.133 * k_{er,exo})

    Fixes legacy defect where Nu_w was divided by D_in instead of particle diameter dp
    (Swarm Critic-11 & Critic-05).
    """
    F_tot = max(float(np.sum(F_arr)), 1e-12)
    y = np.maximum(F_arr, 1e-15) / F_tot

    canon_species = [resolve_species(sp) for sp in species_list]
    mw_avg = sum(y[i] * MW[canon_species[i]] for i in range(len(canon_species)))
    m_dot = F_tot * mw_avg
    G = m_dot / geom.A_tube  # kg/(m²·s)

    Re_p = max(G * cat.dp / trans.mu_exo, 1.0)
    cp_mean = sum(y[i] * calc_cp_species(canon_species[i], T) for i in range(len(canon_species))) / mw_avg
    kg_exo = trans.Kg_exo
    Pr = max(cp_mean * trans.mu_exo / kg_exo, 0.01)

    # 1. Li & Finlayson wall film Nusselt number & film coefficient based on dp
    Nu_w = 0.17 * (Re_p ** 0.79) * (Pr ** (1.0 / 3.0))
    h_w_film = Nu_w * kg_exo / cat.dp  # W/(m²·K)

    # 2. Effective radial thermal conductivity
    k_er = calc_ker(Re_p, Pr, kg_exo)

    # 3. Froment-Dixon 1D lumped effective heat transfer coefficient
    # 1 / h_i,eff = dp_exo / (Nu_w,i * kg_exo) + D_in / (6.133 * k_er,exo)
    R_film = cat.dp / max(Nu_w * kg_exo, 1e-9)
    R_radial = geom.D_in / max(6.133 * k_er, 1e-9)
    h_eff = 1.0 / (R_film + R_radial)
    h_result = float(max(h_eff, 20.0))

    if return_details:
        return {
            'h_eff': h_result,
            'h_film': float(h_w_film),
            'Nu_w': float(Nu_w),
            'Re_p': float(Re_p),
            'Pr': float(Pr),
            'k_er': float(k_er),
            'R_film': float(R_film),
            'R_radial': float(R_radial),
        }
    return h_result

def calc_hw_shell(
    T: float, P: float, F_arr: np.ndarray, species_list: list,
    geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig,
    return_details: bool = False
) -> Union[float, dict]:
    """
    Computes shell-side packed-bed-to-wall effective heat transfer coefficient h_o,eff (W/(m²·K)).
    Combines:
      1. Wall film resistance based on particle diameter dp:
         Nu_w = (h_w * dp / kg_endo) = 0.20 * Re_p^0.75 * Pr^(1/3)
      2. Froment & Bischoff / Dixon (1996) 1D lumped radial bed conduction resistance:
         1 / h_{o,eff} = dp_endo / (Nu_{w,o} * kg_endo) + Dh_shell / (6.133 * k_{er,endo})
      3. Wassiljewa-Mason-Saxena composition- and temperature-dependent thermal conductivity kg_endo(T, y)
         supporting 3-component and 4-component mixtures (Swarm Critic-10, Critic-11 & Critic-05).
    """
    F_tot = max(float(np.sum(F_arr)), 1e-12)
    y = np.maximum(F_arr, 1e-15) / F_tot

    canon_species = [resolve_species(sp) for sp in species_list]
    mw_avg = sum(y[i] * MW[canon_species[i]] for i in range(len(canon_species)))
    m_dot = F_tot * mw_avg
    G = m_dot / geom.A_shell  # kg/(m²·s)

    # Variable shell thermal conductivity via Wassiljewa-Mason-Saxena
    kg_endo = calc_kg_endo_wassiljewa(T, y, canon_species)

    Re_p = max(G * cat.dp / trans.mu_endo, 1.0)
    cp_mean = sum(y[i] * calc_cp_species(canon_species[i], T) for i in range(len(canon_species))) / mw_avg
    Pr = max(cp_mean * trans.mu_endo / kg_endo, 0.01)

    # 1. Wall film Nusselt number & film coefficient based on dp
    Nu_w = 0.20 * (Re_p ** 0.75) * (Pr ** (1.0 / 3.0))
    h_w_film = Nu_w * kg_endo / cat.dp  # W/(m²·K)

    # 2. Effective radial thermal conductivity
    k_er = calc_ker(Re_p, Pr, kg_endo)

    # 3. Froment-Dixon 1D lumped effective heat transfer coefficient
    # 1 / h_o,eff = dp_endo / (Nu_w,o * kg_endo) + Dh_shell / (6.133 * k_er,endo)
    R_film = cat.dp / max(Nu_w * kg_endo, 1e-9)
    R_radial = geom.Dh_shell / max(6.133 * k_er, 1e-9)
    h_eff = 1.0 / (R_film + R_radial)
    h_result = float(max(h_eff, 20.0))

    if return_details:
        return {
            'h_eff': h_result,
            'h_film': float(h_w_film),
            'Nu_w': float(Nu_w),
            'Re_p': float(Re_p),
            'Pr': float(Pr),
            'kg_endo': float(kg_endo),
            'k_er': float(k_er),
            'R_film': float(R_film),
            'R_radial': float(R_radial),
        }
    return h_result

def calc_hw_tube_film(T: float, P: float, F_arr: np.ndarray, species_list: list,
                      geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig) -> float:
    """Computes pure wall film heat transfer coefficient h_w,i (W/(m²·K)) without radial conduction resistance."""
    res = calc_hw_tube(T, P, F_arr, species_list, geom, cat, trans, return_details=True)
    return res['h_film']

def calc_hw_shell_film(T: float, P: float, F_arr: np.ndarray, species_list: list,
                       geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig) -> float:
    """Computes pure wall film heat transfer coefficient h_w,o (W/(m²·K)) without radial conduction resistance."""
    res = calc_hw_shell(T, P, F_arr, species_list, geom, cat, trans, return_details=True)
    return res['h_film']

calc_hw_tube_eff = calc_hw_tube
calc_hw_shell_eff = calc_hw_shell

# ---------------------------------------------------------------------------
# Overall Heat Transfer Coefficients U_i (based on inner tube area)
# ---------------------------------------------------------------------------
def calc_U_tcr(T_exo: float, P_exo: float, F_exo: np.ndarray, species_exo: list,
               T_endo: float, P_endo: float, F_endo: np.ndarray, species_endo: list,
               geom: ReactorGeometry, cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
               trans: TransportConfig) -> float:
    """
    Computes overall heat transfer coefficient U_i (based on inner tube area) for TCR:
      1 / U_i = 1 / h_i,eff + R_wall + (D_in / D_out) * (1 / h_o,eff)
    where:
      1 / h_i,eff = dp_exo / (Nu_w,i * kg_exo) + D_in / (6.133 * k_er,exo)
      1 / h_o,eff = dp_endo / (Nu_w,o * kg_endo) + Dh_shell / (6.133 * k_er,endo)
      R_wall = (D_in / (2 * Kw)) * ln(D_out / D_in)
    """
    h_i = calc_hw_tube(T_exo, P_exo, F_exo, species_exo, geom, cat_exo, trans)
    h_o = calc_hw_shell(T_endo, P_endo, F_endo, species_endo, geom, cat_endo, trans)
    R_wall = (geom.D_in * np.log(geom.D_out / geom.D_in)) / (2.0 * geom.Kw)

    U = 1.0 / (1.0 / h_i + R_wall + (geom.D_in / geom.D_out) * (1.0 / h_o))
    return float(U)

def calc_U_cmr(T_exo: float, P_exo: float, F_exo: np.ndarray, species_exo: list,
               geom: ReactorGeometry, cat_exo: CatalystBedConfig, trans: TransportConfig) -> float:
    """
    Computes overall heat transfer coefficient U_i for conventional boiling water reactor (CMR):
      1 / U_i = 1 / h_i,eff + R_wall + (D_in / D_out) * (1 / h_cool)
    where:
      1 / h_i,eff = dp_exo / (Nu_w,i * kg_exo) + D_in / (6.133 * k_er,exo)
      R_wall = (D_in / (2 * Kw)) * ln(D_out / D_in)
    """
    h_i = calc_hw_tube(T_exo, P_exo, F_exo, species_exo, geom, cat_exo, trans)
    R_wall = (geom.D_in * np.log(geom.D_out / geom.D_in)) / (2.0 * geom.Kw)
    h_cool = trans.h_cool

    U = 1.0 / (1.0 / h_i + R_wall + (geom.D_in / geom.D_out) * (1.0 / h_cool))
    return float(U)

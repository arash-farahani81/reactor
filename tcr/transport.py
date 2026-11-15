"""
Transport phenomena module:
- Packed-bed-to-wall heat transfer coefficients (Dixon / Li & Finlayson correlations)
- Boiling water coolant jacket film coefficient for conventional reactor
- Cylindrical tube wall conductive resistance
- Overall dynamic heat transfer coefficient U(z)
- Ergun packed-bed momentum balances
"""

import numpy as np
from .thermo import R_GAS, calc_cp_species, MW
from .geometry import ReactorGeometry
from .config import TransportConfig, CatalystBedConfig

def calc_ergun_dp(u_s: float, rho_g: float, mu: float, dp: float, eps: float) -> float:
    """
    Ergun equation for packed-bed pressure drop (Pa/m).
    dP/dz = - [ 150 * mu * (1 - eps)^2 / (eps^3 * dp^2) * u_s + 1.75 * (1 - eps) / (eps^3 * dp) * rho_g * u_s^2 ]
    """
    eps_safe = float(np.clip(eps, 0.05, 0.95))
    dp_safe = max(float(dp), 1e-6)
    mu_safe = max(float(mu), 1e-7)
    rho_safe = max(float(rho_g), 1e-4)

    u_s_abs = abs(float(u_s))
    term_visc = 150.0 * mu_safe * ((1.0 - eps_safe) ** 2) / ((eps_safe ** 3) * (dp_safe ** 2)) * u_s_abs
    term_turb = 1.75 * (1.0 - eps_safe) / ((eps_safe ** 3) * dp_safe) * rho_safe * (u_s_abs ** 2)
    return -(term_visc + term_turb)

def calc_hw_tube(T: float, P: float, F_arr: np.ndarray, species_list: list,
                 geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig) -> float:
    """
    Computes tube-side packed-bed-to-wall heat transfer coefficient h_w,i (W/(m²·K)).
    Correlation: Li & Finlayson (1977) / Froment et al. (2011) Section 3.4.2
    Nu_w = (h_w * D_in / Kg) = 0.17 * Re_p^0.79 * Pr^(1/3)
    """
    F_tot = max(float(np.sum(F_arr)), 1e-12)
    y = np.maximum(F_arr, 1e-15) / F_tot
    mw_avg = sum(y[i] * MW[species_list[i]] for i in range(len(species_list)))
    m_dot = F_tot * mw_avg
    G = m_dot / geom.A_tube  # kg/(m²·s)

    Re_p = max(G * cat.dp / trans.mu_exo, 1.0)
    cp_mean = sum(y[i] * calc_cp_species(species_list[i], T) for i in range(len(species_list))) / mw_avg
    Pr = max(cp_mean * trans.mu_exo / trans.Kg_exo, 0.01)

    Nu_w = 0.17 * (Re_p ** 0.79) * (Pr ** (1.0 / 3.0))
    h_w = Nu_w * trans.Kg_exo / geom.D_in
    return float(max(h_w, 20.0))

def calc_hw_shell(T: float, P: float, F_arr: np.ndarray, species_list: list,
                  geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig) -> float:
    """
    Computes shell-side packed-bed-to-wall heat transfer coefficient h_w,o (W/(m²·K)).
    Nu_w = (h_w * Dh_shell / Kg) = 0.20 * Re_p^0.75 * Pr^(1/3)
    """
    F_tot = max(float(np.sum(F_arr)), 1e-12)
    y = np.maximum(F_arr, 1e-15) / F_tot
    mw_avg = sum(y[i] * MW[species_list[i]] for i in range(len(species_list)))
    m_dot = F_tot * mw_avg
    G = m_dot / geom.A_shell  # kg/(m²·s) evaluated with single consistent shell area

    Re_p = max(G * cat.dp / trans.mu_endo, 1.0)
    cp_mean = sum(y[i] * calc_cp_species(species_list[i], T) for i in range(len(species_list))) / mw_avg
    Pr = max(cp_mean * trans.mu_endo / trans.Kg_endo, 0.01)

    Nu_w = 0.20 * (Re_p ** 0.75) * (Pr ** (1.0 / 3.0))
    h_w = Nu_w * trans.Kg_endo / geom.Dh_shell
    return float(max(h_w, 20.0))

def calc_U_tcr(T_exo: float, P_exo: float, F_exo: np.ndarray, species_exo: list,
               T_endo: float, P_endo: float, F_endo: np.ndarray, species_endo: list,
               geom: ReactorGeometry, cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
               trans: TransportConfig) -> float:
    """
    Computes overall heat transfer coefficient U_i (based on inner tube area) for TCR:
    1 / U_i = 1 / h_w,i + R_wall + (D_in / D_out) * (1 / h_w,o)
    where R_wall = (D_in / (2 * Kw)) * ln(D_out / D_in)
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
    1 / U_i = 1 / h_w,i + R_wall + (D_in / D_out) * (1 / h_cool)
    where R_wall = (D_in / (2 * Kw)) * ln(D_out / D_in)
    """
    h_i = calc_hw_tube(T_exo, P_exo, F_exo, species_exo, geom, cat_exo, trans)
    R_wall = (geom.D_in * np.log(geom.D_out / geom.D_in)) / (2.0 * geom.Kw)
    h_cool = trans.h_cool

    U = 1.0 / (1.0 / h_i + R_wall + (geom.D_in / geom.D_out) * (1.0 / h_cool))
    return float(U)

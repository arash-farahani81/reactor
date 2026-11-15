"""
Thermodynamic functions, temperature-dependent heat capacities,
Kirchhoff reaction enthalpies, and equilibrium relations.
All coefficients sourced from NIST Chemistry WebBook & Reid, Prausnitz, Poling.
Single Source of Truth (SSOT) for all enthalpy, heat capacity, and equilibrium calculations.
"""

import numpy as np

R_GAS = 8.314462618  # J/(mol·K), CODATA 2018
T_REF_ENTHALPY = 298.15  # K, standard enthalpy reference temperature
T_REF_VANTHOFF = 500.0   # K, central reference temperature for van 't Hoff anchoring

# Molecular Weights (kg/mol)
MW = {
    # Tube components
    'H2': 2.01588e-3,
    'CO2': 44.0095e-3,
    'CO': 28.0101e-3,
    'CH3OH': 32.04186e-3,
    'H2O': 18.01528e-3,
    'CH4': 16.0425e-3,
    'N2': 28.0134e-3,
    # Shell components
    'C6H12': 84.1595e-3,
    'C6H6': 78.1118e-3,
    'Ar': 39.948e-3,
}

# Standard Molar Enthalpies of Formation at 298.15 K (J/mol) (NIST WebBook)
DELTA_H_F_298 = {
    'H2': 0.0,
    'CO2': -393510.0,
    'CO': -110530.0,
    'CH3OH': -200670.0,  # gas
    'H2O': -241820.0,    # gas
    'CH4': -74870.0,
    'N2': 0.0,
    'C6H12': -123140.0,  # gas
    'C6H6': 82930.0,     # gas
    'Ar': 0.0,
}

# Heat capacity polynomial coefficients: Cp(T) = a0 + a1*T + a2*T^2 + a3*T^3 (J/(mol·K))
# Source: Reid, Prausnitz, Poling (Properties of Gases and Liquids) / NIST WebBook (300 - 800 K)
CP_COEFFS = {
    'H2':     np.array([28.84,    7.65e-5,   3.288e-6, -8.698e-10]),
    'CO2':    np.array([19.80,    0.07344,  -5.602e-5,  1.715e-8]),
    'CO':     np.array([28.16,    1.675e-3,  5.372e-6, -2.222e-9]),
    'CH3OH':  np.array([19.00,    0.09152,  -1.220e-5, -8.039e-9]),
    'H2O':    np.array([32.24,    1.924e-3,  1.055e-5, -3.596e-9]),
    'CH4':    np.array([19.25,    0.05213,   1.197e-5, -1.132e-8]),
    'N2':     np.array([28.88,   -1.570e-3,  8.080e-6, -2.870e-9]),
    'C6H12':  np.array([-50.40,   0.6470,   -3.610e-4,  7.820e-8]),
    'C6H6':   np.array([-33.90,   0.4720,   -2.980e-4,  7.080e-8]),
    'Ar':     np.array([20.786,   0.0,       0.0,       0.0]),
}

def calc_cp_species(species: str, T: float) -> float:
    """Computes molar heat capacity of a pure species at temperature T (K) in J/(mol·K)."""
    c = CP_COEFFS[species]
    return float(c[0] + c[1] * T + c[2] * (T ** 2) + c[3] * (T ** 3))

def calc_enthalpy_species(species: str, T: float) -> float:
    """
    Computes absolute molar enthalpy H_i(T) = Delta H_f,298 + integral(298.15 to T) Cp(T') dT' in J/mol.
    """
    c = CP_COEFFS[species]
    T_ref = T_REF_ENTHALPY
    integral_cp = (
        c[0] * (T - T_ref)
        + 0.5 * c[1] * (T ** 2 - T_ref ** 2)
        + (1.0 / 3.0) * c[2] * (T ** 3 - T_ref ** 3)
        + 0.25 * c[3] * (T ** 4 - T_ref ** 4)
    )
    return float(DELTA_H_F_298[species] + integral_cp)

def calc_dH_rxn_meoh(T: float) -> float:
    """
    Enthalpy of reaction for CO2 + 3 H2 <-> CH3OH + H2O (exothermic) in J/mol.
    Delta H_r(T) = H_CH3OH + H_H2O - H_CO2 - 3*H_H2
    """
    return float(
        calc_enthalpy_species('CH3OH', T)
        + calc_enthalpy_species('H2O', T)
        - calc_enthalpy_species('CO2', T)
        - 3.0 * calc_enthalpy_species('H2', T)
    )

def calc_dH_rxn_rwgs(T: float) -> float:
    """
    Enthalpy of reaction for CO2 + H2 <-> CO + H2O (endothermic) in J/mol.
    Delta H_r(T) = H_CO + H_H2O - H_CO2 - H_H2
    """
    return float(
        calc_enthalpy_species('CO', T)
        + calc_enthalpy_species('H2O', T)
        - calc_enthalpy_species('CO2', T)
        - calc_enthalpy_species('H2', T)
    )

def calc_dH_rxn_ch(T: float) -> float:
    """
    Enthalpy of reaction for C6H12 <-> C6H6 + 3 H2 (endothermic) in J/mol.
    Delta H_r(T) = H_C6H6 + 3*H_H2 - H_C6H12
    """
    return float(
        calc_enthalpy_species('C6H6', T)
        + 3.0 * calc_enthalpy_species('H2', T)
        - calc_enthalpy_species('C6H12', T)
    )

# ---------------------------------------------------------------------------
# Exact van 't Hoff Thermodynamic Integration Framework
# ---------------------------------------------------------------------------
# Strict compliance with van 't Hoff equation: d(ln K)/dT = Delta H_r(T) / (R * T^2).
# Delta H_r(T) is identical to the Kirchhoff enthalpy functions defined above:
#   Delta H_r(T) = A0 + Delta a0*T + (Delta a1/2)*T^2 + (Delta a2/3)*T^3 + (Delta a3/4)*T^4
# Analytical integral of Delta H_r(T) / (R * T^2) dT:
#   I(T) = (1/R) * [-A0/T + Delta a0*ln(T) + (Delta a1/2)*T + (Delta a2/6)*T^2 + (Delta a3/12)*T^3]
#   ln K(T) = I(T) + C
# By the Fundamental Theorem of Calculus: d(ln K)/dT = Delta H_r(T) / (R * T^2) holds to machine precision.

def _build_van_t_hoff_params(stoich: dict[str, float], ln_K_ref: float, T_0: float = 500.0) -> tuple[float, np.ndarray, float]:
    """Computes exact analytical integration parameters A0, Delta Cp, and C for a reaction."""
    dH_298 = sum(DELTA_H_F_298[sp] * coeff for sp, coeff in stoich.items())
    dCp = sum(CP_COEFFS[sp] * coeff for sp, coeff in stoich.items())
    T_r = T_REF_ENTHALPY
    A0 = dH_298 - (
        dCp[0] * T_r
        + 0.5 * dCp[1] * (T_r ** 2)
        + (1.0 / 3.0) * dCp[2] * (T_r ** 3)
        + 0.25 * dCp[3] * (T_r ** 4)
    )
    I_T0 = (1.0 / R_GAS) * (
        -A0 / T_0
        + dCp[0] * np.log(T_0)
        + 0.5 * dCp[1] * T_0
        + (1.0 / 6.0) * dCp[2] * (T_0 ** 2)
        + (1.0 / 12.0) * dCp[3] * (T_0 ** 3)
    )
    C_const = ln_K_ref - I_T0
    return float(A0), dCp, float(C_const)

# Literature anchoring values at T_0 = 500 K:
# Reaction 1 (MeOH): log10(Keq1) = 3066/500 - 10.592 = -4.460 (Graaf 1988 / VBF 1996)
_LN_KEQ1_REF_500 = np.log(10.0) * (3066.0 / 500.0 - 10.592)
# Reaction 2 (RWGS): log10(Keq2) = -2073/500 + 2.029 = -2.117 (Graaf 1988 / VBF 1996)
_LN_KEQ2_REF_500 = np.log(10.0) * (-2073.0 / 500.0 + 2.029)
# Reaction 3 (Dehydrogenation): ln(Kp_ch) = -205400 / (R * 500) + 43.30 (Itoh 1987 / Khademi 2009)
_LN_KP_CH_REF_500 = -205400.0 / (R_GAS * 500.0) + 43.30

_VTH_PARAMS_MEOH = _build_van_t_hoff_params({'CH3OH': 1.0, 'H2O': 1.0, 'CO2': -1.0, 'H2': -3.0}, _LN_KEQ1_REF_500)
_VTH_PARAMS_RWGS = _build_van_t_hoff_params({'CO': 1.0, 'H2O': 1.0, 'CO2': -1.0, 'H2': -1.0}, _LN_KEQ2_REF_500)
_VTH_PARAMS_CH   = _build_van_t_hoff_params({'C6H6': 1.0, 'H2': 3.0, 'C6H12': -1.0}, _LN_KP_CH_REF_500)

def _eval_van_t_hoff(params: tuple[float, np.ndarray, float], T: float) -> float:
    """Evaluates equilibrium constant K(T) from exact analytical van 't Hoff integral."""
    A0, dCp, C_const = params
    T_safe = max(float(T), 100.0)
    ln_K = (1.0 / R_GAS) * (
        -A0 / T_safe
        + dCp[0] * np.log(T_safe)
        + 0.5 * dCp[1] * T_safe
        + (1.0 / 6.0) * dCp[2] * (T_safe ** 2)
        + (1.0 / 12.0) * dCp[3] * (T_safe ** 3)
    ) + C_const
    return float(np.exp(np.clip(ln_K, -50.0, 50.0)))

def calc_keq_meoh(T: float) -> float:
    """
    Thermodynamically consistent equilibrium constant Keq1 (bar^-2) for methanol synthesis:
    CO2 + 3 H2 <-> CH3OH + H2O (Delta n = -2).
    Strictly obeys d(ln Keq1)/dT = Delta H_r,meoh(T) / (R * T^2) across 400-600 K.
    """
    return _eval_van_t_hoff(_VTH_PARAMS_MEOH, T)

def calc_keq_rwgs(T: float) -> float:
    """
    Thermodynamically consistent equilibrium constant Keq2 (dimensionless) for reverse water-gas shift:
    CO2 + H2 <-> CO + H2O (Delta n = 0).
    Strictly obeys d(ln Keq2)/dT = Delta H_r,rwgs(T) / (R * T^2) across 400-600 K.
    """
    return _eval_van_t_hoff(_VTH_PARAMS_RWGS, T)

def calc_kp_dehydrogenation(T: float) -> float:
    """
    Thermodynamically consistent equilibrium constant Kp (bar^3) for cyclohexane dehydrogenation:
    C6H12 <-> C6H6 + 3 H2 (Delta n = +3).
    Strictly obeys d(ln Kp)/dT = Delta H_r,ch(T) / (R * T^2) across 400-600 K.
    """
    return _eval_van_t_hoff(_VTH_PARAMS_CH, T)

# Standardized Aliases
calc_keq1 = calc_keq_meoh
calc_keq2 = calc_keq_rwgs
calc_kp_ch = calc_kp_dehydrogenation
Keq1 = calc_keq_meoh
Keq2 = calc_keq_rwgs
Kp_CH = calc_kp_dehydrogenation

# Exact Analytical Derivative Functions for Verification
def calc_dln_keq1_dT(T: float) -> float:
    """Exact analytical derivative d(ln Keq1)/dT = Delta H_r,meoh(T) / (R * T^2) in 1/K."""
    return float(calc_dH_rxn_meoh(T) / (R_GAS * (T ** 2)))

def calc_dln_keq2_dT(T: float) -> float:
    """Exact analytical derivative d(ln Keq2)/dT = Delta H_r,rwgs(T) / (R * T^2) in 1/K."""
    return float(calc_dH_rxn_rwgs(T) / (R_GAS * (T ** 2)))

def calc_dln_kp_ch_dT(T: float) -> float:
    """Exact analytical derivative d(ln Kp_ch)/dT = Delta H_r,ch(T) / (R * T^2) in 1/K."""
    return float(calc_dH_rxn_ch(T) / (R_GAS * (T ** 2)))

def calc_dln_K_dT(reaction: str, T: float) -> float:
    """Returns exact analytical d(ln K)/dT for 'meoh', 'rwgs', or 'ch'."""
    rxn = reaction.lower()
    if rxn in ('meoh', '1', 'keq1'):
        return calc_dln_keq1_dT(T)
    elif rxn in ('rwgs', '2', 'keq2'):
        return calc_dln_keq2_dT(T)
    elif rxn in ('ch', '3', 'kp_ch', 'dehydro'):
        return calc_dln_kp_ch_dT(T)
    else:
        raise ValueError(f"Unknown reaction {reaction}")

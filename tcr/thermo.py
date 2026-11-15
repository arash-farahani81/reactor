"""
Thermodynamic functions, temperature-dependent heat capacities,
Kirchhoff reaction enthalpies, and equilibrium relations.
All coefficients sourced from NIST Chemistry WebBook & Reid, Prausnitz, Poling.
"""

import numpy as np

R_GAS = 8.314462618  # J/(mol·K), CODATA 2018

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
    T_ref = 298.15
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

def calc_kp_dehydrogenation(T: float) -> float:
    """
    Thermodynamically consistent equilibrium constant Kp (bar^3) for cyclohexane dehydrogenation:
    C6H12 <-> C6H6 + 3 H2 (Delta n = +3).
    ln(Kp) = -Delta H / (R * T) + Delta S / R = -205400 / (R * T) + 43.30.
    """
    ln_kp = -205400.0 / (R_GAS * T) + 43.30
    return float(np.exp(np.clip(ln_kp, -30.0, 30.0)))

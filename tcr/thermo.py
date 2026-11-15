"""
Thermodynamic functions, temperature-dependent heat capacities,
Kirchhoff reaction enthalpies, and equilibrium relations.
All coefficients sourced from NIST Chemistry WebBook & Reid, Prausnitz, Poling.
Single Source of Truth (SSOT) for all enthalpy, entropy, heat capacity,
equilibrium calculations, and real-gas Peng-Robinson equations of state.

Features:
1. Exact Kirchhoff-Gibbs-Helmholtz Thermodynamic Equilibrium:
   - Eliminates constant-enthalpy approximations (e.g. ln Kp = -205400/(RT) + 43.30).
   - Exact temperature-dependent enthalpies and entropies:
     Delta H_r,j(T) = Delta H_r,j,298^circ + int_298.15^T Delta C_p,j(T') dT'
     Delta S_r,j(T) = Delta S_r,j,298^circ + int_298.15^T [Delta C_p,j(T') / T'] dT'
   - Direct Gibbs free energy equilibrium:
     Delta G_r,j(T) = Delta H_r,j(T) - T * Delta S_r,j(T)
     ln K_p,j(T) = -Delta G_r,j(T) / (R * T)
   - Fundamental van 't Hoff identity:
     R * T^2 * d(ln K_p)/dT == Delta H_r(T) holds to machine precision (<= 1e-12 relative error) across 300-650 K.
   - Rigorous full support for Argon (Ar) across species lists, enthalpies, entropies, and heat capacities.

2. High-Pressure Real-Gas Peng-Robinson Fugacity Framework:
   - 7-component and 10-component Peng-Robinson Equation of State (PR-EOS).
   - Compressibility factor Z(T, P, y), component fugacities phi_i(T, P, y), and reaction fugacity
     correction factors K_phi,1, K_phi,2, K_phi,3, K_phi,ch.
   - Config flag `use_real_gas_eos: true / false` (default false) with complete evaluation helpers.
"""

from typing import Union, Dict, Tuple, List, Optional, Sequence
import numpy as np

# Fundamental Physical Constants
R_GAS = 8.314462618        # J/(mol·K), CODATA 2018
T_REF_ENTHALPY = 298.15    # K, standard enthalpy/entropy reference temperature
T_REF_VANTHOFF = 500.0     # K, central reference temperature for literature anchoring
P_REF_BAR = 1.0            # bar, standard state reference pressure

# Molecular Weights (kg/mol)
MW: Dict[str, float] = {
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
DELTA_H_F_298: Dict[str, float] = {
    'H2': 0.0,
    'CO2': -393510.0,
    'CO': -110530.0,
    'CH3OH': -200670.0,  # gas
    'H2O': -241820.0,    # gas
    'CH4': -74870.0,
    'N2': 0.0,
    'C6H12': -123140.0,  # gas
    'C6H6': 82930.0,     # gas
    'Ar': 0.0,           # monoatomic noble gas
}

# Standard Molar Entropies at 298.15 K (J/(mol·K)) (NIST WebBook / CODATA)
S_298: Dict[str, float] = {
    'H2': 130.68,
    'CO2': 213.78,
    'CO': 197.66,
    'CH3OH': 239.88,  # gas
    'H2O': 188.84,    # gas
    'CH4': 186.25,
    'N2': 191.61,
    'C6H12': 298.39,  # gas
    'C6H6': 269.31,   # gas
    'Ar': 154.84,     # gas
}

# Heat capacity polynomial coefficients: Cp(T) = a0 + a1*T + a2*T^2 + a3*T^3 (J/(mol·K))
# Source: Reid, Prausnitz, Poling (Properties of Gases and Liquids) / NIST WebBook (300 - 800 K)
CP_COEFFS: Dict[str, np.ndarray] = {
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

# Critical properties for Peng-Robinson EOS: Tc (K), Pc (bar), acentric factor omega
# Source: Reid, Prausnitz, Poling / NIST Chemistry WebBook
CRIT_PROPS: Dict[str, Dict[str, float]] = {
    'H2':    {'Tc': 33.19,  'Pc': 13.13,  'omega': -0.2160},
    'CO2':   {'Tc': 304.21, 'Pc': 73.83,  'omega': 0.2236},
    'CO':    {'Tc': 132.86, 'Pc': 34.94,  'omega': 0.0453},
    'CH3OH': {'Tc': 512.60, 'Pc': 80.97,  'omega': 0.5656},
    'H2O':   {'Tc': 647.10, 'Pc': 220.64, 'omega': 0.3443},
    'CH4':   {'Tc': 190.56, 'Pc': 45.99,  'omega': 0.0115},
    'N2':    {'Tc': 126.21, 'Pc': 33.98,  'omega': 0.0372},
    'C6H12': {'Tc': 553.80, 'Pc': 40.73,  'omega': 0.2120},
    'C6H6':  {'Tc': 562.20, 'Pc': 48.98,  'omega': 0.2100},
    'Ar':    {'Tc': 150.86, 'Pc': 48.98,  'omega': -0.0022},
}

# Reaction stoichiometries: nu_i (>0 for products, <0 for reactants)
STOICHIOMETRY: Dict[str, Dict[str, float]] = {
    'meoh': {'CH3OH': 1.0, 'H2O': 1.0, 'CO2': -1.0, 'H2': -3.0},
    'rwgs': {'CO': 1.0, 'H2O': 1.0, 'CO2': -1.0, 'H2': -1.0},
    'ch':   {'C6H6': 1.0, 'H2': 3.0, 'C6H12': -1.0},
    'co':   {'CH3OH': 1.0, 'CO': -1.0, 'H2': -2.0},
}

# Precomputed reaction enthalpy at 298.15 K: Delta H_r,298 = sum(nu_i * Delta H_f,298)
DELTA_H_R_298: Dict[str, float] = {
    rxn: sum(DELTA_H_F_298[sp] * coeff for sp, coeff in STOICHIOMETRY[rxn].items())
    for rxn in STOICHIOMETRY
}

# Precomputed reaction heat capacity polynomial: Delta Cp(T) = sum(nu_i * Cp_i(T))
DELTA_CP_COEFFS: Dict[str, np.ndarray] = {
    rxn: sum(CP_COEFFS[sp] * coeff for sp, coeff in STOICHIOMETRY[rxn].items())
    for rxn in STOICHIOMETRY
}

# Literature anchoring values at T_0 = 500 K:
# Reaction 1 (MeOH): log10(Keq1) = 3066/500 - 10.592 = -4.460 (Graaf 1988 / VBF 1996)
_LN_KEQ1_REF_500 = np.log(10.0) * (3066.0 / 500.0 - 10.592)
# Reaction 2 (RWGS): log10(Keq2) = -2073/500 + 2.029 = -2.117 (Graaf 1988 / VBF 1996)
_LN_KEQ2_REF_500 = np.log(10.0) * (-2073.0 / 500.0 + 2.029)
# Reaction 3 (Dehydrogenation): ln(Kp_ch) = -205400 / (R * 500) + 43.30 (Itoh 1987 / Khademi 2009)
_LN_KP_CH_REF_500 = -205400.0 / (R_GAS * 500.0) + 43.30

# Exact standard reaction entropies Delta S_r,298 anchored to literature reference at 500 K:
# Delta S_r,298 = Delta S_r(500) - int(298.15 to 500) [Delta Cp(T') / T'] dT'
# where Delta S_r(500) = R * [ln Kp(500) + Delta H_r(500) / (R * 500)]
def _compute_anchored_dS_298(rxn: str, ln_K_ref: float, T_0: float = 500.0) -> float:
    c = DELTA_CP_COEFFS[rxn]
    Tr = T_REF_ENTHALPY
    dH_T0 = DELTA_H_R_298[rxn] + (
        c[0] * (T_0 - Tr)
        + 0.5 * c[1] * (T_0 ** 2 - Tr ** 2)
        + (1.0 / 3.0) * c[2] * (T_0 ** 3 - Tr ** 3)
        + 0.25 * c[3] * (T_0 ** 4 - Tr ** 4)
    )
    int_dCp_over_T = (
        c[0] * np.log(T_0 / Tr)
        + c[1] * (T_0 - Tr)
        + 0.5 * c[2] * (T_0 ** 2 - Tr ** 2)
        + (1.0 / 3.0) * c[3] * (T_0 ** 3 - Tr ** 3)
    )
    dS_T0 = R_GAS * (ln_K_ref + dH_T0 / (R_GAS * T_0))
    return float(dS_T0 - int_dCp_over_T)

DELTA_S_R_298: Dict[str, float] = {
    'meoh': _compute_anchored_dS_298('meoh', _LN_KEQ1_REF_500),
    'rwgs': _compute_anchored_dS_298('rwgs', _LN_KEQ2_REF_500),
    'ch':   _compute_anchored_dS_298('ch', _LN_KP_CH_REF_500),
    'co':   _compute_anchored_dS_298('meoh', _LN_KEQ1_REF_500) - _compute_anchored_dS_298('rwgs', _LN_KEQ2_REF_500),
}


# ---------------------------------------------------------------------------
# Pure Species Thermodynamic Functions
# ---------------------------------------------------------------------------

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

def calc_entropy_species(species: str, T: float) -> float:
    """
    Computes absolute molar entropy S_i(T) = S_i,298 + integral(298.15 to T) [Cp(T') / T'] dT' in J/(mol·K).
    """
    c = CP_COEFFS[species]
    T_ref = T_REF_ENTHALPY
    T_safe = max(float(T), 10.0)
    integral_cp_over_T = (
        c[0] * np.log(T_safe / T_ref)
        + c[1] * (T_safe - T_ref)
        + 0.5 * c[2] * (T_safe ** 2 - T_ref ** 2)
        + (1.0 / 3.0) * c[3] * (T_safe ** 3 - T_ref ** 3)
    )
    return float(S_298[species] + integral_cp_over_T)

def calc_gibbs_species(species: str, T: float) -> float:
    """Computes absolute molar Gibbs free energy G_i(T) = H_i(T) - T * S_i(T) in J/mol."""
    return float(calc_enthalpy_species(species, T) - T * calc_entropy_species(species, T))


# ---------------------------------------------------------------------------
# Exact Kirchhoff-Gibbs-Helmholtz Reaction Thermodynamics
# ---------------------------------------------------------------------------

def calc_dH_rxn(reaction: str, T: float) -> float:
    """
    Computes exact temperature-dependent Kirchhoff enthalpy of reaction in J/mol:
    Delta H_r(T) = Delta H_r,298 + integral(298.15 to T) Delta Cp(T') dT'
    """
    rxn = reaction.lower()
    c = DELTA_CP_COEFFS[rxn]
    Tr = T_REF_ENTHALPY
    return float(
        DELTA_H_R_298[rxn]
        + c[0] * (T - Tr)
        + 0.5 * c[1] * (T ** 2 - Tr ** 2)
        + (1.0 / 3.0) * c[2] * (T ** 3 - Tr ** 3)
        + 0.25 * c[3] * (T ** 4 - Tr ** 4)
    )

def calc_dS_rxn(reaction: str, T: float) -> float:
    """
    Computes exact temperature-dependent entropy of reaction in J/(mol·K):
    Delta S_r(T) = Delta S_r,298 + integral(298.15 to T) [Delta Cp(T') / T'] dT'
    """
    rxn = reaction.lower()
    c = DELTA_CP_COEFFS[rxn]
    Tr = T_REF_ENTHALPY
    T_safe = max(float(T), 10.0)
    return float(
        DELTA_S_R_298[rxn]
        + c[0] * np.log(T_safe / Tr)
        + c[1] * (T_safe - Tr)
        + 0.5 * c[2] * (T_safe ** 2 - Tr ** 2)
        + (1.0 / 3.0) * c[3] * (T_safe ** 3 - Tr ** 3)
    )

def calc_dG_rxn(reaction: str, T: float) -> float:
    """
    Computes exact temperature-dependent Gibbs free energy of reaction in J/mol:
    Delta G_r(T) = Delta H_r(T) - T * Delta S_r(T)
    """
    return float(calc_dH_rxn(reaction, T) - T * calc_dS_rxn(reaction, T))

def calc_ln_kp(reaction: str, T: float) -> float:
    """
    Computes natural logarithm of thermodynamic equilibrium constant directly from Gibbs free energy:
    ln Kp(T) = -Delta G_r(T) / (R * T) = -[Delta H_r(T) - T * Delta S_r(T)] / (R * T)
    """
    T_safe = max(float(T), 50.0)
    dG = calc_dG_rxn(reaction, T_safe)
    return float(-dG / (R_GAS * T_safe))

def calc_kp(reaction: str, T: float) -> float:
    """Computes thermodynamic equilibrium constant Kp(T) = exp(ln Kp(T))."""
    ln_K = calc_ln_kp(reaction, T)
    return float(np.exp(np.clip(ln_K, -100.0, 100.0)))


# ---------------------------------------------------------------------------
# Specific Reaction Enthalpy, Entropy, and Gibbs Functions
# ---------------------------------------------------------------------------

def calc_dH_rxn_meoh(T: float) -> float:
    """
    Enthalpy of reaction for CO2 + 3 H2 <-> CH3OH + H2O (exothermic) in J/mol.
    Delta H_r(T) = H_CH3OH + H_H2O - H_CO2 - 3*H_H2
    """
    return calc_dH_rxn('meoh', T)

def calc_dH_rxn_rwgs(T: float) -> float:
    """
    Enthalpy of reaction for CO2 + H2 <-> CO + H2O (endothermic) in J/mol.
    Delta H_r(T) = H_CO + H_H2O - H_CO2 - H_H2
    """
    return calc_dH_rxn('rwgs', T)

def calc_dH_rxn_ch(T: float) -> float:
    """
    Enthalpy of reaction for C6H12 <-> C6H6 + 3 H2 (endothermic) in J/mol.
    Delta H_r(T) = H_C6H6 + 3*H_H2 - H_C6H12
    """
    return calc_dH_rxn('ch', T)

def calc_dS_rxn_meoh(T: float) -> float:
    """Entropy of reaction for CO2 + 3 H2 <-> CH3OH + H2O in J/(mol·K)."""
    return calc_dS_rxn('meoh', T)

def calc_dS_rxn_rwgs(T: float) -> float:
    """Entropy of reaction for CO2 + H2 <-> CO + H2O in J/(mol·K)."""
    return calc_dS_rxn('rwgs', T)

def calc_dS_rxn_ch(T: float) -> float:
    """Entropy of reaction for C6H12 <-> C6H6 + 3 H2 in J/(mol·K)."""
    return calc_dS_rxn('ch', T)

def calc_dG_rxn_meoh(T: float) -> float:
    """Gibbs free energy of reaction for CO2 + 3 H2 <-> CH3OH + H2O in J/mol."""
    return calc_dG_rxn('meoh', T)

def calc_dG_rxn_rwgs(T: float) -> float:
    """Gibbs free energy of reaction for CO2 + H2 <-> CO + H2O in J/mol."""
    return calc_dG_rxn('rwgs', T)

def calc_dG_rxn_ch(T: float) -> float:
    """Gibbs free energy of reaction for C6H12 <-> C6H6 + 3 H2 in J/mol."""
    return calc_dG_rxn('ch', T)


# ---------------------------------------------------------------------------
# High-Precision Equilibrium Constant Functions
# ---------------------------------------------------------------------------

def calc_keq_meoh(T: float, P: Optional[float] = None, y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                  use_real_gas_eos: Optional[bool] = None) -> float:
    """
    Thermodynamically consistent equilibrium constant Keq1 (bar^-2) for methanol synthesis:
    CO2 + 3 H2 <-> CH3OH + H2O (Delta n = -2).
    Calculated directly from Gibbs free energy ln Keq1 = -Delta G_r,meoh(T) / (R * T).
    Strictly obeys the van 't Hoff identity d(ln Keq1)/dT = Delta H_r,meoh(T) / (R * T^2).
    If use_real_gas_eos is active (and P, y supplied), corrects by fugacity factor K_phi,1.
    """
    K_ideal = calc_kp('meoh', T)
    is_real = use_real_gas_eos if use_real_gas_eos is not None else USE_REAL_GAS_EOS
    if is_real and P is not None and y is not None:
        K_phi = calc_K_phi_1(T, P, y)
        return float(K_ideal / max(K_phi, 1e-15))
    return float(K_ideal)

def calc_keq_rwgs(T: float, P: Optional[float] = None, y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                  use_real_gas_eos: Optional[bool] = None) -> float:
    """
    Thermodynamically consistent equilibrium constant Keq2 (dimensionless) for reverse water-gas shift:
    CO2 + H2 <-> CO + H2O (Delta n = 0).
    Calculated directly from Gibbs free energy ln Keq2 = -Delta G_r,rwgs(T) / (R * T).
    Strictly obeys the van 't Hoff identity d(ln Keq2)/dT = Delta H_r,rwgs(T) / (R * T^2).
    If use_real_gas_eos is active (and P, y supplied), corrects by fugacity factor K_phi,2.
    """
    K_ideal = calc_kp('rwgs', T)
    is_real = use_real_gas_eos if use_real_gas_eos is not None else USE_REAL_GAS_EOS
    if is_real and P is not None and y is not None:
        K_phi = calc_K_phi_2(T, P, y)
        return float(K_ideal / max(K_phi, 1e-15))
    return float(K_ideal)

def calc_kp_dehydrogenation(T: float, P: Optional[float] = None, y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                            use_real_gas_eos: Optional[bool] = None) -> float:
    """
    Thermodynamically consistent equilibrium constant Kp (bar^3) for cyclohexane dehydrogenation:
    C6H12 <-> C6H6 + 3 H2 (Delta n = +3).
    Calculated directly from Gibbs free energy ln Kp = -Delta G_r,ch(T) / (R * T).
    Strictly obeys the van 't Hoff identity d(ln Kp)/dT = Delta H_r,ch(T) / (R * T^2).
    If use_real_gas_eos is active (and P, y supplied), corrects by fugacity factor K_phi,ch.
    """
    K_ideal = calc_kp('ch', T)
    is_real = use_real_gas_eos if use_real_gas_eos is not None else USE_REAL_GAS_EOS
    if is_real and P is not None and y is not None:
        K_phi = calc_K_phi_ch(T, P, y)
        return float(K_ideal / max(K_phi, 1e-15))
    return float(K_ideal)

# Standardized Aliases
calc_keq1 = calc_keq_meoh
calc_keq2 = calc_keq_rwgs
calc_kp_ch = calc_kp_dehydrogenation
Keq1 = calc_keq_meoh
Keq2 = calc_keq_rwgs
Kp_CH = calc_kp_dehydrogenation


# ---------------------------------------------------------------------------
# Exact Analytical Derivatives & van 't Hoff Verification
# ---------------------------------------------------------------------------

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

def verify_van_t_hoff_identity(T_min: float = 300.0, T_max: float = 650.0,
                               n_points: int = 100) -> Dict[str, float]:
    """
    Verifies that the exact van 't Hoff identity R * T^2 * d(ln K_p)/dT == Delta H_r(T)
    holds to machine precision (<= 1e-12 relative error) across 300-650 K.
    Uses complex-step differentiation (step size 1e-30) to compute d(ln K_p)/dT
    without any finite-difference truncation error.
    """
    T_vals = np.linspace(T_min, T_max, n_points)
    max_rel_errors: Dict[str, float] = {}
    h = 1e-30  # complex-step perturbation

    for rxn in ('meoh', 'rwgs', 'ch'):
        max_err = 0.0
        c = DELTA_CP_COEFFS[rxn]
        Tr = T_REF_ENTHALPY
        for T in T_vals:
            T_c = complex(T, h)
            # Complex step analytical enthalpy & entropy
            dH_c = DELTA_H_R_298[rxn] + (
                c[0] * (T_c - Tr)
                + 0.5 * c[1] * (T_c ** 2 - Tr ** 2)
                + (1.0 / 3.0) * c[2] * (T_c ** 3 - Tr ** 3)
                + 0.25 * c[3] * (T_c ** 4 - Tr ** 4)
            )
            dS_c = DELTA_S_R_298[rxn] + (
                c[0] * np.log(T_c / Tr)
                + c[1] * (T_c - Tr)
                + 0.5 * c[2] * (T_c ** 2 - Tr ** 2)
                + (1.0 / 3.0) * c[3] * (T_c ** 3 - Tr ** 3)
            )
            ln_Kp_c = -(dH_c - T_c * dS_c) / (R_GAS * T_c)
            dln_Kp_dT = ln_Kp_c.imag / h

            dH_expected = calc_dH_rxn(rxn, T)
            dH_computed = R_GAS * (T ** 2) * dln_Kp_dT
            rel_err = abs(dH_computed - dH_expected) / abs(dH_expected)
            if rel_err > max_err:
                max_err = rel_err
        max_rel_errors[rxn] = float(max_err)
    return max_rel_errors


# ---------------------------------------------------------------------------
# High-Pressure Real-Gas Peng-Robinson Equation of State (PR-EOS)
# ---------------------------------------------------------------------------

USE_REAL_GAS_EOS: bool = False

def set_use_real_gas_eos(enabled: bool) -> None:
    """Sets global flag for Peng-Robinson real gas fugacity corrections."""
    global USE_REAL_GAS_EOS
    USE_REAL_GAS_EOS = bool(enabled)

def get_use_real_gas_eos() -> bool:
    """Returns current state of Peng-Robinson real gas fugacity flag."""
    return USE_REAL_GAS_EOS

# Standard species order for vector/matrix conversions
SPECIES_TUBE: Tuple[str, ...] = ('H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2')
SPECIES_SHELL: Tuple[str, ...] = ('C6H12', 'C6H6', 'H2', 'Ar')

def calc_pr_pure_params(species: str, T: float) -> Tuple[float, float]:
    """
    Calculates pure component Peng-Robinson parameters a_i(T) [Pa·m^6/mol^2] and b_i [m^3/mol].
    a_i(T) = 0.45724 * (R * T_c)^2 / P_c * [1 + m * (1 - sqrt(T_r))]^2
    b_i    = 0.07780 * R * T_c / P_c
    """
    props = CRIT_PROPS[species]
    Tc = props['Tc']
    Pc_Pa = props['Pc'] * 1e5
    om = props['omega']
    Tr = max(float(T), 1.0) / Tc
    m = 0.37464 + 1.54226 * om - 0.26992 * (om ** 2)
    alpha = (1.0 + m * (1.0 - np.sqrt(Tr))) ** 2
    a_i = 0.45724 * ((R_GAS * Tc) ** 2) / Pc_Pa * alpha
    b_i = 0.07780 * (R_GAS * Tc) / Pc_Pa
    return float(a_i), float(b_i)

def _solve_cubic_pr(A: float, B: float, phase: str = 'gas') -> float:
    """
    Solves the cubic Peng-Robinson equation for compressibility factor Z:
    Z^3 - (1 - B)*Z^2 + (A - 3*B^2 - 2*B)*Z - (A*B - B^2 - B^3) = 0
    Returns the physical root (largest for gas, smallest for liquid).
    """
    p = -(1.0 - B)
    q = A - 3.0 * (B ** 2) - 2.0 * B
    r = -(A * B - (B ** 2) - (B ** 3))

    a1 = q - (p ** 2) / 3.0
    a0 = r - (p * q) / 3.0 + 2.0 * (p ** 3) / 27.0
    D = (a0 / 2.0) ** 2 + (a1 / 3.0) ** 3

    if D <= 0:
        r_pol = np.sqrt(max(-(a1 ** 3) / 27.0, 1e-30))
        cos_phi = np.clip(-a0 / (2.0 * r_pol), -1.0, 1.0)
        phi = np.arccos(cos_phi)
        m = 2.0 * np.sqrt(max(-a1 / 3.0, 0.0))
        x0 = m * np.cos(phi / 3.0)
        x1 = m * np.cos((phi + 2.0 * np.pi) / 3.0)
        x2 = m * np.cos((phi + 4.0 * np.pi) / 3.0)
        roots = [float(x0 - p / 3.0), float(x1 - p / 3.0), float(x2 - p / 3.0)]
    else:
        sqrt_D = np.sqrt(D)
        u = np.cbrt(-a0 / 2.0 + sqrt_D)
        v = np.cbrt(-a0 / 2.0 - sqrt_D)
        roots = [float(u + v - p / 3.0)]

    valid_roots = [z for z in roots if z > B]
    if not valid_roots:
        return 1.0
    if phase.lower() == 'liquid':
        return min(valid_roots)
    return max(valid_roots)

def calc_pr_eos(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                species_list: Optional[Sequence[str]] = None,
                phase: str = 'gas') -> Tuple[float, Dict[str, float]]:
    """
    Computes compressibility factor Z and component fugacity coefficients phi_i using Peng-Robinson EOS.
    P can be provided in Pa (> 500) or bar (<= 500).
    y can be provided as a dict {species: mole_fraction} or array-like.
    """
    P_Pa = float(P) * 1e5 if float(P) <= 500.0 else float(P)
    T_K = max(float(T), 10.0)

    if isinstance(y, dict):
        comps = list(y.keys())
        y_arr = np.array([max(float(y[c]), 0.0) for c in comps], dtype=float)
    else:
        y_arr = np.array(y, dtype=float).ravel()
        if species_list is not None:
            comps = list(species_list)
        elif len(y_arr) == 7:
            comps = list(SPECIES_TUBE)
        elif len(y_arr) == 4:
            comps = list(SPECIES_SHELL)
        elif len(y_arr) == 3:
            comps = ['C6H12', 'C6H6', 'H2']
        else:
            raise ValueError(f"Cannot infer species list for array of length {len(y_arr)}; provide species_list.")

    y_sum = np.sum(y_arr)
    if y_sum > 0:
        y_norm = y_arr / y_sum
    else:
        y_norm = np.ones(len(comps)) / len(comps)

    n = len(comps)
    a_pure = np.zeros(n)
    b_pure = np.zeros(n)
    for i, c in enumerate(comps):
        a_pure[i], b_pure[i] = calc_pr_pure_params(c, T_K)

    # Van der Waals quadratic mixing rules (kij = 0)
    a_matrix = np.sqrt(np.outer(a_pure, a_pure))
    a_mix = float(y_norm @ a_matrix @ y_norm)
    b_mix = float(np.sum(y_norm * b_pure))

    A = a_mix * P_Pa / ((R_GAS * T_K) ** 2)
    B = b_mix * P_Pa / (R_GAS * T_K)

    Z = _solve_cubic_pr(A, B, phase=phase)

    # Rigorous analytical fugacity coefficients
    SQRT2 = 1.4142135623730951
    phi: Dict[str, float] = {}
    for i, c in enumerate(comps):
        sum_ay = float(np.sum(y_norm * a_matrix[i, :]))
        ln_phi_i = (
            (b_pure[i] / b_mix) * (Z - 1.0)
            - np.log(max(Z - B, 1e-15))
            - (A / (2.0 * SQRT2 * B))
            * (2.0 * sum_ay / a_mix - b_pure[i] / b_mix)
            * np.log((Z + (1.0 + SQRT2) * B) / max(Z + (1.0 - SQRT2) * B, 1e-15))
        )
        phi[c] = float(np.exp(ln_phi_i))

    return Z, phi

def calc_compressibility_factor_pr(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                                   species_list: Optional[Sequence[str]] = None,
                                   phase: str = 'gas') -> float:
    """Computes Peng-Robinson compressibility factor Z(T, P, y)."""
    Z, _ = calc_pr_eos(T, P, y, species_list=species_list, phase=phase)
    return Z

def calc_fugacity_coefficients_pr(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                                  species_list: Optional[Sequence[str]] = None,
                                  phase: str = 'gas') -> Dict[str, float]:
    """Computes Peng-Robinson component fugacity coefficients phi_i(T, P, y)."""
    _, phi = calc_pr_eos(T, P, y, species_list=species_list, phase=phase)
    return phi

def calc_K_phi_1(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                 species_list: Optional[Sequence[str]] = None) -> float:
    """
    Reaction 1 (CO2 + 3 H2 <-> CH3OH + H2O):
    K_phi,1 = (phi_MeOH * phi_H2O) / (phi_CO2 * phi_H2^3)
    """
    phi = calc_fugacity_coefficients_pr(T, P, y, species_list=species_list)
    return float((phi['CH3OH'] * phi['H2O']) / (phi['CO2'] * (phi['H2'] ** 3)))

def calc_K_phi_2(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                 species_list: Optional[Sequence[str]] = None) -> float:
    """
    Reaction 2 (CO2 + H2 <-> CO + H2O):
    K_phi,2 = (phi_CO * phi_H2O) / (phi_CO2 * phi_H2)
    """
    phi = calc_fugacity_coefficients_pr(T, P, y, species_list=species_list)
    return float((phi['CO'] * phi['H2O']) / (phi['CO2'] * phi['H2']))

def calc_K_phi_3(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                 species_list: Optional[Sequence[str]] = None) -> float:
    """
    Reaction 3 (CO + 2 H2 <-> CH3OH):
    K_phi,3 = phi_MeOH / (phi_CO * phi_H2^2)
    """
    phi = calc_fugacity_coefficients_pr(T, P, y, species_list=species_list)
    return float(phi['CH3OH'] / (phi['CO'] * (phi['H2'] ** 2)))

def calc_K_phi_ch(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
                  species_list: Optional[Sequence[str]] = None) -> float:
    """
    Reaction 4 (C6H12 <-> C6H6 + 3 H2):
    K_phi,ch = (phi_C6H6 * phi_H2^3) / phi_C6H12
    """
    phi = calc_fugacity_coefficients_pr(T, P, y, species_list=species_list)
    return float((phi['C6H6'] * (phi['H2'] ** 3)) / phi['C6H12'])

def calc_K_phi(T: float, P: float, y: Union[Dict[str, float], np.ndarray, Sequence[float]],
               species_list: Optional[Sequence[str]] = None) -> Dict[str, float]:
    """Computes all reaction fugacity correction factors K_phi for tube or shell mixtures."""
    phi = calc_fugacity_coefficients_pr(T, P, y, species_list=species_list)
    res: Dict[str, float] = {}
    if 'CH3OH' in phi and 'H2O' in phi and 'CO2' in phi and 'H2' in phi:
        res['K_phi_1'] = float((phi['CH3OH'] * phi['H2O']) / (phi['CO2'] * (phi['H2'] ** 3)))
        res['K_phi_meoh'] = res['K_phi_1']
    if 'CO' in phi and 'H2O' in phi and 'CO2' in phi and 'H2' in phi:
        res['K_phi_2'] = float((phi['CO'] * phi['H2O']) / (phi['CO2'] * phi['H2']))
        res['K_phi_rwgs'] = res['K_phi_2']
    if 'CH3OH' in phi and 'CO' in phi and 'H2' in phi:
        res['K_phi_3'] = float(phi['CH3OH'] / (phi['CO'] * (phi['H2'] ** 2)))
        res['K_phi_co'] = res['K_phi_3']
    if 'C6H12' in phi and 'C6H6' in phi and 'H2' in phi:
        res['K_phi_ch'] = float((phi['C6H6'] * (phi['H2'] ** 3)) / phi['C6H12'])
    return res

def calc_effective_keq_meoh(T: float, P: float = 76.98e5,
                            y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                            use_real_gas_eos: Optional[bool] = None) -> float:
    """Helper to evaluate effective Keq1 with optional PR-EOS fugacity correction."""
    if y is None:
        y = {'H2': 0.659, 'CO2': 0.094, 'CO': 0.046, 'CH3OH': 0.005, 'H2O': 0.0004, 'CH4': 0.1026, 'N2': 0.093}
    return calc_keq_meoh(T, P=P, y=y, use_real_gas_eos=use_real_gas_eos)

def calc_effective_keq_rwgs(T: float, P: float = 76.98e5,
                            y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                            use_real_gas_eos: Optional[bool] = None) -> float:
    """Helper to evaluate effective Keq2 with optional PR-EOS fugacity correction."""
    if y is None:
        y = {'H2': 0.659, 'CO2': 0.094, 'CO': 0.046, 'CH3OH': 0.005, 'H2O': 0.0004, 'CH4': 0.1026, 'N2': 0.093}
    return calc_keq_rwgs(T, P=P, y=y, use_real_gas_eos=use_real_gas_eos)

def calc_effective_kp_ch(T: float, P: float = 2.5e5,
                         y: Optional[Union[Dict[str, float], np.ndarray]] = None,
                         use_real_gas_eos: Optional[bool] = None) -> float:
    """Helper to evaluate effective Kp_ch with optional PR-EOS fugacity correction."""
    if y is None:
        y = {'C6H12': 0.95, 'C6H6': 0.02, 'H2': 0.03, 'Ar': 0.0}
    return calc_kp_dehydrogenation(T, P=P, y=y, use_real_gas_eos=use_real_gas_eos)


# ---------------------------------------------------------------------------
# Backwards Compatibility Interface
# ---------------------------------------------------------------------------

def _build_van_t_hoff_params(stoich: Dict[str, float], ln_K_ref: float,
                             T_0: float = 500.0) -> Tuple[float, np.ndarray, float]:
    """Legacy helper for analytical integration parameters."""
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

_VTH_PARAMS_MEOH = _build_van_t_hoff_params(STOICHIOMETRY['meoh'], _LN_KEQ1_REF_500)
_VTH_PARAMS_RWGS = _build_van_t_hoff_params(STOICHIOMETRY['rwgs'], _LN_KEQ2_REF_500)
_VTH_PARAMS_CH   = _build_van_t_hoff_params(STOICHIOMETRY['ch'], _LN_KP_CH_REF_500)

def _eval_van_t_hoff(params: Tuple[float, np.ndarray, float], T: float) -> float:
    """Legacy van 't Hoff evaluator."""
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

"""
Rigorous chemical kinetics module for:
1. Production Mode:
   - Methanol Synthesis & RWGS: Vanden Bussche & Froment (1996) LHHW kinetics over CuO/ZnO/Al2O3.
   - Cyclohexane Dehydrogenation: Peer-reviewed Itoh (1987, AIChE J. 33, 1576-1578) / Kusakabe et al. (2003) /
     Khademi et al. (2009, Eq. 8 & Table 4) rate law in algebraically rationalized, singularity-free form
     with exact Kirchhoff K_p(T) equilibrium from tcr.thermo.
2. Benchmark Mode:
   - Methanol Synthesis: Graaf et al. (1988) dual-site model (Khademi et al., 2009, Tables 1-3).
   - Cyclohexane Dehydrogenation: Itoh (1987) / Khademi Table 4 model with support for both
     literal printed Table 4 mode (typo Kp = 4.89e35 * exp(-3190/T) Pa^3) and Kirchhoff-corrected mode.
3. Catalyst Pellet Diffusion & Effectiveness Factors:
   - Analytical spherical pellet effectiveness factor: eta(Phi) = (3 / Phi) * (coth(Phi) - 1 / Phi).
   - Dynamic axial Thiele modulus Phi(z) and eta(z) calculators for tube and shell pellets.
   - Full support for constant eta sensitivity runs and dynamic Thiele modulus evaluation.

Strict Microscopic Reversibility:
All reversible rate equations formulate the driving force as (1 - eta_eq), where eta_eq = Q / K is the
approach to equilibrium. At exact chemical equilibrium partial pressures (eta_eq == 1.0), rates evaluate
to strictly 0.0 (<= 1e-15 mol/(kg_cat·s)).
"""

import numpy as np
from .thermo import (
    R_GAS, calc_kp_dehydrogenation, calc_keq_meoh, calc_keq_rwgs,
    calc_keq1, calc_keq2, calc_kp_ch, Keq1, Keq2, Kp_CH
)

# Reversibility tolerance for approach to equilibrium
_REV_TOL = 1e-14
_RATE_ZERO_TOL = 1e-15

# Global override for benchmark sub-mode ("khademi_literal" vs "khademi_thermo")
_ITOH_BENCHMARK_MODE_OVERRIDE: str | None = None


def set_itoh_benchmark_mode(mode: str | None) -> None:
    """Sets a global override for Itoh benchmark mode ('khademi_literal' or 'khademi_thermo')."""
    global _ITOH_BENCHMARK_MODE_OVERRIDE
    _ITOH_BENCHMARK_MODE_OVERRIDE = mode


def get_itoh_benchmark_mode() -> str | None:
    """Returns the current global override for Itoh benchmark mode."""
    return _ITOH_BENCHMARK_MODE_OVERRIDE


# ---------------------------------------------------------------------------
# Thiele Modulus & Effectiveness Factor Calculator (Swarm Critic-06)
# ---------------------------------------------------------------------------

def calc_effectiveness_factor_spherical(Phi: float | np.ndarray) -> float | np.ndarray:
    """
    Evaluates the analytical internal effectiveness factor for a spherical catalyst pellet:
        eta = (3 / Phi) * (coth(Phi) - 1 / Phi)

    Implementation features:
    - Small Phi limit (Phi < 1e-3): series expansion eta = 1 - Phi^2 / 15 + 2 * Phi^4 / 315
      to prevent numerical division-by-zero and catastrophic cancellation.
    - Large Phi limit (Phi > 35): asymptotic approximation eta = (3 / Phi) * (1 - 1 / Phi)
      to avoid overflow in hyperbolic functions.
    - Intermediate range: exact hyperbolic formulation via 1.0 / tanh(Phi).
    - Strictly clipped to [0.0, 1.0].
    """
    phi_arr = np.asarray(Phi, dtype=float)
    scalar_input = phi_arr.ndim == 0
    phi_flat = np.atleast_1d(np.maximum(phi_arr, 0.0))
    eta = np.zeros_like(phi_flat)

    # 1. Small Phi limit
    mask_small = phi_flat < 1e-3
    if np.any(mask_small):
        p = phi_flat[mask_small]
        eta[mask_small] = 1.0 - (p ** 2) / 15.0 + 2.0 * (p ** 4) / 315.0

    # 2. Large Phi limit
    mask_large = phi_flat > 35.0
    if np.any(mask_large):
        p = phi_flat[mask_large]
        eta[mask_large] = (3.0 / p) * (1.0 - 1.0 / p)

    # 3. Intermediate range
    mask_mid = ~mask_small & ~mask_large
    if np.any(mask_mid):
        p = phi_flat[mask_mid]
        coth_p = 1.0 / np.tanh(p)
        eta[mask_mid] = (3.0 / p) * (coth_p - 1.0 / p)

    eta = np.clip(eta, 0.0, 1.0)
    return float(eta[0]) if scalar_input else eta


# Standard Aliases
calc_thiele_eta = calc_effectiveness_factor_spherical
calc_eta_spherical = calc_effectiveness_factor_spherical


def calc_thiele_modulus_general(r_vol: float, C_reactant: float, Rp: float,
                                Deff: float = 2.5e-6) -> float:
    """
    Evaluates the generalized Thiele modulus for an intrinsic volumetric reaction rate:
        Phi = Rp * sqrt( max(r_vol, 0) / (Deff * max(C_reactant, 1e-15)) )
    """
    r_safe = max(float(r_vol), 0.0)
    c_safe = max(float(C_reactant), 1e-15)
    k_pseudo1 = r_safe / c_safe
    return float(Rp * np.sqrt(k_pseudo1 / max(Deff, 1e-12)))


def calc_thiele_tube(T_exo: float, P_exo: float, y_exo: np.ndarray,
                     dp: float = 0.00547, rho_b: float = 1140.0, eps: float = 0.39,
                     Deff: float = 2.5e-6) -> tuple[float, float]:
    """
    Calculates dynamic Thiele modulus Phi and effectiveness factor eta for tube catalyst (CuO/ZnO/Al2O3):
    - Rp = dp / 2
    - rho_p = rho_b / (1 - eps)
    - Evaluates intrinsic methanol synthesis rate (eta=1.0)
    - Returns (Phi, eta).
    """
    T = max(float(T_exo), 200.0)
    P_Pa = max(float(P_exo) * 1e5 if P_exo < 1000.0 else float(P_exo), 1e4)
    y_safe = np.maximum(y_exo, 1e-15)

    Rp = max(dp, 1e-5) / 2.0
    rho_p = rho_b / max(1.0 - eps, 0.1)

    # Intrinsic rate with eta=1.0
    r_meoh_int, _ = kinetics_vbf_production(T, P_Pa, y_safe, eta_exo=1.0)
    r_vol = max(r_meoh_int * rho_p, 0.0)

    # Key limiting reactant: CO2 concentration
    C_CO2 = max(y_safe[1] * P_Pa / (R_GAS * T), 1e-15)
    Phi = calc_thiele_modulus_general(r_vol, C_CO2, Rp, Deff=Deff)
    eta = calc_effectiveness_factor_spherical(Phi)
    return float(Phi), float(eta)


def calc_thiele_shell(T_endo: float, P_endo: float, y_endo: np.ndarray,
                      dp: float = 0.00355, rho_b: float = 1077.0, eps: float = 0.40,
                      Deff: float = 2.5e-6, rho_cat_endo: float = 1077.0) -> tuple[float, float]:
    """
    Calculates dynamic Thiele modulus Phi and effectiveness factor eta for shell catalyst (Pt/Al2O3):
    - Rp = dp / 2
    - rho_p = rho_b / (1 - eps)
    - Evaluates intrinsic cyclohexane dehydrogenation rate (eta=1.0)
    - Returns (Phi, eta).
    """
    T = max(float(T_endo), 200.0)
    P_Pa = max(float(P_endo) * 1e5 if P_endo < 1000.0 else float(P_endo), 1e4)
    y_safe = np.maximum(y_endo, 1e-15)

    Rp = max(dp, 1e-5) / 2.0
    rho_p = rho_b / max(1.0 - eps, 0.1)

    # Intrinsic rate with eta=1.0 in mass units
    r_ch_int = kinetics_ch_production(T, P_Pa, y_safe, eta_endo=1.0, rho_cat_endo=rho_cat_endo)
    r_vol = max(r_ch_int * rho_p, 0.0)

    # Key reactant: Cyclohexane concentration
    C_CH = max(y_safe[0] * P_Pa / (R_GAS * T), 1e-15)
    Phi = calc_thiele_modulus_general(r_vol, C_CH, Rp, Deff=Deff)
    eta = calc_effectiveness_factor_spherical(Phi)
    return float(Phi), float(eta)


# ---------------------------------------------------------------------------
# Production Kinetics
# ---------------------------------------------------------------------------

def kinetics_vbf_production(T_exo: float, P_exo: float, y_exo: np.ndarray,
                            eta_exo: float | str | None = 0.35) -> tuple[float, float]:
    """
    Vanden Bussche & Froment (1996) LHHW rate equations for methanol synthesis from CO2
    and reverse water-gas shift (RWGS).
    Units: r_meoh, r_rwgs in mol/(kg_cat·s).
    y_exo: [H2, CO2, CO, CH3OH, H2O, CH4, N2]

    eta_exo:
    - float: constant effectiveness factor (default 0.35)
    - 'dynamic': evaluates dynamic Thiele effectiveness factor eta(z)
    - callable: evaluates custom eta function
    """
    P_Pa = max(float(P_exo) * 1e5 if P_exo < 1000.0 else float(P_exo), 1e4)
    P_bar = P_Pa / 1e5
    y_safe = np.maximum(y_exo, 1e-15)
    p = y_safe * P_bar

    pH2   = max(p[0], 1e-15)
    pCO2  = max(p[1], 1e-15)
    pCO   = max(p[2], 1e-15)
    pMeOH = max(p[3], 1e-15)
    pH2O  = max(p[4], 1e-15)
    T = max(T_exo, 273.15)

    # Rate constants (mol/(kg_cat·s·bar^n))
    k1 = 1.07 * np.exp(36696.0 / (R_GAS * T))
    k2 = 1.22e10 * np.exp(-94765.0 / (R_GAS * T))

    # Adsorption coefficients (VBF 1996, positive exponent for exothermic H2 adsorption)
    sqrt_KH2 = 0.499 * np.exp(17197.0 / (R_GAS * T))
    KH2O = 6.62e-11 * np.exp(124119.0 / (R_GAS * T))
    K_ratio = 3453.38

    # Thermodynamic equilibrium constants strictly obeying van 't Hoff equation: d(ln K)/dT = Delta H_r(T) / (R * T^2)
    # Passed P and y to enable real-gas Peng-Robinson fugacity corrections K_phi when active
    keq1 = calc_keq1(T, P=P_bar, y=y_safe)  # bar^-2
    keq2 = calc_keq2(T, P=P_bar, y=y_safe)  # dimensionless

    # Reversible driving forces via (1 - eta) formulation ensuring microscopic reversibility:
    eta1 = (pMeOH * pH2O) / (max(keq1, 1e-30) * pCO2 * (pH2 ** 3))
    df1 = 0.0 if abs(1.0 - eta1) <= _REV_TOL else pCO2 * pH2 * (1.0 - eta1)

    eta2 = (pCO * pH2O) / (max(keq2, 1e-30) * pCO2 * pH2)
    df2 = 0.0 if abs(1.0 - eta2) <= _REV_TOL else pCO2 * (1.0 - eta2)

    # Langmuir-Hinshelwood denominator
    denom = 1.0 + K_ratio * (pH2O / pH2) + sqrt_KH2 * np.sqrt(pH2) + KH2O * pH2O

    # Effectiveness factor evaluation
    if eta_exo is None:
        eff_eta = 0.35
    elif isinstance(eta_exo, str) and eta_exo.lower() == "dynamic":
        _, eff_eta = calc_thiele_tube(T, P_Pa, y_safe)
    elif callable(eta_exo):
        eff_eta = float(eta_exo(T, P_Pa, y_safe))
    else:
        eff_eta = float(eta_exo)

    r_meoh = eff_eta * k1 * df1 / (denom ** 3)
    r_rwgs = eff_eta * k2 * df2 / denom

    if abs(r_meoh) <= _RATE_ZERO_TOL:
        r_meoh = 0.0
    if abs(r_rwgs) <= _RATE_ZERO_TOL:
        r_rwgs = 0.0

    return float(r_meoh), float(r_rwgs)


def kinetics_ch_production(T_endo: float, P_endo: float, y_endo: np.ndarray,
                           eta_endo: float | str | None = 1.0,
                           rho_cat_endo: float = 1077.0) -> float:
    """
    Reversible catalytic cyclohexane dehydrogenation over Pt/Al2O3:
        C6H12 <-> C6H6 + 3 H2

    Unified on peer-reviewed Itoh (1987, AIChE J. 33, 1576–1578) / Kusakabe et al. (2003) /
    Khademi et al. (2009, Eq. 8 & Table 4) rate law in algebraically rationalized,
    singularity-free form:

        r_CH = eta_endo(z) * k(T) * [K_p(T) * P_CH - P_B * P_H2^3] / [P_H2^3 + K_B(T) * K_p(T) * P_CH]

    Parameters:
        k(T)   = 0.221 * exp(-4270 / T) mol / (m^3_cat * s * Pa)
        K_B(T) = 2.03e-10 * exp(6270 / T) Pa^-1
        K_p(T) = calc_kp_dehydrogenation(T) * 1e15 Pa^3 (exact Kirchhoff thermo)
        rho_cat_endo = 1077.0 kg_cat / m^3 (catalyst bed bulk density, converting rate to mol/(kg_cat·s))

    eta_endo:
        - float: constant effectiveness factor (default 1.0)
        - 'dynamic': evaluates dynamic Thiele effectiveness factor eta(z)
        - callable: evaluates custom eta function

    Units:
        r_ch in mol/(kg_cat·s).
        y_endo: [C6H12, C6H6, H2, (inert Ar/N2)]
    """
    P_Pa = max(float(P_endo) * 1e5 if P_endo < 1000.0 else float(P_endo), 1e4)
    y_safe = np.maximum(y_endo, 1e-15)
    p_Pa = y_safe * P_Pa

    P_CH = max(p_Pa[0], 1e-15)
    P_B  = max(p_Pa[1], 1e-15)
    P_H2 = max(p_Pa[2], 1e-15)
    T = max(T_endo, 200.0)

    # Itoh (1987) / Khademi Table 4 rate constant: k in mol/(m^3_cat·s·Pa)
    k_vol = 0.221 * np.exp(-4270.0 / T)
    k_mass = k_vol / rho_cat_endo if rho_cat_endo > 0.0 else k_vol

    # Benzene adsorption constant: K_B in Pa^-1
    K_B = 2.03e-10 * np.exp(6270.0 / T)

    # Exact Kirchhoff thermodynamic equilibrium constant in Pa^3: Kp_bar3 * (1e5)^3
    P_bar = P_Pa / 1e5
    Kp_bar3 = calc_kp_dehydrogenation(T, P=P_bar, y=y_safe)
    Kp_Pa3 = Kp_bar3 * 1e15

    # Reversible driving force via approach to equilibrium ensuring microscopic reversibility
    eta_eq = (P_B * (P_H2 ** 3)) / max(Kp_Pa3 * P_CH, 1e-30)
    if abs(1.0 - eta_eq) <= _REV_TOL:
        num = 0.0
    else:
        num = Kp_Pa3 * P_CH - P_B * (P_H2 ** 3)

    # Algebraically rationalized, singularity-free denominator (no division by P_H2^3)
    denom = (P_H2 ** 3) + K_B * Kp_Pa3 * P_CH
    denom = max(denom, 1e-30)

    # Effectiveness factor evaluation
    if eta_endo is None:
        eff_eta = 1.0
    elif isinstance(eta_endo, str) and eta_endo.lower() == "dynamic":
        _, eff_eta = calc_thiele_shell(T, P_Pa, y_safe, rho_cat_endo=rho_cat_endo)
    elif callable(eta_endo):
        eff_eta = float(eta_endo(T, P_Pa, y_safe))
    else:
        eff_eta = float(eta_endo)

    r_ch = eff_eta * k_mass * (num / denom)
    if abs(r_ch) <= _RATE_ZERO_TOL:
        r_ch = 0.0
    return float(r_ch)


# ---------------------------------------------------------------------------
# Benchmark Kinetics (Khademi et al., 2009)
# ---------------------------------------------------------------------------

def kinetics_graaf_benchmark(T_exo: float, P_exo: float, y_exo: np.ndarray,
                             eta_exo: float | str | None = 1.0) -> tuple[float, float, float]:
    """
    Graaf et al. (1988) kinetics from Khademi et al. (2009), Tables 1-3.
    Reactions:
      1. CO + 2 H2 <-> CH3OH
      2. CO2 + H2  <-> CO + H2O
      3. CO2 + 3 H2 <-> CH3OH + H2O
    Units: rates in mol/(kg_cat·s).
    """
    P_Pa = max(float(P_exo) * 1e5 if P_exo < 1000.0 else float(P_exo), 1e4)
    P_bar = P_Pa / 1e5
    y_safe = np.maximum(y_exo, 1e-15)
    f = y_safe * P_bar  # approximate fugacity with partial pressure in bar

    fH2   = max(f[0], 1e-15)
    fCO2  = max(f[1], 1e-15)
    fCO   = max(f[2], 1e-15)
    fMeOH = max(f[3], 1e-15)
    fH2O  = max(f[4], 1e-15)
    T = max(T_exo, 273.15)

    # Kinetic rate constants (Graaf 1988 / Khademi Table 2)
    k1 = 4.39e7 * np.exp(-77200.0 / (R_GAS * T))
    k2 = 7.12e8 * np.exp(-106900.0 / (R_GAS * T))
    k3 = 1.07e11 * np.exp(-124900.0 / (R_GAS * T))

    # Adsorption constants
    K_CO  = 2.16e-5 * np.exp(46800.0 / (R_GAS * T))
    K_CO2 = 7.05e-7 * np.exp(61700.0 / (R_GAS * T))
    K_H2O_over_sqrt_KH2 = 6.37e-9 * np.exp(84000.0 / (R_GAS * T))

    # Equilibrium constants (Khademi Table 3, in bar^-2, dimensionless, bar^-2)
    Kp1 = 10.0 ** (5139.0 / T - 12.621)
    Kp2 = 10.0 ** (-2073.0 / T + 2.029)
    Kp3 = 10.0 ** (3066.0 / T - 10.592)

    denom1 = 1.0 + K_CO * fCO + K_CO2 * fCO2
    denom2 = np.sqrt(fH2) + K_H2O_over_sqrt_KH2 * fH2O
    denom = denom1 * denom2

    # Reversible driving forces via (1 - eta) formulation ensuring microscopic reversibility
    eta1 = fMeOH / (max(Kp1, 1e-30) * fCO * (fH2 ** 2))
    df1 = 0.0 if abs(1.0 - eta1) <= _REV_TOL else fCO * (fH2 ** 1.5) * (1.0 - eta1)

    eta2 = (fCO * fH2O) / (max(Kp2, 1e-30) * fCO2 * fH2)
    df2 = 0.0 if abs(1.0 - eta2) <= _REV_TOL else fCO2 * fH2 * (1.0 - eta2)

    eta3 = (fMeOH * fH2O) / (max(Kp3, 1e-30) * fCO2 * (fH2 ** 3))
    df3 = 0.0 if abs(1.0 - eta3) <= _REV_TOL else fCO2 * (fH2 ** 1.5) * (1.0 - eta3)

    if eta_exo is None:
        eff_eta = 1.0
    elif isinstance(eta_exo, str) and eta_exo.lower() == "dynamic":
        _, eff_eta = calc_thiele_tube(T, P_Pa, y_safe)
    elif callable(eta_exo):
        eff_eta = float(eta_exo(T, P_Pa, y_safe))
    else:
        eff_eta = float(eta_exo)

    r1 = eff_eta * k1 * K_CO * df1 / denom
    r2 = eff_eta * k2 * K_CO2 * df2 / denom
    r3 = eff_eta * k3 * K_CO2 * df3 / denom

    if abs(r1) <= _RATE_ZERO_TOL:
        r1 = 0.0
    if abs(r2) <= _RATE_ZERO_TOL:
        r2 = 0.0
    if abs(r3) <= _RATE_ZERO_TOL:
        r3 = 0.0

    return float(r1), float(r2), float(r3)


def kinetics_itoh_benchmark(T_endo: float, P_endo: float, y_endo: np.ndarray,
                            mode: str = "khademi_literal",
                            eta_endo: float | str | None = 1.0,
                            rho_cat_endo: float = 1077.0) -> float:
    """
    Itoh (1987) / Khademi et al. (2009, Eq. 8 & Table 4) dehydrogenation kinetics in
    algebraically rationalized, singularity-free form:
        r_CH = eta_endo * k(T) * [K_p(T) * P_CH - P_B * P_H2^3] / [P_H2^3 + K_B(T) * K_p(T) * P_CH]

    Parameters:
        k(T)   = 0.221 * exp(-4270 / T) mol / (m^3_cat * s * Pa)
        K_B(T) = 2.03e-10 * exp(6270 / T) Pa^-1
        rho_cat_endo = 1077.0 kg/m^3 (converts rate to mol/(kg_cat·s))

    mode:
      - 'khademi_literal': uses printed Table 4 formula: Kp = 4.89e35 * exp(-3190 / T) Pa^3
      - 'khademi_thermo' (or 'kirchhoff'): uses exact Kirchhoff thermodynamic Kp converted to Pa^3
    """
    effective_mode = _ITOH_BENCHMARK_MODE_OVERRIDE if _ITOH_BENCHMARK_MODE_OVERRIDE is not None else mode

    P_Pa = max(float(P_endo) * 1e5 if P_endo < 1000.0 else float(P_endo), 1e4)
    y_safe = np.maximum(y_endo, 1e-15)
    p_Pa = y_safe * P_Pa

    P_CH = max(p_Pa[0], 1e-15)
    P_B  = max(p_Pa[1], 1e-15)
    P_H2 = max(p_Pa[2], 1e-15)
    T = max(T_endo, 200.0)

    # Khademi Table 4 rate constant: k in mol/(m^3_cat·s·Pa)
    k_vol = 0.221 * np.exp(-4270.0 / T)
    k_mass = k_vol / rho_cat_endo if rho_cat_endo > 0.0 else k_vol

    # Benzene adsorption constant: K_B in Pa^-1
    K_B = 2.03e-10 * np.exp(6270.0 / T)

    if effective_mode in ("khademi_literal", "literal"):
        # Literal formula as printed in Khademi (2009) Table 4 (contains typo: -3190/T)
        Kp_Pa3 = 4.89e35 * np.exp(-3190.0 / T)
    else:
        # Thermodynamically consistent Kirchhoff Kp in Pa^3: Kp_bar3 * (1e5)^3
        P_bar = P_Pa / 1e5
        Kp_bar3 = calc_kp_dehydrogenation(T, P=P_bar, y=y_safe)
        Kp_Pa3 = Kp_bar3 * 1e15

    # Reversible driving force via approach to equilibrium
    eta_eq = (P_B * (P_H2 ** 3)) / max(Kp_Pa3 * P_CH, 1e-30)
    if abs(1.0 - eta_eq) <= _REV_TOL:
        num = 0.0
    else:
        num = Kp_Pa3 * P_CH - P_B * (P_H2 ** 3)

    # Algebraically rationalized denominator
    denom = (P_H2 ** 3) + K_B * Kp_Pa3 * P_CH
    denom = max(denom, 1e-30)

    # Effectiveness factor evaluation
    if eta_endo is None:
        eff_eta = 1.0
    elif isinstance(eta_endo, str) and eta_endo.lower() == "dynamic":
        _, eff_eta = calc_thiele_shell(T, P_Pa, y_safe, rho_cat_endo=rho_cat_endo)
    elif callable(eta_endo):
        eff_eta = float(eta_endo(T, P_Pa, y_safe))
    else:
        eff_eta = float(eta_endo)

    r_ch = eff_eta * k_mass * (num / denom)
    if abs(r_ch) <= _RATE_ZERO_TOL:
        r_ch = 0.0
    return float(r_ch)


# ---------------------------------------------------------------------------
# Axial Thiele Profile Evaluator
# ---------------------------------------------------------------------------

def calc_axial_thiele_profile(sol, geom, cat_tube=None, cat_shell=None,
                              Deff_tube: float = 2.5e-6, Deff_shell: float = 2.5e-6) -> dict:
    """
    Evaluates continuous axial profiles of Thiele modulus Phi(z) and effectiveness factor eta(z)
    for tube and shell along the reactor axis from a solved reactor solution object.
    """
    z_grid = getattr(sol, 't', None)
    if z_grid is None:
        raise ValueError("Solution object must contain axial grid attribute 't'")

    dp_tube = getattr(cat_tube, 'dp', 0.00547) if cat_tube else 0.00547
    rho_b_tube = getattr(cat_tube, 'rho_b', 1140.0) if cat_tube else 1140.0
    eps_tube = getattr(cat_tube, 'eps', 0.39) if cat_tube else 0.39

    dp_shell = getattr(cat_shell, 'dp', 0.00355) if cat_shell else 0.00355
    rho_b_shell = getattr(cat_shell, 'rho_b', 1077.0) if cat_shell else 1077.0
    eps_shell = getattr(cat_shell, 'eps', 0.40) if cat_shell else 0.40

    m = len(z_grid)
    phi_tube = np.zeros(m)
    eta_tube = np.zeros(m)
    phi_shell = np.zeros(m)
    eta_shell = np.zeros(m)

    num_states = sol.y.shape[0]
    is_tcr = num_states >= 14
    len_endo = (num_states - 11) if is_tcr else 0

    for j in range(m):
        # Tube state
        F_t = np.maximum(sol.y[0:7, j], 1e-15)
        y_t = F_t / np.sum(F_t)
        T_t = sol.y[7 + len_endo, j] if is_tcr else sol.y[7, j]
        P_t = sol.y[9 + len_endo, j] if is_tcr else sol.y[8, j]

        phi_tube[j], eta_tube[j] = calc_thiele_tube(
            T_t, P_t, y_t, dp=dp_tube, rho_b=rho_b_tube, eps=eps_tube, Deff=Deff_tube
        )

        if is_tcr:
            # Shell state
            F_s = np.maximum(sol.y[7:7+len_endo, j], 1e-15)
            y_s = F_s / np.sum(F_s)
            T_s = sol.y[8 + len_endo, j]
            P_s = sol.y[10 + len_endo, j]

            phi_shell[j], eta_shell[j] = calc_thiele_shell(
                T_s, P_s, y_s, dp=dp_shell, rho_b=rho_b_shell, eps=eps_shell,
                Deff=Deff_shell, rho_cat_endo=rho_b_shell
            )

    return {
        "z": z_grid,
        "Phi_tube": phi_tube,
        "eta_tube": eta_tube,
        "Phi_shell": phi_shell if is_tcr else None,
        "eta_shell": eta_shell if is_tcr else None,
        "eta_tube_mean": float(np.mean(eta_tube)),
        "eta_tube_min": float(np.min(eta_tube)),
        "tube_eta_mean": float(np.mean(eta_tube)),
        "tube_eta_min": float(np.min(eta_tube)),
        "eta_shell_mean": float(np.mean(eta_shell)) if is_tcr else None,
        "eta_shell_min": float(np.min(eta_shell)) if is_tcr else None,
        "shell_eta_mean": float(np.mean(eta_shell)) if is_tcr else None,
        "shell_eta_min": float(np.min(eta_shell)) if is_tcr else None,
    }

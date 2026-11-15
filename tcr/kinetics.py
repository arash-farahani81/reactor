"""
Rigorous chemical kinetics module for:
1. Production Mode:
   - Methanol Synthesis & RWGS: Vanden Bussche & Froment (1996) LHHW kinetics over CuO/ZnO/Al2O3.
   - Cyclohexane Dehydrogenation: Reversible LHHW kinetics over Pt/Al2O3 with van 't Hoff equilibrium.
2. Benchmark Mode:
   - Methanol Synthesis: Graaf et al. (1988) dual-site model (Khademi et al., 2009, Tables 1-3).
   - Cyclohexane Dehydrogenation: Itoh (1987) / Khademi Table 4 model (both literal and thermo Kp).
"""

import numpy as np
from .thermo import R_GAS, calc_kp_dehydrogenation

def kinetics_vbf_production(T_exo: float, P_exo: float, y_exo: np.ndarray,
                            eta_exo: float = 0.35) -> tuple[float, float]:
    """
    Vanden Bussche & Froment (1996) LHHW rate equations for methanol synthesis from CO2
    and reverse water-gas shift (RWGS).
    Units: r_meoh, r_rwgs in mol/(kg_cat·s).
    y_exo: [H2, CO2, CO, CH3OH, H2O, CH4, N2]
    """
    P_bar = max(P_exo, 1e4) / 1e5
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

    # Thermodynamic equilibrium constants
    log10_Keq1 = 3066.0 / T - 10.592
    Keq1 = 10.0 ** log10_Keq1  # bar^-2
    log10_Keq2 = -2073.0 / T + 2.029
    Keq2 = 10.0 ** log10_Keq2  # dimensionless

    # Driving forces
    df1 = pCO2 * pH2 - (1.0 / Keq1) * (pMeOH * pH2O / (pH2 ** 2))
    df2 = pCO2 - (1.0 / Keq2) * (pCO * pH2O / pH2)

    # Langmuir-Hinshelwood denominator
    denom = 1.0 + K_ratio * (pH2O / pH2) + sqrt_KH2 * np.sqrt(pH2) + KH2O * pH2O

    r_meoh = eta_exo * k1 * df1 / (denom ** 3)
    r_rwgs = eta_exo * k2 * df2 / denom

    return float(r_meoh), float(r_rwgs)

def kinetics_ch_production(T_endo: float, P_endo: float, y_endo: np.ndarray,
                           eta_endo: float = 1.0) -> float:
    """
    Reversible Langmuir-Hinshelwood catalytic cyclohexane dehydrogenation over Pt/Al2O3:
    C6H12 <-> C6H6 + 3 H2
    Units: r_ch in mol/(kg_cat·s).
    y_endo: [C6H12, C6H6, H2, (inert)]
    """
    P_bar = max(P_endo, 1e4) / 1e5
    y_safe = np.maximum(y_endo, 1e-15)
    p = y_safe * P_bar

    pCH = max(p[0], 1e-15)
    pB  = max(p[1], 1e-15)
    pH2 = max(p[2], 1e-15)
    T = max(T_endo, 273.15)

    # Arrhenius rate constant (mol/(kg_cat·s·bar))
    k_ch0 = 2.45e4
    Ea_ch = 68500.0  # J/mol
    k_ch = k_ch0 * np.exp(-Ea_ch / (R_GAS * T))

    # Thermodynamically consistent van 't Hoff equilibrium constant (bar^3)
    Kp = calc_kp_dehydrogenation(T)

    df_ch = pCH - (pB * (pH2 ** 3)) / max(Kp, 1e-10)

    # Benzene product inhibition
    K_ads_B = 0.85 * np.exp(12000.0 / (R_GAS * T))
    denom_ch = 1.0 + K_ads_B * pB

    r_ch = eta_endo * k_ch * df_ch / denom_ch
    return float(r_ch)

def kinetics_graaf_benchmark(T_exo: float, P_exo: float, y_exo: np.ndarray,
                             eta_exo: float = 1.0) -> tuple[float, float, float]:
    """
    Graaf et al. (1988) kinetics from Khademi et al. (2009), Tables 1-3.
    Reactions:
      1. CO + 2 H2 <-> CH3OH
      2. CO2 + H2  <-> CO + H2O
      3. CO2 + 3 H2 <-> CH3OH + H2O
    Units: rates in mol/(kg_cat·s).
    """
    P_bar = max(P_exo, 1e4) / 1e5
    y_safe = np.maximum(y_exo, 1e-15)
    f = y_safe * P_bar  # approximate fugacity with partial pressure in bar

    fH2   = max(f[0], 1e-15)
    fCO2  = max(f[1], 1e-15)
    fCO   = max(f[2], 1e-15)
    fMeOH = max(f[3], 1e-15)
    fH2O  = max(f[4], 1e-15)
    T = max(T_exo, 273.15)

    # Kinetic rate constants (Graaf 1988 / Khademi Table 2)
    # k = A * exp(-E/RT)
    k1 = 4.39e7 * np.exp(-77200.0 / (R_GAS * T))
    k2 = 7.12e8 * np.exp(-106900.0 / (R_GAS * T))
    k3 = 1.07e11 * np.exp(-124900.0 / (R_GAS * T))

    # Adsorption constants
    K_CO  = 2.16e-5 * np.exp(46800.0 / (R_GAS * T))
    K_CO2 = 7.05e-7 * np.exp(61700.0 / (R_GAS * T))
    K_H2O_over_sqrt_KH2 = 6.37e-9 * np.exp(84000.0 / (R_GAS * T))

    # Equilibrium constants (Khademi Table 3, in bar^-2, dimensionless, bar^-2)
    # log10(Kp1) = 5139/T - 12.621
    # log10(Kp2) = -2073/T + 2.029
    # log10(Kp3) = 3066/T - 10.592
    Kp1 = 10.0 ** (5139.0 / T - 12.621)
    Kp2 = 10.0 ** (-2073.0 / T + 2.029)
    Kp3 = 10.0 ** (3066.0 / T - 10.592)

    denom1 = 1.0 + K_CO * fCO + K_CO2 * fCO2
    denom2 = np.sqrt(fH2) + K_H2O_over_sqrt_KH2 * fH2O
    denom = denom1 * denom2

    df1 = fCO * (fH2 ** 1.5) - fMeOH / (np.sqrt(fH2) * Kp1)
    r1 = eta_exo * k1 * K_CO * df1 / denom

    df2 = fCO2 * fH2 - (fCO * fH2O) / Kp2
    r2 = eta_exo * k2 * K_CO2 * df2 / denom

    df3 = fCO2 * (fH2 ** 1.5) - (fMeOH * fH2O) / ((fH2 ** 1.5) * Kp3)
    r3 = eta_exo * k3 * K_CO2 * df3 / denom

    return float(r1), float(r2), float(r3)

def kinetics_itoh_benchmark(T_endo: float, P_endo: float, y_endo: np.ndarray,
                            mode: str = "khademi_literal", eta_endo: float = 1.0) -> float:
    """
    Itoh (1987) dehydrogenation kinetics from Khademi et al. (2009), Table 4.
    r = k * [P_CH - (P_B * P_H2^3) / Kp]
    k = 4.09e5 * exp(-72400 / RT) mol/(kg_cat·s·Pa) -> converted to bar
    mode:
      - 'khademi_literal': uses printed Table 4 formula: Kp = 4.89e35 * exp(-3190 / T) Pa^3
      - 'khademi_thermo': uses van 't Hoff formula converted to Pa^3
    """
    P_Pa = max(P_endo, 1e4)
    y_safe = np.maximum(y_endo, 1e-15)
    p_Pa = y_safe * P_Pa

    P_CH = max(p_Pa[0], 1e-15)
    P_B  = max(p_Pa[1], 1e-15)
    P_H2 = max(p_Pa[2], 1e-15)
    T = max(T_endo, 273.15)

    # k in mol/(kg_cat·s·Pa) (Khademi Table 4)
    k = 4.09e5 * np.exp(-72400.0 / (R_GAS * T))

    if mode == "khademi_literal":
        # Literal formula as printed in Khademi (2009) Table 4 (contains typo: -3190/T)
        Kp_Pa3 = 4.89e35 * np.exp(-3190.0 / T)
    else:
        # Thermodynamically consistent Kp in Pa^3: Kp_bar3 * (1e5)^3
        Kp_bar3 = calc_kp_dehydrogenation(T)
        Kp_Pa3 = Kp_bar3 * 1e15

    df = P_CH - (P_B * (P_H2 ** 3)) / max(Kp_Pa3, 1e-30)
    r_ch = eta_endo * k * df
    return float(r_ch)

"""
=============================================================================
Thermally Coupled Shell-and-Tube Reactor Simulation
=============================================================================
Exothermic side (Tube):  Methanol Synthesis over CuO/ZnO/Al2O3
Endothermic side (Shell): Cyclohexane Dehydrogenation over Pt/Al2O3

Based on:
  - Khademi, Setoodeh, Rahimpour, Jahanmiri (2009), Int. J. Hydrogen Energy
  - Vanden Bussche & Froment (1996), J. Catalysis, 161, 1-10
  
Author: Reactor Design Course Project
=============================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
import matplotlib.pyplot as plt
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: PHYSICAL CONSTANTS & REACTOR GEOMETRY
# ═══════════════════════════════════════════════════════════════════════════════

R_gas = 8.314          # J/(mol·K)  Universal gas constant

# --- Reactor Geometry ---
L = 7.022              # m          Reactor length
D_tube_in = 0.038      # m          Tube inner diameter (methanol side)
D_tube_out = 0.042     # m          Tube outer diameter
D_shell_in = 0.06      # m          Shell inner diameter (effective per tube)
N_tubes = 2962         # -          Number of tubes (industrial scale, but we model per-tube)

A_tube = np.pi / 4 * D_tube_in**2        # m²  Cross-sectional area of tube
A_shell = np.pi / 4 * (D_shell_in**2 - D_tube_out**2)  # m²  Annular shell area per tube
A_shell_flow = A_shell * 6.5             # m²  Effective bundle pitch flow area per tube

# --- Catalyst Properties ---
# Methanol synthesis (tube side)
dp_exo = 5.47e-3       # m          Catalyst particle diameter
rho_b_exo = 1140.0     # kg/m³      Bulk catalyst density
eps_exo = 0.39          # -          Bed void fraction

# Cyclohexane dehydrogenation (shell side)
dp_endo = 3.55e-3      # m          Catalyst particle diameter
rho_b_endo = 1077.0    # kg/m³      Bulk catalyst density
eps_endo = 0.39         # -          Bed void fraction

# --- Heat Transfer ---
U_overall = 666.0      # W/(m²·K)   Dynamic overall heat transfer coefficient (Khademi Table 5)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: FEED CONDITIONS
# ═══════════════════════════════════════════════════════════════════════════════

# --- Tube Side (Methanol Synthesis) ---
# Feed composition (mole fractions)
# Components: [H2, CO2, CO, CH3OH, H2O, CH4, N2]
y_exo_0 = np.array([0.6590, 0.0940, 0.0460, 0.0050, 0.0004, 0.1026, 0.0930])
comp_exo_names = ['H₂', 'CO₂', 'CO', 'CH₃OH', 'H₂O', 'CH₄', 'N₂']

T_exo_0 = 503.0        # K          Inlet temperature (tube side)
P_exo_0 = 76.98e5      # Pa         Inlet pressure (76.98 bar)
F_exo_total_0 = 0.64   # mol/s      Total molar flow rate per tube

# Initial molar flow rates (mol/s)
F_exo_0 = y_exo_0 * F_exo_total_0

# --- Shell Side (Cyclohexane Dehydrogenation) ---
# Components: [C6H12 (cyclohexane), C6H6 (benzene), H2]
y_endo_0 = np.array([1.0, 0.0, 0.0])  # Pure cyclohexane feed
comp_endo_names = ['C₆H₁₂', 'C₆H₆', 'H₂']

T_endo_0 = 423.0       # K          Inlet temperature (shell side)
P_endo_0 = 1.01325e5   # Pa         Inlet pressure (1.01325 bar, atmospheric)
F_endo_total_0 = 0.111 # mol/s      Total molar flow rate per tube (Khademi et al. 2009)

# Initial molar flow rates (mol/s)
F_endo_0 = y_endo_0 * F_endo_total_0


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: THERMODYNAMIC & KINETIC DATA
# ═══════════════════════════════════════════════════════════════════════════════

# --- Cp values (simplified, J/(mol·K)) ---
Cp_dict_exo = {
    'H2': 29.0, 'CO2': 37.1, 'CO': 29.1,
    'CH3OH': 44.0, 'H2O': 33.6, 'CH4': 35.3, 'N2': 29.1
}
Cp_dict_endo = {
    'C6H12': 156.0, 'C6H6': 82.0, 'H2': 29.0
}

# Average molecular weights (g/mol)
MW_exo = np.array([2.016, 44.01, 28.01, 32.04, 18.015, 16.04, 28.013])
MW_endo = np.array([84.16, 78.11, 2.016])

# Heats of reaction (J/mol)
dH_methanol_CO2 = -49500.0     # CO2 + 3H2 -> CH3OH + H2O
dH_RWGS = 41200.0              # CO2 + H2 -> CO + H2O
dH_cyclohexane = 205400.0      # C6H12 -> C6H6 + 3H2

# Gas viscosity (Pa·s) - approximate
mu_exo = 1.5e-5    # Pa·s
mu_endo = 1.0e-5   # Pa·s


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: KINETIC MODELS
# ═══════════════════════════════════════════════════════════════════════════════

def methanol_kinetics(T, P_total, y):
    """
    Methanol synthesis kinetics (Vanden Bussche & Froment, 1996)
    
    Reactions:
      R1: CO2 + 3H2 <-> CH3OH + H2O   (Methanol from CO2)
      R2: CO2 + H2  <-> CO + H2O       (Reverse Water-Gas Shift)
    
    Parameters
    ----------
    T : float
        Temperature (K)
    P_total : float
        Total pressure (Pa)
    y : array
        Mole fractions [H2, CO2, CO, CH3OH, H2O, CH4, N2]
    
    Returns
    -------
    r1, r2 : float
        Reaction rates (mol/(kg_cat·s))
    """
    # Partial pressures (bar) — kinetics use bar
    P_bar = P_total / 1e5
    p = y * P_bar  # partial pressures in bar
    
    pH2   = max(p[0], 1e-10)
    pCO2  = max(p[1], 1e-10)
    pCO   = max(p[2], 1e-10)
    pMeOH = max(p[3], 1e-10)
    pH2O  = max(p[4], 1e-10)
    
    # --- Rate constants (Arrhenius form: k = A * exp(B/RT)) ---
    # Vanden Bussche & Froment (1996)
    k1  = 1.07   * np.exp(36696.0 / (R_gas * T))      # Methanol synthesis
    k2  = 1.22e10 * np.exp(-94765.0 / (R_gas * T))     # RWGS
    
    # --- Adsorption constants ---
    sqrt_KH2 = 0.499 * np.exp(17197.0 / (R_gas * T))   # sqrt(KH2) parameter (positive exponent: exothermic adsorption)
    KH2O = 6.62e-11 * np.exp(124119.0 / (R_gas * T))  # KH2O adsorption parameter
    
    # --- Equilibrium constants ---
    # Methanol synthesis equilibrium: CO2 + 3H2 <-> CH3OH + H2O
    # log10(Keq1) correlation (bar^-2)
    Keq1 = 10.0**(3066.0 / T - 10.592)
    
    # RWGS equilibrium: CO2 + H2 <-> CO + H2O (Endothermic: Keq2 increases with T)
    Keq2 = 10.0**(-2073.0 / T + 2.029)
    
    # --- Driving force terms ---
    # Methanol synthesis
    driving1 = pCO2 * pH2 - (1.0 / Keq1) * (pMeOH * pH2O / (pH2**2))
    
    # RWGS
    driving2 = pCO2 - (1.0 / Keq2) * (pCO * pH2O / pH2)
    
    # --- Exact Vanden Bussche & Froment (1996) Denominator ---
    denom = 1.0 + 3453.38 * (pH2O / pH2) + sqrt_KH2 * np.sqrt(pH2) + KH2O * pH2O
    
    # Catalyst internal diffusion effectiveness factor (5.47 mm pellets)
    eta_exo = 0.35
    
    # --- Reaction rates (mol / kg_cat / s) ---
    r1 = eta_exo * k1 * driving1 / (denom**3)   # Methanol synthesis from CO2
    r2 = eta_exo * k2 * driving2 / denom       # RWGS
    
    return r1, r2


def cyclohexane_kinetics(T, P_total, y):
    """
    Cyclohexane dehydrogenation kinetics (Temkin/Langmuir reversible model)
    Based on Khademi et al. (2009) and Sinfelt (1960)
    
    Reaction: C6H12 <-> C6H6 + 3H2
    
    Parameters
    ----------
    T : float
        Temperature (K)
    P_total : float
        Total pressure (Pa)
    y : array
        Mole fractions [C6H12, C6H6, H2]
    
    Returns
    -------
    r_CH : float
        Rate of cyclohexane consumption (mol/(kg_cat·s)) - fully reversible
    """
    # Partial pressures (bar)
    P_bar = P_total / 1e5
    p = y * P_bar
    pCH = max(p[0], 1e-10)   # Cyclohexane
    pBz = max(p[1], 1e-10)   # Benzene
    pH2 = max(p[2], 1e-10)   # Hydrogen
    
    # Rate constant (Arrhenius): mol/(kg_cat·s·bar)
    k_CH = 2.45e4 * np.exp(-68500.0 / (R_gas * T))
    
    # Thermodynamic equilibrium constant Kp (bar³)
    Kp = 10.0**(-18900.0 / T + 30.55)
    
    # Benzene adsorption constant (bar^-1)
    K_Bz = 0.85 * np.exp(12000.0 / (R_gas * T))
    
    # Reversible driving force and Temkin denominator (no artificial clamping)
    df_CH = pCH - (pBz * (pH2**3)) / max(Kp, 1e-10)
    denom_CH = 1.0 + K_Bz * pBz
    r_CH = k_CH * df_CH / denom_CH
    
    return r_CH


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: ERGUN EQUATION (PRESSURE DROP)
# ═══════════════════════════════════════════════════════════════════════════════

def ergun_pressure_drop(u_s, rho_g, mu_g, dp, eps):
    """
    Ergun equation for pressure drop in packed bed.
    
    Returns dP/dz (Pa/m) — negative value (pressure decreases)
    """
    term1 = 150.0 * mu_g * (1.0 - eps)**2 / (eps**3 * dp**2) * u_s
    term2 = 1.75 * (1.0 - eps) / (eps**3 * dp) * rho_g * u_s**2
    return -(term1 + term2)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: ODE SYSTEM (CO-CURRENT FLOW)
# ═══════════════════════════════════════════════════════════════════════════════

def coupled_reactor_odes(z, state):
    """
    System of ODEs for co-current thermally coupled reactor.
    
    State vector:
      [F_H2, F_CO2, F_CO, F_CH3OH, F_H2O, F_CH4, F_N2,  (tube: 0-6)
       F_C6H12, F_C6H6, F_H2_endo,                        (shell: 7-9)
       T_exo, T_endo,                                       (temperatures: 10-11)
       P_exo, P_endo]                                        (pressures: 12-13)
    """
    # ── Unpack state variables ──
    F_exo = state[0:7]     # Molar flows tube side (mol/s)
    F_endo = state[7:10]   # Molar flows shell side (mol/s)
    T_exo = state[10]      # Temperature tube (K)
    T_endo = state[11]     # Temperature shell (K)
    P_exo = state[12]      # Pressure tube (Pa)
    P_endo = state[13]     # Pressure shell (Pa)
    
    # ── Ensure non-negative flows and reasonable conditions ──
    F_exo = np.maximum(F_exo, 1e-15)
    F_endo = np.maximum(F_endo, 1e-15)
    T_exo = max(T_exo, 300.0)
    T_endo = max(T_endo, 300.0)
    P_exo = max(P_exo, 1e4)
    P_endo = max(P_endo, 1e4)
    
    # ── Mole fractions ──
    F_exo_total = np.sum(F_exo)
    F_endo_total = np.sum(F_endo)
    y_exo = F_exo / F_exo_total
    y_endo = F_endo / F_endo_total
    
    # ── Kinetics ──
    # Tube side: methanol synthesis
    r1, r2 = methanol_kinetics(T_exo, P_exo, y_exo)
    
    # Shell side: cyclohexane dehydrogenation
    r_CH = cyclohexane_kinetics(T_endo, P_endo, y_endo)
    
    # ── Material balances (tube side) ──
    # R1: CO2 + 3H2 -> CH3OH + H2O
    # R2: CO2 + H2  -> CO + H2O  (RWGS)
    dF_exo = np.zeros(7)
    dF_exo[0] = A_tube * rho_b_exo * (-3.0 * r1 - r2)     # H2
    dF_exo[1] = A_tube * rho_b_exo * (-r1 - r2)            # CO2
    dF_exo[2] = A_tube * rho_b_exo * (r2)                  # CO
    dF_exo[3] = A_tube * rho_b_exo * (r1)                  # CH3OH
    dF_exo[4] = A_tube * rho_b_exo * (r1 + r2)             # H2O
    dF_exo[5] = 0.0                                         # CH4 (inert)
    dF_exo[6] = 0.0                                         # N2  (inert)
    
    # ── Material balances (shell side) ──
    # C6H12 -> C6H6 + 3H2
    dF_endo = np.zeros(3)
    dF_endo[0] = A_shell * rho_b_endo * (-r_CH)             # C6H12
    dF_endo[1] = A_shell * rho_b_endo * (r_CH)              # C6H6
    dF_endo[2] = A_shell * rho_b_endo * (3.0 * r_CH)        # H2
    
    # ── Heat capacities (flow-weighted) ──
    Cp_exo_arr = np.array([29.0, 37.1, 29.1, 44.0, 33.6, 35.3, 29.1])
    Cp_endo_arr = np.array([156.0, 82.0, 29.0])
    
    sum_FCp_exo = np.dot(F_exo, Cp_exo_arr)    # mol/s * J/(mol·K) = W/K
    sum_FCp_endo = np.dot(F_endo, Cp_endo_arr)
    
    # ── Heat transfer between tube and shell ──
    # Q = U * a_H * (T_exo - T_endo)   per unit length
    # a_H = π * D_tube_in (inner perimeter per tube)
    a_H = np.pi * D_tube_in  # m²/m (heat transfer area per unit length)
    Q_per_length = U_overall * a_H * (T_exo - T_endo)  # W/m
    
    # ── Energy balances ──
    # Tube side (exothermic): gains heat from reaction, loses heat to shell
    heat_rxn_exo = A_tube * rho_b_exo * ((-dH_methanol_CO2) * r1 + (-dH_RWGS) * r2)
    dT_exo = (heat_rxn_exo - Q_per_length) / sum_FCp_exo
    
    # Shell side (endothermic): absorbs heat from reaction, gains heat from tube
    heat_rxn_endo = A_shell * rho_b_endo * (-dH_cyclohexane) * r_CH
    dT_endo = (heat_rxn_endo + Q_per_length) / sum_FCp_endo
    
    # ── Pressure drop (Ergun equation) ──
    # Tube side
    MW_avg_exo = np.dot(y_exo, MW_exo) / 1000.0  # kg/mol
    rho_g_exo = P_exo * MW_avg_exo / (R_gas * T_exo)  # kg/m³
    u_s_exo = F_exo_total * R_gas * T_exo / (P_exo * A_tube)  # Superficial velocity
    dP_exo = ergun_pressure_drop(u_s_exo, rho_g_exo, mu_exo, dp_exo, eps_exo)
    
    # Shell side
    MW_avg_endo = np.dot(y_endo, MW_endo) / 1000.0  # kg/mol
    rho_g_endo = P_endo * MW_avg_endo / (R_gas * T_endo)  # kg/m³
    u_s_endo = F_endo_total * R_gas * T_endo / (P_endo * A_shell_flow)
    dP_endo = ergun_pressure_drop(u_s_endo, rho_g_endo, mu_endo, dp_endo, eps_endo)
    
    # ── Assemble derivatives ──
    dstate = np.zeros(14)
    dstate[0:7] = dF_exo
    dstate[7:10] = dF_endo
    dstate[10] = dT_exo
    dstate[11] = dT_endo
    dstate[12] = dP_exo
    dstate[13] = dP_endo
    
    return dstate


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: SOLVE AND PLOT
# ═══════════════════════════════════════════════════════════════════════════════

def solve_cocurrent():
    """Solve the co-current thermally coupled reactor model."""
    
    # Initial state vector
    state0 = np.zeros(14)
    state0[0:7] = F_exo_0          # Tube molar flows
    state0[7:10] = F_endo_0        # Shell molar flows
    state0[10] = T_exo_0           # Tube temperature
    state0[11] = T_endo_0          # Shell temperature
    state0[12] = P_exo_0           # Tube pressure
    state0[13] = P_endo_0          # Shell pressure
    
    # Solve IVP
    z_span = (0.0, L)
    z_eval = np.linspace(0.0, L, 500)
    
    print("=" * 70)
    print("  THERMALLY COUPLED REACTOR - CO-CURRENT SIMULATION")
    print("=" * 70)
    print(f"\n  Reactor Length: {L:.3f} m")
    print(f"  Tube ID/OD: {D_tube_in*1000:.1f}/{D_tube_out*1000:.1f} mm")
    print(f"  Number of Tubes: {N_tubes}")
    print(f"\n  --- Tube Side (Methanol Synthesis) ---")
    print(f"  T_inlet = {T_exo_0:.1f} K ({T_exo_0-273.15:.1f} °C)")
    print(f"  P_inlet = {P_exo_0/1e5:.2f} bar")
    print(f"  F_total = {F_exo_total_0:.4f} mol/s per tube")
    print(f"\n  --- Shell Side (Cyclohexane Dehydrogenation) ---")
    print(f"  T_inlet = {T_endo_0:.1f} K ({T_endo_0-273.15:.1f} °C)")
    print(f"  P_inlet = {P_endo_0/1e5:.2f} bar")
    print(f"  F_total = {F_endo_total_0:.4f} mol/s per tube")
    print(f"\n  Solving ODEs...")
    
    sol = solve_ivp(
        coupled_reactor_odes,
        z_span,
        state0,
        method='Radau',      # Stiff solver for coupled kinetics
        t_eval=z_eval,
        rtol=1e-8,
        atol=1e-10,
        max_step=0.05
    )
    
    if not sol.success:
        print(f"\n  ⚠ Solver warning: {sol.message}")
    else:
        print("  ✓ Solution converged successfully!")
    
    return sol


def compute_results(sol):
    """Compute conversions, yields, and other results from the solution."""
    
    z = sol.t
    F_exo = sol.y[0:7, :]
    F_endo = sol.y[7:10, :]
    T_exo = sol.y[10, :]
    T_endo = sol.y[11, :]
    P_exo = sol.y[12, :]
    P_endo = sol.y[13, :]
    
    # --- Conversions & Yields ---
    # CO2 conversion in methanol synthesis
    CO2_in = F_exo[1, 0]
    CO2_out = F_exo[1, :]
    X_CO2 = (CO2_in - CO2_out) / CO2_in * 100.0
    
    # CO conversion
    CO_in = F_exo[2, 0]
    CO_out = F_exo[2, :]
    X_CO = (CO_in - CO_out) / CO_in * 100.0
    
    # Total carbon oxide conversion (CO + CO2)
    X_COx = ((CO2_in + CO_in) - (CO2_out + CO_out)) / (CO2_in + CO_in) * 100.0
    
    # Methanol yield based on total carbon oxides fed (Khademi et al. 2009, Eq. 28)
    MeOH_in = F_exo[3, 0]
    MeOH_out = F_exo[3, :]
    Y_MeOH = (MeOH_out - MeOH_in) / (CO2_in + CO_in) * 100.0
    
    # Cyclohexane conversion
    X_CH = (F_endo[0, 0] - F_endo[0, :]) / F_endo[0, 0] * 100.0
    
    # Hydrogen production from dehydrogenation
    F_H2_produced = F_endo[2, :] - F_endo[2, 0]
    
    # --- Print final results ---
    print(f"\n{'=' * 70}")
    print(f"  SIMULATION RESULTS (at z = {L:.2f} m)")
    print(f"{'=' * 70}")
    print(f"\n  --- Tube Side (Methanol Synthesis) ---")
    print(f"  CO₂ Conversion:     {X_CO2[-1]:.2f}%")
    print(f"  CO Conversion:      {X_CO[-1]:.2f}%")
    print(f"  Total COx Conv:     {X_COx[-1]:.2f}%")
    print(f"  Methanol Yield:     {Y_MeOH[-1]:.2f}% (Eq. 28: based on CO+CO₂)")
    print(f"  T_outlet:           {T_exo[-1]:.1f} K ({T_exo[-1]-273.15:.1f} °C)")
    print(f"  P_outlet:           {P_exo[-1]/1e5:.2f} bar")
    print(f"  F_CH3OH (outlet):   {F_exo[3,-1]*1000:.3f} mmol/s")
    
    print(f"\n  --- Shell Side (Cyclohexane Dehydrogenation) ---")
    print(f"  C₆H₁₂ Conversion:  {X_CH[-1]:.2f}%")
    print(f"  T_outlet:           {T_endo[-1]:.1f} K ({T_endo[-1]-273.15:.1f} °C)")
    print(f"  P_outlet:           {P_endo[-1]/1e5:.4f} bar")
    print(f"  H₂ produced:        {F_H2_produced[-1]*1000:.3f} mmol/s")
    
    print(f"\n  --- Temperature Change ---")
    print(f"  ΔT_tube:  {T_exo[-1]-T_exo[0]:+.1f} K")
    print(f"  ΔT_shell: {T_endo[-1]-T_endo[0]:+.1f} K")
    
    return {
        'z': z,
        'F_exo': F_exo, 'F_endo': F_endo,
        'T_exo': T_exo, 'T_endo': T_endo,
        'P_exo': P_exo, 'P_endo': P_endo,
        'X_CO2': X_CO2, 'X_CO': X_CO, 'X_COx': X_COx, 'Y_MeOH': Y_MeOH,
        'X_CH': X_CH, 'F_H2_produced': F_H2_produced
    }


def plot_results(results, save_path='reactor_results.png'):
    """Generate publication-quality plots of reactor performance."""
    
    z = results['z']
    
    # --- Setup figure ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Thermally Coupled Shell-and-Tube Reactor\n'
                 'Methanol Synthesis (Tube) + Cyclohexane Dehydrogenation (Shell)',
                 fontsize=14, fontweight='bold')
    
    colors = {
        'exo': '#E63946',    # Red for exothermic
        'endo': '#457B9D',   # Blue for endothermic
        'meoh': '#2A9D8F',   # Teal for methanol
        'h2': '#F4A261',     # Orange for H2
        'co2': '#264653',    # Dark for CO2
    }
    
    # --- Plot 1: Temperature Profiles ---
    ax = axes[0, 0]
    ax.plot(z, results['T_exo'], '-', color=colors['exo'], linewidth=2.5,
            label='Tube (Methanol Synthesis)')
    ax.plot(z, results['T_endo'], '-', color=colors['endo'], linewidth=2.5,
            label='Shell (C₆H₁₂ Dehydrogenation)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Temperature (K)', fontsize=11)
    ax.set_title('(a) Temperature Profiles', fontsize=12)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 2: Conversion Profiles ---
    ax = axes[0, 1]
    ax.plot(z, results['X_CO2'], '-', color=colors['co2'], linewidth=2.5,
            label='CO₂ Conversion (Tube)')
    ax.plot(z, results['X_CH'], '--', color=colors['endo'], linewidth=2.5,
            label='C₆H₁₂ Conversion (Shell)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(b) Conversion Profiles', fontsize=12)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 3: Methanol & Benzene Production ---
    ax = axes[0, 2]
    ax.plot(z, results['F_exo'][3, :] * 1000, '-', color=colors['meoh'], linewidth=2.5,
            label='CH₃OH (Tube)')
    ax2 = ax.twinx()
    ax2.plot(z, results['F_endo'][1, :] * 1000, '--', color=colors['endo'], linewidth=2.5,
             label='C₆H₆ (Shell)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('F_CH₃OH (mmol/s)', fontsize=11, color=colors['meoh'])
    ax2.set_ylabel('F_C₆H₆ (mmol/s)', fontsize=11, color=colors['endo'])
    ax.set_title('(c) Product Molar Flow Rates', fontsize=12)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='center right')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 4: Pressure Drop ---
    ax = axes[1, 0]
    ax.plot(z, results['P_exo'] / 1e5, '-', color=colors['exo'], linewidth=2.5,
            label='Tube Pressure')
    ax3 = ax.twinx()
    ax3.plot(z, results['P_endo'] / 1e5, '--', color=colors['endo'], linewidth=2.5,
             label='Shell Pressure')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('P_tube (bar)', fontsize=11, color=colors['exo'])
    ax3.set_ylabel('P_shell (bar)', fontsize=11, color=colors['endo'])
    ax.set_title('(d) Pressure Profiles', fontsize=12)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax3.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='center right')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 5: Tube-side Molar Flows ---
    ax = axes[1, 1]
    labels_plot = ['H₂', 'CO₂', 'CO', 'CH₃OH', 'H₂O']
    colors_flow = ['#E63946', '#264653', '#457B9D', '#2A9D8F', '#F4A261']
    for i in range(5):
        ax.plot(z, results['F_exo'][i, :] * 1000, '-', color=colors_flow[i],
                linewidth=2, label=labels_plot[i])
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Molar Flow (mmol/s)', fontsize=11)
    ax.set_title('(e) Tube-Side Molar Flows', fontsize=12)
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.grid(True, alpha=0.3)
    
    # --- Plot 6: H2 Production (Shell) ---
    ax = axes[1, 2]
    ax.plot(z, results['F_H2_produced'] * 1000, '-', color=colors['h2'], linewidth=2.5)
    ax.fill_between(z, 0, results['F_H2_produced'] * 1000, alpha=0.2, color=colors['h2'])
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('H₂ Produced (mmol/s)', fontsize=11)
    ax.set_title('(f) Hydrogen Production (Shell Side)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  ✓ Plot saved to: {save_path}")
    return save_path


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8: SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def sensitivity_analysis():
    """
    Perform parametric sensitivity analysis on key operating variables.
    Variables: Inlet temperature (tube), Feed ratio (shell), U_overall
    """
    
    print(f"\n{'=' * 70}")
    print("  PARAMETRIC SENSITIVITY ANALYSIS")
    print(f"{'=' * 70}")
    
    global T_exo_0, F_endo_total_0, U_overall, F_endo_0
    
    # Save original values
    T_exo_orig = T_exo_0
    F_endo_orig = F_endo_total_0
    U_orig = U_overall
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Sensitivity Analysis — Thermally Coupled Reactor',
                 fontsize=14, fontweight='bold')
    
    # --- Effect of Tube Inlet Temperature ---
    T_range = np.arange(490, 535, 5)
    X_CO2_final = []
    X_CH_final = []
    
    print("\n  [1/3] Varying tube inlet temperature...")
    for T in T_range:
        T_exo_0 = T
        sol = solve_ivp(coupled_reactor_odes, (0, L),
                       np.concatenate([F_exo_0, F_endo_0, [T_exo_0, T_endo_0, P_exo_0, P_endo_0]]),
                       method='Radau', t_eval=[L], rtol=1e-6, atol=1e-8, max_step=0.1)
        if sol.success:
            X_CO2_final.append((F_exo_0[1] - sol.y[1, -1]) / F_exo_0[1] * 100)
            X_CH_final.append((F_endo_0[0] - sol.y[7, -1]) / F_endo_0[0] * 100)
        else:
            X_CO2_final.append(np.nan)
            X_CH_final.append(np.nan)
    
    T_exo_0 = T_exo_orig  # Reset
    
    ax = axes[0]
    ax.plot(T_range, X_CO2_final, 'o-', color='#E63946', linewidth=2, label='CO₂ Conv.')
    ax.plot(T_range, X_CH_final, 's--', color='#457B9D', linewidth=2, label='C₆H₁₂ Conv.')
    ax.set_xlabel('Tube Inlet Temperature (K)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(a) Effect of Tube Inlet Temperature', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # --- Effect of Shell Feed Flow Rate ---
    F_range = np.arange(0.02, 0.22, 0.02)
    X_CO2_final2 = []
    X_CH_final2 = []
    
    print("  [2/3] Varying shell feed flow rate...")
    for F_val in F_range:
        F_endo_total_0 = F_val
        F_endo_0_temp = np.array([F_val, 0.0, 0.0])
        sol = solve_ivp(coupled_reactor_odes, (0, L),
                       np.concatenate([F_exo_0, F_endo_0_temp, [T_exo_orig, T_endo_0, P_exo_0, P_endo_0]]),
                       method='Radau', t_eval=[L], rtol=1e-6, atol=1e-8, max_step=0.1)
        if sol.success:
            X_CO2_final2.append((F_exo_0[1] - sol.y[1, -1]) / F_exo_0[1] * 100)
            X_CH_final2.append((F_endo_0_temp[0] - sol.y[7, -1]) / F_endo_0_temp[0] * 100)
        else:
            X_CO2_final2.append(np.nan)
            X_CH_final2.append(np.nan)
    
    F_endo_total_0 = F_endo_orig  # Reset
    F_endo_0 = y_endo_0 * F_endo_total_0
    
    ax = axes[1]
    ax.plot(F_range * 1000, X_CO2_final2, 'o-', color='#E63946', linewidth=2, label='CO₂ Conv.')
    ax.plot(F_range * 1000, X_CH_final2, 's--', color='#457B9D', linewidth=2, label='C₆H₁₂ Conv.')
    ax.set_xlabel('Shell Feed Flow (mmol/s)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(b) Effect of Shell Feed Flow Rate', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # --- Effect of Overall Heat Transfer Coefficient ---
    U_range = np.arange(300, 950, 60)
    X_CO2_final3 = []
    X_CH_final3 = []
    T_exo_out3 = []
    
    print("  [3/3] Varying heat transfer coefficient...")
    for U_val in U_range:
        U_overall = U_val
        sol = solve_ivp(coupled_reactor_odes, (0, L),
                       np.concatenate([F_exo_0, F_endo_0, [T_exo_orig, T_endo_0, P_exo_0, P_endo_0]]),
                       method='Radau', t_eval=[L], rtol=1e-6, atol=1e-8, max_step=0.1)
        if sol.success:
            X_CO2_final3.append((F_exo_0[1] - sol.y[1, -1]) / F_exo_0[1] * 100)
            X_CH_final3.append((F_endo_0[0] - sol.y[7, -1]) / F_endo_0[0] * 100)
            T_exo_out3.append(sol.y[10, -1])
        else:
            X_CO2_final3.append(np.nan)
            X_CH_final3.append(np.nan)
            T_exo_out3.append(np.nan)
    
    U_overall = U_orig  # Reset
    
    ax = axes[2]
    ax.plot(U_range, X_CO2_final3, 'o-', color='#E63946', linewidth=2, label='CO₂ Conv.')
    ax.plot(U_range, X_CH_final3, 's--', color='#457B9D', linewidth=2, label='C₆H₁₂ Conv.')
    ax.set_xlabel('U (W/m²·K)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(c) Effect of Heat Transfer Coefficient', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sensitivity_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("\n  ✓ Sensitivity analysis plot saved to: sensitivity_analysis.png")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9: COMPARISON — COUPLED vs. CONVENTIONAL REACTOR
# ═══════════════════════════════════════════════════════════════════════════════

def conventional_reactor_odes(z, state):
    """
    Conventional methanol reactor with constant-temperature coolant (no coupling).
    
    State: [F_H2, F_CO2, F_CO, F_CH3OH, F_H2O, F_CH4, F_N2, T_exo, P_exo]
    """
    F_exo = state[0:7]
    T_exo = state[7]
    P_exo = state[8]
    
    F_exo = np.maximum(F_exo, 1e-15)
    T_exo = max(T_exo, 300.0)
    P_exo = max(P_exo, 1e4)
    
    F_exo_total = np.sum(F_exo)
    y_exo = F_exo / F_exo_total
    
    r1, r2 = methanol_kinetics(T_exo, P_exo, y_exo)
    
    dF_exo = np.zeros(7)
    dF_exo[0] = A_tube * rho_b_exo * (-3.0 * r1 - r2)
    dF_exo[1] = A_tube * rho_b_exo * (-r1 - r2)
    dF_exo[2] = A_tube * rho_b_exo * (r2)
    dF_exo[3] = A_tube * rho_b_exo * (r1)
    dF_exo[4] = A_tube * rho_b_exo * (r1 + r2)
    dF_exo[5] = 0.0
    dF_exo[6] = 0.0
    
    Cp_exo_arr = np.array([29.0, 37.1, 29.1, 44.0, 33.6, 35.3, 29.1])
    sum_FCp = np.dot(F_exo, Cp_exo_arr)
    
    T_coolant = 503.0  # K — boiling water coolant
    a_H = np.pi * D_tube_in
    Q_per_length = U_overall * a_H * (T_exo - T_coolant)
    
    heat_rxn = A_tube * rho_b_exo * ((-dH_methanol_CO2) * r1 + (-dH_RWGS) * r2)
    dT = (heat_rxn - Q_per_length) / sum_FCp
    
    MW_avg = np.dot(y_exo, MW_exo) / 1000.0
    rho_g = P_exo * MW_avg / (R_gas * T_exo)
    u_s = F_exo_total * R_gas * T_exo / (P_exo * A_tube)
    dP = ergun_pressure_drop(u_s, rho_g, mu_exo, dp_exo, eps_exo)
    
    dstate = np.zeros(9)
    dstate[0:7] = dF_exo
    dstate[7] = dT
    dstate[8] = dP
    
    return dstate


def compare_configurations():
    """Compare thermally coupled reactor with conventional reactor."""
    
    print(f"\n{'=' * 70}")
    print("  COMPARISON: COUPLED vs. CONVENTIONAL REACTOR")
    print(f"{'=' * 70}")
    
    # Solve coupled reactor
    state0_coupled = np.concatenate([F_exo_0, F_endo_0,
                                     [T_exo_0, T_endo_0, P_exo_0, P_endo_0]])
    z_eval = np.linspace(0, L, 500)
    
    sol_coupled = solve_ivp(coupled_reactor_odes, (0, L), state0_coupled,
                           method='Radau', t_eval=z_eval, rtol=1e-8, atol=1e-10,
                           max_step=0.05)
    
    # Solve conventional reactor
    state0_conv = np.concatenate([F_exo_0, [T_exo_0, P_exo_0]])
    sol_conv = solve_ivp(conventional_reactor_odes, (0, L), state0_conv,
                        method='Radau', t_eval=z_eval, rtol=1e-8, atol=1e-10,
                        max_step=0.05)
    
    # Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Coupled vs. Conventional Methanol Reactor Comparison',
                 fontsize=14, fontweight='bold')
    
    # Temperature comparison
    ax = axes[0]
    ax.plot(sol_coupled.t, sol_coupled.y[10, :], '-', color='#E63946', linewidth=2.5,
            label='Coupled (Tube)')
    ax.plot(sol_coupled.t, sol_coupled.y[11, :], '--', color='#457B9D', linewidth=2.5,
            label='Coupled (Shell)')
    ax.plot(sol_conv.t, sol_conv.y[7, :], '-.', color='#2A9D8F', linewidth=2.5,
            label='Conventional (with coolant)')
    ax.axhline(y=503, color='gray', linestyle=':', alpha=0.5, label='Coolant T')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Temperature (K)', fontsize=11)
    ax.set_title('(a) Temperature Profiles', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # CO2 conversion comparison
    ax = axes[1]
    X_coupled = (F_exo_0[1] - sol_coupled.y[1, :]) / F_exo_0[1] * 100
    X_conv = (F_exo_0[1] - sol_conv.y[1, :]) / F_exo_0[1] * 100
    ax.plot(sol_coupled.t, X_coupled, '-', color='#E63946', linewidth=2.5,
            label='Coupled Reactor')
    ax.plot(sol_conv.t, X_conv, '--', color='#2A9D8F', linewidth=2.5,
            label='Conventional Reactor')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('CO₂ Conversion (%)', fontsize=11)
    ax.set_title('(b) CO₂ Conversion Comparison', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('comparison_coupled_vs_conventional.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n  Coupled Reactor — CO₂ Conv: {X_coupled[-1]:.2f}%, T_out: {sol_coupled.y[10,-1]:.1f} K")
    print(f"  Conventional   — CO₂ Conv: {X_conv[-1]:.2f}%, T_out: {sol_conv.y[7,-1]:.1f} K")
    print("\n  ✓ Comparison plot saved to: comparison_coupled_vs_conventional.png")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 10: MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    
    # 1. Solve co-current model
    sol = solve_cocurrent()
    
    # 2. Compute and display results
    results = compute_results(sol)
    
    # 3. Generate main plots
    plot_results(results)
    
    # 4. Run sensitivity analysis
    sensitivity_analysis()
    
    # 5. Compare with conventional reactor
    compare_configurations()
    
    print(f"\n{'=' * 70}")
    print("  ALL SIMULATIONS COMPLETE ✓")
    print(f"{'=' * 70}")
    print("\n  Output files:")
    print("    - reactor_results.png")
    print("    - sensitivity_analysis.png")
    print("    - comparison_coupled_vs_conventional.png")

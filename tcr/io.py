"""
Input/Output and serialization module:
- Exports master simulation_results.json
- Exports simulation_summary.md
- Computes process thermal duties (preheating, vaporization) from thermodynamics
"""

import os
import json
import time
import numpy as np
from .thermo import calc_cp_species, MW
from .geometry import ReactorGeometry
from .config import AppConfig

def calc_process_thermal_duties(cfg: AppConfig, geom: ReactorGeometry, res_knee: dict) -> dict:
    """
    Computes rigorous thermal duties for feed preheating and vaporization (MW).
    Replaces arbitrary hardcoded numbers with physical thermodynamic integration.
    """
    T_ambient = 298.15  # K (25 °C)

    # 1. Syngas feed preheating: ambient (298.15 K) to T_exo_0
    T_tube_in = float(res_knee.get('decision_variables', {}).get('T_exo_0_K', cfg.feed_tube.T0))
    F_tube_tot = float(res_knee.get('decision_variables', {}).get('F_exo_mol_s', cfg.feed_tube.F_total0))
    y_tube = cfg.feed_tube.y0
    species_tube = ['H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2']

    preheat_tube_W_per_tube = 0.0
    for idx, sp in enumerate(species_tube):
        F_sp = F_tube_tot * y_tube[idx]
        # trapezoidal numerical integration of Cp(T) from T_ambient to T_tube_in
        T_arr = np.linspace(T_ambient, T_tube_in, 50)
        Cp_arr = np.array([calc_cp_species(sp, T) for T in T_arr])
        delta_h_sp = np.trapezoid(Cp_arr, T_arr)
        preheat_tube_W_per_tube += F_sp * delta_h_sp

    preheat_tube_MW = preheat_tube_W_per_tube * geom.N_tubes * 1e-6

    # 2. Cyclohexane liquid preheat & vaporization: ambient to normal boiling point (353.9 K) + Delta H_vap (30.1 kJ/mol) + vapor preheat to T_endo_0
    T_shell_in = float(res_knee.get('decision_variables', {}).get('T_endo_0_K', cfg.feed_shell.T0))
    F_shell_tot = float(res_knee.get('decision_variables', {}).get('F_endo_mol_s', cfg.feed_shell.F_total0))

    # Cyclohexane physical properties:
    T_boil_ch = 353.93  # K (80.78 °C, NIST)
    delta_H_vap_ch = 30100.0  # J/mol at normal boiling point (NIST)
    Cp_liq_ch = 156.0  # J/(mol·K) liquid mean heat capacity

    # Liquid preheat: 298.15 -> 353.93 K
    q_liq_preheat = Cp_liq_ch * (T_boil_ch - T_ambient)
    # Vaporization:
    q_vap = delta_H_vap_ch
    # Vapor superheat: 353.93 K -> T_shell_in
    T_arr_vap = np.linspace(T_boil_ch, max(T_shell_in, T_boil_ch), 50)
    Cp_vap_arr = np.array([calc_cp_species('C6H12', T) for T in T_arr_vap])
    q_superheat = np.trapezoid(Cp_vap_arr, T_arr_vap)

    q_ch_total_J_per_mol = q_liq_preheat + q_vap + q_superheat
    preheat_shell_MW = F_shell_tot * q_ch_total_J_per_mol * geom.N_tubes * 1e-6

    return {
        "syngas_preheat_duty_MW": round(float(preheat_tube_MW), 2),
        "cyclohexane_vaporization_and_preheat_duty_MW": round(float(preheat_shell_MW), 2),
        "total_feed_preheat_duty_MW": round(float(preheat_tube_MW + preheat_shell_MW), 2),
        "reaction_zone_external_cooling_utility_kW": 0.0
    }

def export_simulation_results(data: dict, json_path: str = "simulation_results.json"):
    """Saves master dictionary to JSON."""
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

def export_simulation_summary(data: dict, md_path: str = "simulation_summary.md"):
    """Generates markdown summary table directly from simulation dictionary."""
    conv = data.get("conventional_baseline", {})
    conv_opt = data.get("conventional_optimized", {})
    co = data.get("cocurrent_baseline", {})
    ct = data.get("countercurrent_baseline", {})
    knee = data.get("tcr_nsga2_optimum", {})

    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Thermally Coupled Reactor (TCR) Simulation Performance Summary\n\n")
        f.write("Generated from Single Source of Truth architecture (`tcr/` package).\n\n")
        f.write("| Performance Indicator | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")
        f.write(f"| **Methanol Single-Pass Yield (%)** | {conv.get('yield_khademi_eq28', 0.0):.2f}% | {conv_opt.get('yield_khademi_eq28', 0.0):.2f}% | {co.get('yield_khademi_eq28', 0.0):.2f}% | {ct.get('yield_khademi_eq28', 0.0):.2f}% | **{knee.get('yield_khademi_eq28', 0.0):.2f}%** |\n")
        f.write(f"| **Methanol Production (mol/s)** | {conv.get('MeOH_rate_net_plant', 0.0):.2f} | {conv_opt.get('MeOH_rate_net_plant', 0.0):.2f} | {co.get('MeOH_rate_net_plant', 0.0):.2f} | {ct.get('MeOH_rate_net_plant', 0.0):.2f} | **{knee.get('MeOH_rate_net_plant', 0.0):.2f}** |\n")
        f.write(f"| **Methanol Production (ton/day)** | {conv.get('daily_meoh_ton', 0.0):.2f} | {conv_opt.get('daily_meoh_ton', 0.0):.2f} | {co.get('daily_meoh_ton', 0.0):.2f} | {ct.get('daily_meoh_ton', 0.0):.2f} | **{knee.get('daily_meoh_ton', 0.0):.2f}** |\n")
        f.write(f"| **CO₂ Conversion (%)** | {conv.get('X_CO2', 0.0):.2f}% | {conv_opt.get('X_CO2', 0.0):.2f}% | {co.get('X_CO2', 0.0):.2f}% | {ct.get('X_CO2', 0.0):.2f}% | **{knee.get('X_CO2', 0.0):.2f}%** |\n")
        f.write(f"| **Total COx Conversion (%)** | {conv.get('X_COx', 0.0):.2f}% | {conv_opt.get('X_COx', 0.0):.2f}% | {co.get('X_COx', 0.0):.2f}% | {ct.get('X_COx', 0.0):.2f}% | **{knee.get('X_COx', 0.0):.2f}%** |\n")
        f.write(f"| **Cyclohexane Conversion (%)** | N/A | N/A | {co.get('X_CH', 0.0):.2f}% | {ct.get('X_CH', 0.0):.2f}% | **{knee.get('X_CH', 0.0):.2f}%** |\n")
        f.write(f"| **Benzene Production (mol/s)** | 0.00 | 0.00 | {co.get('Benzene_rate_plant', 0.0):.2f} | {ct.get('Benzene_rate_plant', 0.0):.2f} | **{knee.get('Benzene_rate_plant', 0.0):.2f}** |\n")
        f.write(f"| **Benzene Production (ton/day)** | 0.00 | 0.00 | {co.get('daily_benzene_ton', 0.0):.2f} | {ct.get('daily_benzene_ton', 0.0):.2f} | **{knee.get('daily_benzene_ton', 0.0):.2f}** |\n")
        f.write(f"| **Coproduced H₂ Rate (mol/s)** | 0.00 | 0.00 | {co.get('H2_rate_plant', 0.0):.2f} | {ct.get('H2_rate_plant', 0.0):.2f} | **{knee.get('H2_rate_plant', 0.0):.2f}** |\n")
        f.write(f"| **Coproduced H₂ Rate (ton/day)** | 0.00 | 0.00 | {co.get('daily_h2_ton', 0.0):.2f} | {ct.get('daily_h2_ton', 0.0):.2f} | **{knee.get('daily_h2_ton', 0.0):.2f}** |\n")
        f.write(f"| **Peak Bed Temperature (K)** | {conv.get('T_max_exo', 0.0):.2f} | {conv_opt.get('T_max_exo', 0.0):.2f} | {co.get('T_max_exo', 0.0):.2f} | {ct.get('T_max_exo', 0.0):.2f} | **{knee.get('T_max_exo', 0.0):.2f}** |\n")
        f.write(f"| **Tube Pressure Drop (bar)** | {conv.get('dP_tube', 0.0):.2f} | {conv_opt.get('dP_tube', 0.0):.2f} | {co.get('dP_tube', 0.0):.2f} | {ct.get('dP_tube', 0.0):.2f} | **{knee.get('dP_tube', 0.0):.2f}** |\n")
        f.write(f"| **Shell Pressure Drop (bar)** | N/A | N/A | {co.get('dP_shell', 0.0):.2f} | {ct.get('dP_shell', 0.0):.2f} | **{knee.get('dP_shell', 0.0):.2f}** |\n")
        f.write(f"| **Reaction Cooling Utility Required** | Boiling Water ({conv.get('Q_wall_total_MW', 0.0):.2f} MW) | Boiling Water ({conv_opt.get('Q_wall_total_MW', 0.0):.2f} MW) | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** |\n")

"""
========================================================================================
Single Entry Point (SSOT): Thermally Coupled Shell-and-Tube Reactor (TCR) Pipeline
========================================================================================
Executes:
1. Conventional Methanol Reactor (CMR) Baseline & Optimized
2. TCR Co-Current & Counter-Current Rigorous Simulations
3. NSGA-II Multi-Objective Pareto Optimization
4. Khademi et al. (2009) Benchmark Comparison & Typo Analysis
5. Sensitivity Analysis
6. High-Resolution Publication Figure Generation (300 DPI)
7. Single Source of Truth Serialization to simulation_results.json & simulation_summary.md
8. Automated Manuscript Compilation
========================================================================================
"""

import sys
import os
import time
import numpy as np

from tcr.config import load_config
from tcr.geometry import ReactorGeometry
from tcr.solvers import solve_conventional, solve_cocurrent, solve_countercurrent
from tcr.metrics import evaluate_run_metrics
from tcr.optimize import optimize_cmr, run_nsga2_tcr
from tcr.benchmark import run_khademi_benchmark
from tcr.plots import (
    plot_multiphysics_profiles, plot_pareto_frontier, plot_reactor_results,
    plot_comparison_coupled_vs_conventional, plot_sensitivity_analysis
)
from tcr.io import (
    calc_process_thermal_duties, export_simulation_results, export_simulation_summary
)

def main():
    print("=" * 85)
    print("  THERMALLY COUPLED REACTOR (TCR) RIGOROUS ENGINEERING SIMULATION PIPELINE")
    print("=" * 85)

    t_start_total = time.time()
    cfg = load_config("config.yaml")
    geom = ReactorGeometry(cfg.geometry)

    # ──────────────────────────────────────────────────────────────────────────
    # Step 1: Conventional Reactor Baseline (CMR)
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[1/7] Simulating Conventional Jacketed Reactor Baseline...")
    t0 = time.time()
    sol_cmr_base = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    res_cmr_base = evaluate_run_metrics(sol_cmr_base, geom, is_conventional=True)
    print(f"  ✓ CMR Baseline solved in {time.time()-t0:.2f}s | MeOH Yield: {res_cmr_base['yield_khademi_eq28']:.2f}% | "
          f"Rate: {res_cmr_base['MeOH_rate_net_plant']:.2f} mol/s ({res_cmr_base['daily_meoh_ton']:.2f} t/d) | T_max: {res_cmr_base['T_max_exo']:.1f} K")

    # ──────────────────────────────────────────────────────────────────────────
    # Step 2: Conventional Reactor Fair Optimization (CMR_opt)
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[2/7] Optimizing Conventional Reactor (CMR) over (T_exo_0, T_coolant)...")
    t0 = time.time()
    res_cmr_opt = optimize_cmr(cfg)
    print(f"  ✓ CMR Optimized solved in {time.time()-t0:.2f}s | Optimal T_in: {res_cmr_opt['T_exo_0']:.1f} K, T_cool: {res_cmr_opt['T_coolant']:.1f} K | "
          f"MeOH Yield: {res_cmr_opt['yield_khademi_eq28']:.2f}% | Rate: {res_cmr_opt['MeOH_rate_net_plant']:.2f} mol/s | T_max: {res_cmr_opt['T_max_exo']:.1f} K")

    # ──────────────────────────────────────────────────────────────────────────
    # Step 3: TCR Co-Current Simulation
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[3/7] Simulating TCR Co-Current Baseline...")
    t0 = time.time()
    sol_co_base = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    res_co_base = evaluate_run_metrics(sol_co_base, geom, mode="cocurrent")
    print(f"  ✓ TCR Co-Current solved in {time.time()-t0:.2f}s | MeOH Yield: {res_co_base['yield_khademi_eq28']:.2f}% | "
          f"CH Conv: {res_co_base['X_CH']:.2f}% | Benzene: {res_co_base['Benzene_rate_plant']:.2f} mol/s | H2: {res_co_base['H2_rate_plant']:.2f} mol/s | T_max: {res_co_base['T_max_exo']:.1f} K")

    # ──────────────────────────────────────────────────────────────────────────
    # Step 4: TCR Counter-Current Simulation (Picard Fixed-Point Iteration)
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[4/7] Simulating TCR Counter-Current BVP via Under-Relaxed Picard Iteration...")
    t0 = time.time()
    sol_ct_base, ct_info = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    res_ct_base = evaluate_run_metrics(sol_ct_base, geom, mode="countercurrent")
    print(f"  ✓ TCR Counter-Current converged in {time.time()-t0:.2f}s ({ct_info['iterations']} iters, res_T: {ct_info['final_res_T']:.2e}) | "
          f"MeOH Yield: {res_ct_base['yield_khademi_eq28']:.2f}% | CH Conv: {res_ct_base['X_CH']:.2f}% | Benzene: {res_ct_base['Benzene_rate_plant']:.2f} mol/s | H2: {res_ct_base['H2_rate_plant']:.2f} mol/s | T_max: {res_ct_base['T_max_exo']:.1f} K")

    # ──────────────────────────────────────────────────────────────────────────
    # Step 5: NSGA-II Multi-Objective Optimization
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[5/7] Executing NSGA-II Multi-Objective Pareto Optimization...")
    opt_res = run_nsga2_tcr(cfg)
    knee_results = opt_res['knee_results']
    knee_vars = knee_results['decision_variables']
    print(f"  ✓ NSGA-II completed in {opt_res['opt_time_s']:.2f}s | Knee Point: T_exo={knee_vars['T_exo_0_K']:.1f} K, T_endo={knee_vars['T_endo_0_K']:.1f} K, "
          f"F_endo={knee_vars['F_endo_mol_s']:.4f} mol/s, P_endo={knee_vars.get('P_endo_0_bar', 2.0):.2f} bar, F_exo={knee_vars['F_exo_mol_s']:.4f} mol/s | "
          f"MeOH Yield: {knee_results['yield_khademi_eq28']:.2f}% | Rate: {knee_results['MeOH_rate_net_plant']:.2f} mol/s ({knee_results['daily_meoh_ton']:.2f} t/d) | "
          f"Benzene: {knee_results['Benzene_rate_plant']:.2f} mol/s | H2: {knee_results['H2_rate_plant']:.2f} mol/s | T_max: {knee_results['T_max_exo']:.1f} K")

    # ──────────────────────────────────────────────────────────────────────────
    # Step 6: Benchmark & Sensitivity Analysis
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[6/7] Running Khademi et al. (2009) Benchmark & Parametric Sensitivity...")
    bench_data = run_khademi_benchmark(cfg)
    print(f"  ✓ Khademi Benchmark: Published Yield={bench_data['khademi_published']['tcr_yield_pct']:.2f}%, "
          f"Simulated Thermo Yield={bench_data['simulation_tcr_thermo']['yield_pct']:.2f}%")

    # Parametric sweeps
    sens_results = {
        'T_exo_sweep': {
            'T_vals': [490.0, 495.0, 500.0, 503.15, 510.0, 515.0, 520.0],
            'X_CO2': [],
            'X_CH': []
        },
        'F_endo_sweep': {
            'F_vals': [0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080],
            'X_CO2': [],
            'X_CH': []
        },
        'eta_sweep': {
            'eta_vals': [0.20, 0.35, 0.50, 0.70, 1.00],
            'X_CO2': [],
            'yield_meoh': []
        }
    }

    # Sweep 1: T_exo_0
    for t_val in sens_results['T_exo_sweep']['T_vals']:
        sol_tmp = solve_cocurrent(
            geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
            cfg.feed_tube.y0, t_val, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
            cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
            kinetics_mode="production_vbf"
        )
        res_tmp = evaluate_run_metrics(sol_tmp, geom, mode="cocurrent")
        sens_results['T_exo_sweep']['X_CO2'].append(res_tmp['X_CO2'])
        sens_results['T_exo_sweep']['X_CH'].append(res_tmp['X_CH'])

    # Sweep 2: F_endo_0
    for f_val in sens_results['F_endo_sweep']['F_vals']:
        sol_tmp = solve_cocurrent(
            geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
            cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
            cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, f_val,
            kinetics_mode="production_vbf"
        )
        res_tmp = evaluate_run_metrics(sol_tmp, geom, mode="cocurrent")
        sens_results['F_endo_sweep']['X_CO2'].append(res_tmp['X_CO2'])
        sens_results['F_endo_sweep']['X_CH'].append(res_tmp['X_CH'])

    # Sweep 3: eta_exo
    for eta_val in sens_results['eta_sweep']['eta_vals']:
        cat_tmp = CatalystBedConfig(dp=cfg.cat_tube.dp, rho_b=cfg.cat_tube.rho_b, eps=cfg.cat_tube.eps, eta=eta_val)
        sol_tmp = solve_cocurrent(
            geom, cat_tmp, cfg.cat_shell, cfg.transport, cfg.solver,
            cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
            cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
            kinetics_mode="production_vbf"
        )
        res_tmp = evaluate_run_metrics(sol_tmp, geom, mode="cocurrent")
        sens_results['eta_sweep']['X_CO2'].append(res_tmp['X_CO2'])
        sens_results['eta_sweep']['yield_meoh'].append(res_tmp['yield_khademi_eq28'])

    # ──────────────────────────────────────────────────────────────────────────
    # Step 7: Plotting & Master Artifact Export
    # ──────────────────────────────────────────────────────────────────────────
    print("\n[7/7] Generating All 5 Master Figures & Serializing SSOT Deliverables...")
    plot_multiphysics_profiles(sol_co_base, sol_ct_base, sol_cmr_base, geom, "tcr_multiphysics_profiles.png")
    plot_pareto_frontier(opt_res, "nsga2_pareto_front.png")
    plot_reactor_results(sol_co_base, geom, "reactor_results.png")
    plot_comparison_coupled_vs_conventional(sol_co_base, sol_cmr_base, "comparison_coupled_vs_conventional.png")
    plot_sensitivity_analysis(sens_results, "sensitivity_analysis.png")
    print("  ✓ All 5 PNG figures successfully generated.")

    # Process thermal duties
    process_duties = calc_process_thermal_duties(cfg, geom, knee_results)

    # Attach verified equal-throughput syngas feed rate and shell exit pressures
    res_cmr_base['F_exo_total_0'] = float(cfg.feed_tube.F_total0)
    res_cmr_opt['F_exo_total_0'] = float(cfg.feed_tube.F_total0)
    res_co_base['F_exo_total_0'] = float(cfg.feed_tube.F_total0)
    res_ct_base['F_exo_total_0'] = float(cfg.feed_tube.F_total0)
    knee_results['F_exo_total_0'] = float(cfg.feed_tube.F_total0)
    res_co_base['P_shell_out'] = round(float((cfg.feed_shell.P0 / 1e5) - res_co_base['dP_shell']), 4)
    res_ct_base['P_shell_out'] = round(float((cfg.feed_shell.P0 / 1e5) - res_ct_base['dP_shell']), 4)

    # Master SSOT data structure
    master_ssot = {
        "metadata": {
            "title": "Thermally Coupled Shell-and-Tube Reactor (TCR) Single Source of Truth",
            "version": "2.0.0",
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "execution_time_total_s": round(time.time() - t_start_total, 2),
            "N_tubes": geom.N_tubes,
            "shell_inlet_pressure_bar": cfg.feed_shell.P0 / 1e5
        },
        "conventional_baseline": res_cmr_base,
        "conventional_optimized": res_cmr_opt,
        "cocurrent_baseline": res_co_base,
        "countercurrent_baseline": res_ct_base,
        "tcr_nsga2_optimum": knee_results,
        "benchmark_khademi_2009": bench_data,
        "process_thermal_duties": process_duties,
        "sensitivity_analysis": sens_results
    }

    export_simulation_results(master_ssot, "simulation_results.json")
    export_simulation_summary(master_ssot, "simulation_summary.md")
    print("  ✓ Exported simulation_results.json and simulation_summary.md")

    # If build_manuscript.py exists, compile manuscript
    if os.path.exists("build_manuscript.py"):
        import subprocess
        subprocess.run([sys.executable, "build_manuscript.py"], check=True)
        print("  ✓ Compiled TCR_Methanol_Cyclohexane_Manuscript.md from template")

    print("\n" + "=" * 85)
    print(f"  PIPELINE COMPLETED SUCCESSFULLY IN {time.time() - t_start_total:.2f} SECONDS! ✓")
    print("=" * 85)

if __name__ == "__main__":
    from tcr.config import CatalystBedConfig
    main()

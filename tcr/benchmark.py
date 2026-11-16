"""
Benchmark module against Khademi et al. (2009).
Implements exact operating conditions (Tables 7-9) under two sub-modes:
1. 'khademi_literal': Graaf (1988) kinetics + literal printed Table 4 dehydrogenation
   with published typo Kp = 4.89e35 * exp(-3190/T) Pa^3.
2. 'khademi_thermo': Graaf (1988) kinetics + exact Kirchhoff thermodynamically consistent Kp.

Also evaluates:
- Intraparticle diffusion effectiveness factor sensitivities for tube (eta_exo) and shell (eta_endo).
- Axial Thiele modulus Phi(z) and dynamic effectiveness factor eta(z) profiles.
"""

import numpy as np
from .geometry import ReactorGeometry
from .config import AppConfig, CatalystBedConfig, SolverConfig
from .solvers import solve_conventional, solve_cocurrent
from .metrics import evaluate_run_metrics
from .kinetics import set_itoh_benchmark_mode, calc_axial_thiele_profile


def run_khademi_benchmark(cfg: AppConfig) -> dict:
    """
    Runs the Khademi et al. (2009) benchmark comparison:
    Tube: T0 = 527.8 K, P0 = 76.98 bar, F0 = 0.640 mol/s
    Shell: T0 = 423.0 K, P0 = 1.01325 bar, F0 = 0.111 mol/s (10% C6H12 + 90% Ar)
    """
    geom = ReactorGeometry(cfg.geometry)
    sol_cfg = SolverConfig(
        z_points=600, rtol=1e-6, atol=1e-8,
        picard_max_iter=100, picard_tol_T=1e-3, picard_tol_F=1e-6, picard_alpha=0.70
    )

    feed_t = cfg.feed_tube_bench
    feed_s = cfg.feed_shell_bench

    # 1. Conventional Reactor Benchmark (Graaf 1988, eta_exo = 1.0 as assumed in Khademi pseudo-homogeneous)
    cat_exo_10 = CatalystBedConfig(dp=cfg.cat_tube.dp, rho_b=cfg.cat_tube.rho_b, eps=cfg.cat_tube.eps, eta=1.0)
    sol_cmr_bench = solve_conventional(
        geom, cat_exo_10, cfg.transport, sol_cfg,
        feed_t.y0, feed_t.T0, feed_t.P0, feed_t.F_total0,
        T_coolant=503.15, kinetics_mode="benchmark_graaf"
    )
    res_cmr_bench = evaluate_run_metrics(sol_cmr_bench, geom, is_conventional=True)

    # 2. Co-Current TCR Benchmark — Sub-mode: khademi_thermo (Kirchhoff-corrected mode)
    set_itoh_benchmark_mode("khademi_thermo")
    try:
        sol_tcr_thermo = solve_cocurrent(
            geom, cat_exo_10, cfg.cat_shell, cfg.transport, sol_cfg,
            feed_t.y0, feed_t.T0, feed_t.P0, feed_t.F_total0,
            feed_s.y0, feed_s.T0, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_tcr_thermo = evaluate_run_metrics(sol_tcr_thermo, geom, mode="cocurrent")
    finally:
        set_itoh_benchmark_mode(None)

    # 3. Co-Current TCR Benchmark — Sub-mode: khademi_literal (Literal Table 4 mode with printed typo)
    set_itoh_benchmark_mode("khademi_literal")
    try:
        sol_tcr_literal = solve_cocurrent(
            geom, cat_exo_10, cfg.cat_shell, cfg.transport, sol_cfg,
            feed_t.y0, feed_t.T0, feed_t.P0, feed_t.F_total0,
            feed_s.y0, feed_s.T0, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_tcr_literal = evaluate_run_metrics(sol_tcr_literal, geom, mode="cocurrent")
    finally:
        set_itoh_benchmark_mode(None)

    # 4. Effect of tube intraparticle diffusion effectiveness factor eta_exo in {1.0, 0.7, 0.5, 0.35}
    eta_exo_sweep = {}
    for eta_val in [1.0, 0.7, 0.5, 0.35]:
        cat_var = CatalystBedConfig(dp=cfg.cat_tube.dp, rho_b=cfg.cat_tube.rho_b, eps=cfg.cat_tube.eps, eta=eta_val)
        sol_eta = solve_cocurrent(
            geom, cat_var, cfg.cat_shell, cfg.transport, sol_cfg,
            feed_t.y0, feed_t.T0, feed_t.P0, feed_t.F_total0,
            feed_s.y0, feed_s.T0, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_eta = evaluate_run_metrics(sol_eta, geom, mode="cocurrent")
        eta_exo_sweep[f"eta_{eta_val}"] = {
            "yield_khademi_eq28": res_eta["yield_khademi_eq28"],
            "X_CO2": res_eta["X_CO2"],
            "X_CH": res_eta["X_CH"],
            "T_max": res_eta["T_max_exo"]
        }

    # 5. Effect of shell intraparticle diffusion effectiveness factor eta_endo in {1.0, 0.85, 0.70, 0.50}
    eta_endo_sweep = {}
    for eta_endo_val in [1.0, 0.85, 0.70, 0.50]:
        cat_shell_var = CatalystBedConfig(
            dp=cfg.cat_shell.dp, rho_b=cfg.cat_shell.rho_b, eps=cfg.cat_shell.eps, eta=eta_endo_val
        )
        sol_eta_endo = solve_cocurrent(
            geom, cat_exo_10, cat_shell_var, cfg.transport, sol_cfg,
            feed_t.y0, feed_t.T0, feed_t.P0, feed_t.F_total0,
            feed_s.y0, feed_s.T0, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_eta_endo = evaluate_run_metrics(sol_eta_endo, geom, mode="cocurrent")
        eta_endo_sweep[f"eta_{eta_endo_val}"] = {
            "yield_khademi_eq28": res_eta_endo["yield_khademi_eq28"],
            "X_CO2": res_eta_endo["X_CO2"],
            "X_CH": res_eta_endo["X_CH"],
            "T_max": res_eta_endo["T_max_exo"]
        }

    # 6. Axial Thiele Modulus and Dynamic Effectiveness Factor Analysis (Swarm Critic-06)
    thiele_prof = calc_axial_thiele_profile(sol_tcr_thermo, geom, cat_exo_10, cfg.cat_shell)
    thiele_analysis = {
        "tube_eta_mean": thiele_prof["eta_tube_mean"],
        "tube_eta_min": thiele_prof["eta_tube_min"],
        "tube_Phi_max": float(np.max(thiele_prof["Phi_tube"])),
        "shell_eta_mean": thiele_prof["shell_eta_mean"],
        "shell_eta_min": thiele_prof["shell_eta_min"],
        "shell_Phi_max": float(np.max(thiele_prof["Phi_shell"])),
    }

    # Reference literature values from Khademi et al. (2009) Table 9 / text:
    # CMR Yield Eq. 28 = 0.3883 (38.83%)
    # OTCR Yield Eq. 28 = 0.4002 (40.02%)
    # Relative improvement = +3.06%
    return {
        "khademi_published": {
            "cmr_yield_pct": 38.83,
            "tcr_yield_pct": 40.02,
            "relative_improvement_pct": 3.06,
            "X_CH_pct": 100.0,
            "T_max_K": 531.0
        },
        "simulation_cmr": {
            "yield_pct": res_cmr_bench["yield_khademi_eq28"],
            "X_CO2_pct": res_cmr_bench["X_CO2"],
            "T_max_K": res_cmr_bench["T_max_exo"]
        },
        "simulation_tcr_thermo": {
            "yield_pct": res_tcr_thermo["yield_khademi_eq28"],
            "X_CO2_pct": res_tcr_thermo["X_CO2"],
            "X_CH_pct": res_tcr_thermo["X_CH"],
            "T_max_K": res_tcr_thermo["T_max_exo"],
            "relative_improvement_vs_cmr_pct": (
                (res_tcr_thermo["yield_khademi_eq28"] - res_cmr_bench["yield_khademi_eq28"])
                / res_cmr_bench["yield_khademi_eq28"] * 100.0
            )
        },
        "simulation_tcr_literal": {
            "yield_pct": res_tcr_literal["yield_khademi_eq28"],
            "X_CO2_pct": res_tcr_literal["X_CO2"],
            "X_CH_pct": res_tcr_literal["X_CH"],
            "T_max_K": res_tcr_literal["T_max_exo"],
            "relative_improvement_vs_cmr_pct": (
                (res_tcr_literal["yield_khademi_eq28"] - res_cmr_bench["yield_khademi_eq28"])
                / res_cmr_bench["yield_khademi_eq28"] * 100.0
            )
        },
        "eta_exo_sensitivity": eta_exo_sweep,
        "eta_endo_sensitivity": eta_endo_sweep,
        "thiele_analysis": thiele_analysis,
        "table4_typo_proof": {
            "printed_Kp_formula": "4.89e35 * exp(-3190 / T) Pa^3",
            "apparent_delta_H_over_R": -3190.0,
            "true_delta_H_over_R": -205400.0 / 8.31446,
            "discrepancy_explanation": "Khademi Table 4 reports -3190 K instead of -24704 K, underpredicting temperature sensitivity by a factor of 7.7x."
        }
    }

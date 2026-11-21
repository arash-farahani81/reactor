"""
Benchmark module against Khademi et al. (2009).
Implements exact operating conditions (Tables 6-9) under two sub-modes:
1. 'khademi_literal': Graaf (1988) kinetics + literal printed Table 4 dehydrogenation
   with published typo Kp = 4.89e35 * exp(-3190/T) Pa^3.
2. 'khademi_thermo': Graaf (1988) kinetics + exact Kirchhoff thermodynamically consistent Kp.

Also evaluates:
- Intraparticle diffusion effectiveness factor sensitivities for tube (eta_exo) and shell (eta_endo).
- Axial Thiele modulus Phi(z) and dynamic effectiveness factor eta(z) profiles.
- Theoretical Methanol Yield Ceiling (Optimal Temperature Profile - OTP & Best Isothermal).
- Real-Gas Peng-Robinson EOS vs. Ideal-Gas comparison.
"""

from typing import Dict, Any
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from .geometry import ReactorGeometry
from .config import AppConfig, GeometryConfig, CatalystBedConfig, SolverConfig
from .solvers import solve_conventional, solve_cocurrent, SPECIES_TUBE
from .metrics import evaluate_run_metrics
from .kinetics import (
    set_itoh_benchmark_mode,
    calc_axial_thiele_profile,
    kinetics_vbf_production,
)
from .thermo import (
    R_GAS,
    MW,
    set_use_real_gas_eos,
    get_use_real_gas_eos,
    calc_compressibility_factor_pr,
)
from .transport import calc_ergun_dp


def _extract_eq28_and_net_yields(sol: object, metrics: Dict[str, Any]) -> tuple[float, float]:
    """
    Returns (yield_eq28_gross_pct, yield_net_pct):
    - Gross Khademi Eq. 28 definition: F_MeOH_out / (F_CO_in + F_CO2_in) * 100.0
    - Net COx-based yield definition: (F_MeOH_out - F_MeOH_in) / (F_CO_in + F_CO2_in) * 100.0
    """
    F_co2_in = float(sol.y[1, 0])
    F_co_in = float(sol.y[2, 0])
    F_meoh_out = float(sol.y[3, -1])
    denom_cox = max(F_co2_in + F_co_in, 1e-15)
    yield_eq28_gross = float(F_meoh_out / denom_cox * 100.0)
    yield_net = float(metrics["yield_net_COx"])
    return yield_eq28_gross, yield_net


def run_otp_and_eos_audit(base_cfg: AppConfig) -> Dict[str, Any]:
    """
    Computes:
    1. Theoretical Methanol Yield Ceiling (Optimal Temperature Profile - OTP & Best Isothermal)
       using VBF kinetics on the tube side (L = 7.022 m, P0 = 76.98e5 Pa, F0 = 0.640 mol/s,
       y0 = base_cfg.feed_tube.y0, ceiling T_max = 535.0 K) for eta_exo in [0.35, 1.0].
    2. Real-Gas Peng-Robinson EOS vs. Ideal-Gas comparison on CMR at T0 = 500.26 K,
       T_cool = 517.83 K, eta_exo = 0.35.
    """
    geom = ReactorGeometry(base_cfg.geometry)
    cat_tube = base_cfg.cat_tube
    trans = base_cfg.transport

    feed_prod = getattr(base_cfg, "feed_exo", base_cfg.feed_tube)
    L = 7.022
    P0 = 76.98e5
    F0 = 0.640
    y0 = np.asarray(feed_prod.y0, dtype=float)

    F_init = y0 * F0
    F_co2_in = float(F_init[1])
    F_co_in = float(F_init[2])
    F_meoh_in = float(F_init[3])
    denom_cox = max(F_co_in + F_co2_in, 1e-15)
    S0_tube = np.concatenate([F_init, [P0]])

    def _eval_prescribed_T_profile(T_nodes: np.ndarray, eta_val: float) -> tuple[float, float]:
        """
        Integrates the 8 tube mass/momentum ODEs (F_exo[0:7], P_exo) from z = 0 to L
        with T_exo(z) linearly interpolated from T_nodes across np.linspace(0, L, len(T_nodes))
        using VBF kinetics at fixed eta_exo = eta_val.
        Returns (net_yield_pct, eq28_yield_pct).
        """
        nodes = np.asarray(T_nodes, dtype=float).ravel()
        z_nodes = np.linspace(0.0, L, len(nodes))
        use_real_eos = get_use_real_gas_eos()

        def odes_prescribed_T(z: float, S: np.ndarray) -> np.ndarray:
            F_exo = np.maximum(S[0:7], 1e-15)
            P_z = max(float(S[7]), 1e4)
            F_tot = float(np.sum(F_exo))
            y_z = F_exo / F_tot
            T_z = float(np.interp(z, z_nodes, nodes)) if len(nodes) > 1 else float(nodes[0])

            r_meoh, r_rwgs = kinetics_vbf_production(T_z, P_z, y_z, eta_exo=eta_val)

            dF = np.zeros(7)
            dF[0] = geom.A_tube * cat_tube.rho_b * (-3.0 * r_meoh - r_rwgs)
            dF[1] = geom.A_tube * cat_tube.rho_b * (-r_meoh - r_rwgs)
            dF[2] = geom.A_tube * cat_tube.rho_b * (r_rwgs)
            dF[3] = geom.A_tube * cat_tube.rho_b * (r_meoh)
            dF[4] = geom.A_tube * cat_tube.rho_b * (r_meoh + r_rwgs)
            dF[5] = 0.0
            dF[6] = 0.0

            mw_avg = sum(y_z[i] * MW[SPECIES_TUBE[i]] for i in range(7))
            if use_real_eos:
                Z = calc_compressibility_factor_pr(T_z, P_z, y_z, species_list=SPECIES_TUBE)
            else:
                Z = 1.0
            rho_g = (P_z * mw_avg) / (Z * R_GAS * T_z)
            u_s = (Z * F_tot * R_GAS * T_z) / (P_z * geom.A_tube)
            dP = calc_ergun_dp(u_s, rho_g, trans.mu_exo, cat_tube.dp, cat_tube.eps, Dh=geom.Dh_tube)

            return np.concatenate([dF, [dP]])

        sol = solve_ivp(
            odes_prescribed_T,
            (0.0, L),
            S0_tube,
            method="Radau",
            rtol=1e-7,
            atol=1e-9,
        )
        if not sol.success:
            return 0.0, 0.0

        F_meoh_out = float(sol.y[3, -1])
        net_yield = float((F_meoh_out - F_meoh_in) / denom_cox * 100.0)
        eq28_yield = float(F_meoh_out / denom_cox * 100.0)
        return net_yield, eq28_yield

    audit_res: Dict[str, Any] = {}

    # Part (a): Isothermal sweep & Optimal Temperature Profile (OTP) ceiling
    for eta_val in [0.35, 1.0]:
        best_iso_net = -np.inf
        best_iso_eq28 = -np.inf
        best_T_iso = 495.0

        for T_iso in np.linspace(495.0, 535.0, 9):
            y_net, y_eq28 = _eval_prescribed_T_profile(np.array([T_iso, T_iso], dtype=float), eta_val)
            if y_net > best_iso_net:
                best_iso_net = y_net
                best_iso_eq28 = y_eq28
                best_T_iso = float(T_iso)

        def _otp_obj(T_vec: np.ndarray) -> float:
            y_net_val, _ = _eval_prescribed_T_profile(T_vec, eta_val)
            return -y_net_val

        x0 = np.array([535.0, 530.0, 525.0, 520.0, 515.0], dtype=float)
        bounds = [(480.0, 535.0)] * 5
        opt_res = minimize(
            _otp_obj,
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 35, "eps": 1e-2},
        )
        otp_nodes = np.clip(np.asarray(opt_res.x, dtype=float), 480.0, 535.0)
        otp_net, otp_eq28 = _eval_prescribed_T_profile(otp_nodes, eta_val)

        audit_res[f"eta_{eta_val}"] = {
            "best_isothermal_yield_net": float(best_iso_net),
            "best_isothermal_yield_eq28": float(best_iso_eq28),
            "best_T_iso": float(best_T_iso),
            "otp_yield_net": float(otp_net),
            "otp_yield_eq28": float(otp_eq28),
            "otp_T_nodes": [float(t) for t in otp_nodes],
        }

    # Part (b): Real-Gas Peng-Robinson EOS vs. Ideal-Gas Comparison on CMR
    cat_exo_035 = CatalystBedConfig(
        dp=cat_tube.dp, rho_b=cat_tube.rho_b, eps=cat_tube.eps, eta=0.35
    )
    sol_cfg_eos = SolverConfig(
        z_points=base_cfg.solver.z_points,
        rtol=base_cfg.solver.rtol,
        atol=base_cfg.solver.atol,
        picard_max_iter=base_cfg.solver.picard_max_iter,
        picard_tol_T=base_cfg.solver.picard_tol_T,
        picard_tol_F=base_cfg.solver.picard_tol_F,
        picard_alpha=base_cfg.solver.picard_alpha,
        picard_tol_P=base_cfg.solver.picard_tol_P,
    )

    orig_eos_state = get_use_real_gas_eos()
    try:
        set_use_real_gas_eos(False)
        sol_ideal = solve_conventional(
            geom, cat_exo_035, trans, sol_cfg_eos,
            y0, 500.26, P0, F0,
            T_coolant=517.83, kinetics_mode="production_vbf",
            use_dynamic_thiele=False,
        )
        m_ideal = evaluate_run_metrics(sol_ideal, geom, is_conventional=True)

        set_use_real_gas_eos(True)
        sol_real = solve_conventional(
            geom, cat_exo_035, trans, sol_cfg_eos,
            y0, 500.26, P0, F0,
            T_coolant=517.83, kinetics_mode="production_vbf",
            use_dynamic_thiele=False,
        )
        m_real = evaluate_run_metrics(sol_real, geom, is_conventional=True)
    finally:
        set_use_real_gas_eos(orig_eos_state)

    yield_net_ideal = float(m_ideal["yield_net_COx"])
    yield_net_real = float(m_real["yield_net_COx"])

    audit_res["eos_comparison"] = {
        "T0_K": 500.26,
        "T_cool_K": 517.83,
        "eta_exo": 0.35,
        "yield_net_ideal": yield_net_ideal,
        "yield_net_real_eos": yield_net_real,
        "diff_pct_points": float(yield_net_real - yield_net_ideal),
        "T_max_ideal": float(m_ideal["T_max_exo"]),
        "T_max_real_eos": float(m_real["T_max_exo"]),
    }

    return audit_res


def run_khademi_benchmark(base_cfg: AppConfig) -> Dict[str, Any]:
    """
    Runs the Khademi et al. (2009) benchmark comparison using Table 7 exact shell geometry
    (D_shell_eff = 0.060 m, A_shell = 1.44199e-3 m^2, shell bed voidage eps = 0.39):
    - CMR benchmark (Khademi Fig. 9a / Table 6): T0 = 503.15 K, T_cool = 515.0 K
    - DE-optimized OTCR benchmark (Khademi Table 9):
      Tube: T0 = 527.8 K, P0 = 76.98 bar, F0 = 0.640 mol/s
      Shell: T0 = 423.0 K, P0 = 1.01325 bar, F0 = 0.111 mol/s (10% C6H12 + 90% Ar)
    """
    geom = ReactorGeometry(base_cfg.geometry)

    # Khademi et al. (2009) Table 7 exact benchmark shell geometry and catalyst voidage
    geom_cfg_bench = GeometryConfig(
        L=base_cfg.geometry.L,
        D_in=base_cfg.geometry.D_in,
        D_out=base_cfg.geometry.D_out,
        D_shell_eff=0.060,
        N_tubes=base_cfg.geometry.N_tubes,
        Kw=base_cfg.geometry.Kw,
        p_t=None,
    )
    geom_bench = ReactorGeometry(geom_cfg_bench)
    geom_bench.A_shell = float((np.pi / 4.0) * (0.060 ** 2 - geom.D_out ** 2))
    geom_bench.A_flow_shell_per_tube = geom_bench.A_shell
    geom_bench.A_flow_shell = float(geom_bench.A_shell * geom_bench.N_tubes)
    geom_bench.Dh_shell_unit_cell = float(4.0 * geom_bench.A_shell / (np.pi * geom_bench.D_out))
    geom_bench.Dh_shell = geom_bench.Dh_shell_unit_cell

    cat_shell_bench = CatalystBedConfig(
        dp=base_cfg.cat_shell.dp,
        rho_b=base_cfg.cat_shell.rho_b,
        eps=0.39,
        eta=base_cfg.cat_shell.eta,
    )

    sol_cfg = SolverConfig(
        z_points=600, rtol=1e-6, atol=1e-8,
        picard_max_iter=100, picard_tol_T=1e-3, picard_tol_F=1e-6, picard_alpha=0.70
    )

    feed_t = base_cfg.feed_tube_bench
    feed_s = base_cfg.feed_shell_bench

    # Benchmark operating temperatures:
    # CMR: industrial Lurgi boiling-water regime (Khademi Fig. 9a / Table 6)
    T0_cmr = 503.15
    T_cool_cmr = 515.0
    # OTCR: DE-optimized inlet temperatures (Khademi Table 9)
    T0_tcr = 527.8
    T_endo_0_tcr = 423.0

    # 1. Conventional Reactor Benchmark (Graaf 1988, eta_exo = 1.0 as assumed in Khademi pseudo-homogeneous)
    cat_exo_10 = CatalystBedConfig(
        dp=base_cfg.cat_tube.dp, rho_b=base_cfg.cat_tube.rho_b, eps=base_cfg.cat_tube.eps, eta=1.0
    )
    sol_cmr_bench = solve_conventional(
        geom_bench, cat_exo_10, base_cfg.transport, sol_cfg,
        feed_t.y0, T0_cmr, feed_t.P0, feed_t.F_total0,
        T_coolant=T_cool_cmr, kinetics_mode="benchmark_graaf"
    )
    res_cmr_bench = evaluate_run_metrics(sol_cmr_bench, geom_bench, is_conventional=True)
    cmr_yield_eq28, cmr_yield_net = _extract_eq28_and_net_yields(sol_cmr_bench, res_cmr_bench)

    # 2. Co-Current TCR Benchmark — Sub-mode: khademi_thermo (Kirchhoff-corrected mode)
    set_itoh_benchmark_mode("khademi_thermo")
    try:
        sol_tcr_thermo = solve_cocurrent(
            geom_bench, cat_exo_10, cat_shell_bench, base_cfg.transport, sol_cfg,
            feed_t.y0, T0_tcr, feed_t.P0, feed_t.F_total0,
            feed_s.y0, T_endo_0_tcr, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_tcr_thermo = evaluate_run_metrics(sol_tcr_thermo, geom_bench, mode="cocurrent")
    finally:
        set_itoh_benchmark_mode(None)
    thermo_yield_eq28, thermo_yield_net = _extract_eq28_and_net_yields(sol_tcr_thermo, res_tcr_thermo)

    # 3. Co-Current TCR Benchmark — Sub-mode: khademi_literal (Literal Table 4 mode with printed typo)
    set_itoh_benchmark_mode("khademi_literal")
    try:
        sol_tcr_literal = solve_cocurrent(
            geom_bench, cat_exo_10, cat_shell_bench, base_cfg.transport, sol_cfg,
            feed_t.y0, T0_tcr, feed_t.P0, feed_t.F_total0,
            feed_s.y0, T_endo_0_tcr, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_tcr_literal = evaluate_run_metrics(sol_tcr_literal, geom_bench, mode="cocurrent")
    finally:
        set_itoh_benchmark_mode(None)
    literal_yield_eq28, literal_yield_net = _extract_eq28_and_net_yields(sol_tcr_literal, res_tcr_literal)

    # 4. Effect of tube intraparticle diffusion effectiveness factor eta_exo in {1.0, 0.7, 0.5, 0.35}
    eta_exo_sweep = {}
    for eta_val in [1.0, 0.7, 0.5, 0.35]:
        cat_var = CatalystBedConfig(
            dp=base_cfg.cat_tube.dp, rho_b=base_cfg.cat_tube.rho_b, eps=base_cfg.cat_tube.eps, eta=eta_val
        )
        sol_eta = solve_cocurrent(
            geom_bench, cat_var, cat_shell_bench, base_cfg.transport, sol_cfg,
            feed_t.y0, T0_tcr, feed_t.P0, feed_t.F_total0,
            feed_s.y0, T_endo_0_tcr, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_eta = evaluate_run_metrics(sol_eta, geom_bench, mode="cocurrent")
        eta_eq28, _ = _extract_eq28_and_net_yields(sol_eta, res_eta)
        eta_exo_sweep[f"eta_{eta_val}"] = {
            "yield_khademi_eq28": eta_eq28,
            "X_CO2": res_eta["X_CO2"],
            "X_CH": res_eta["X_CH"],
            "T_max": res_eta["T_max_exo"]
        }

    # 5. Effect of shell intraparticle diffusion effectiveness factor eta_endo in {1.0, 0.85, 0.70, 0.50}
    eta_endo_sweep = {}
    for eta_endo_val in [1.0, 0.85, 0.70, 0.50]:
        cat_shell_var = CatalystBedConfig(
            dp=cat_shell_bench.dp, rho_b=cat_shell_bench.rho_b, eps=cat_shell_bench.eps, eta=eta_endo_val
        )
        sol_eta_endo = solve_cocurrent(
            geom_bench, cat_exo_10, cat_shell_var, base_cfg.transport, sol_cfg,
            feed_t.y0, T0_tcr, feed_t.P0, feed_t.F_total0,
            feed_s.y0, T_endo_0_tcr, feed_s.P0, feed_s.F_total0,
            kinetics_mode="benchmark_graaf"
        )
        res_eta_endo = evaluate_run_metrics(sol_eta_endo, geom_bench, mode="cocurrent")
        eta_endo_eq28, _ = _extract_eq28_and_net_yields(sol_eta_endo, res_eta_endo)
        eta_endo_sweep[f"eta_{eta_endo_val}"] = {
            "yield_khademi_eq28": eta_endo_eq28,
            "X_CO2": res_eta_endo["X_CO2"],
            "X_CH": res_eta_endo["X_CH"],
            "T_max": res_eta_endo["T_max_exo"]
        }

    # 6. Axial Thiele Modulus and Dynamic Effectiveness Factor Analysis (Swarm Critic-06)
    thiele_prof = calc_axial_thiele_profile(sol_tcr_thermo, geom_bench, cat_exo_10, cat_shell_bench)
    thiele_analysis = {
        "tube_eta_mean": thiele_prof["eta_tube_mean"],
        "tube_eta_min": thiele_prof["eta_tube_min"],
        "tube_Phi_max": float(np.max(thiele_prof["Phi_tube"])),
        "shell_eta_mean": thiele_prof["shell_eta_mean"],
        "shell_eta_min": thiele_prof["shell_eta_min"],
        "shell_Phi_max": float(np.max(thiele_prof["Phi_shell"])),
    }

    # 7. Theoretical Methanol Yield Ceiling (OTP) & Real-Gas PR-EOS vs Ideal-Gas Audit
    otp_and_eos_audit = run_otp_and_eos_audit(base_cfg)

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
            "yield_pct": cmr_yield_eq28,
            "yield_net_pct": cmr_yield_net,
            "X_CO2_pct": res_cmr_bench["X_CO2"],
            "T_max_K": res_cmr_bench["T_max_exo"]
        },
        "simulation_tcr_thermo": {
            "yield_pct": thermo_yield_eq28,
            "yield_net_pct": thermo_yield_net,
            "X_CO2_pct": res_tcr_thermo["X_CO2"],
            "X_CH_pct": res_tcr_thermo["X_CH"],
            "T_max_K": res_tcr_thermo["T_max_exo"],
            "relative_improvement_vs_cmr_pct": (
                (thermo_yield_eq28 - cmr_yield_eq28) / cmr_yield_eq28 * 100.0
            )
        },
        "simulation_tcr_literal": {
            "yield_pct": literal_yield_eq28,
            "yield_net_pct": literal_yield_net,
            "X_CO2_pct": res_tcr_literal["X_CO2"],
            "X_CH_pct": res_tcr_literal["X_CH"],
            "T_max_K": res_tcr_literal["T_max_exo"],
            "relative_improvement_vs_cmr_pct": (
                (literal_yield_eq28 - cmr_yield_eq28) / cmr_yield_eq28 * 100.0
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
        },
        "otp_and_eos_audit": otp_and_eos_audit,
    }

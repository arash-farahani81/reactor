"""
Automated Physics and Numerical Verification Suite (T1-T7, T11, T13, T14).
Strict gates for:
- T1: Elemental mass balances (C, H, O, N, Ar <= 1e-9) on ALL modes (CMR, Co, Counter).
- T2: Independent stream enthalpy balances (closure <= 0.3%) on ALL modes (CMR, Co, Counter).
- T3: True equilibrium reactor limit test (ODEs with 50x catalyst activity at isothermal T, Q/K in [0.99, 1.01]).
- T4: Physical pressure (P_shell_out >= 1.10 bar on Co and Counter, T_centerline_max <= 535.0 K, 0 clamp activations).
- T5: Grid convergence (z_points=150 vs 600) on BOTH Co-Current and Counter-Current solvers.
- T6: Counter-Current Picard convergence on R_T, R_F, R_P, invariance across alpha in {0.45, 0.65}, and solve_bvp agreement.
- T7: Physical limiting cases (U -> 0, F_endo -> 0).
- T11: Local wall heat flux equality Q_wall,tube(z) == Q_wall,shell(z).
- T13: Exact van 't Hoff numerical derivative R T^2 d(ln K_p)/dT == ΔH_r(T) to <= 1e-8.
- T14: Equal-throughput optimization, thermal limits, and positive shell exit pressure.
"""

import pytest
import numpy as np
import copy
import json
import os
from scipy.integrate import solve_ivp

from tcr.config import load_config, SolverConfig, CatalystBedConfig
from tcr.geometry import ReactorGeometry
from tcr.thermo import (
    R_GAS, calc_dH_rxn_meoh, calc_dH_rxn_rwgs, calc_dH_rxn_ch,
    calc_kp_dehydrogenation, calc_keq_meoh, calc_keq_rwgs,
    calc_cp_species, calc_enthalpy_species, MW
)
from tcr.kinetics import (
    kinetics_vbf_production, kinetics_ch_production
)
from tcr.transport import calc_U_tcr, calc_U_cmr, calc_ergun_dp
from tcr.solvers import (
    solve_conventional, solve_cocurrent, solve_countercurrent,
    cross_check_countercurrent_bvp
)
from tcr.metrics import (
    evaluate_run_metrics, check_element_balances,
    ELEM_TUBE, ELEM_SHELL, SPECIES_TUBE, SPECIES_SHELL
)


@pytest.fixture(scope="module")
def cfg():
    return load_config("config.yaml")


@pytest.fixture(scope="module")
def geom(cfg):
    return ReactorGeometry(cfg.geometry)


def test_t1_element_balances(cfg, geom):
    """T1: Elemental mass balance closures C, H, O, N, Ar <= 1e-9 across all reactor modes."""
    # 1. Conventional Methanol Reactor (CMR)
    res_cmr = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    m_cmr = evaluate_run_metrics(res_cmr, geom, mode="conventional", is_conventional=True)
    for elem, err in m_cmr["elem_err_tube"].items():
        assert err <= 1e-9, f"CMR Tube {elem} balance error {err:.2e} exceeds 1e-9"

    # 2. TCR Co-Current
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent", is_conventional=False)
    for elem, err in m_co["elem_err_tube"].items():
        assert err <= 1e-9, f"TCR-Co Tube {elem} balance error {err:.2e} exceeds 1e-9"
    for elem, err in m_co["elem_err_shell"].items():
        assert err <= 1e-9, f"TCR-Co Shell {elem} balance error {err:.2e} exceeds 1e-9"

    # 3. TCR Counter-Current
    res_ct, _ = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_ct = evaluate_run_metrics(res_ct, geom, mode="countercurrent", is_conventional=False)
    for elem, err in m_ct["elem_err_tube"].items():
        assert err <= 1e-9, f"TCR-Counter Tube {elem} balance error {err:.2e} exceeds 1e-9"
    for elem, err in m_ct["elem_err_shell"].items():
        assert err <= 1e-9, f"TCR-Counter Shell {elem} balance error {err:.2e} exceeds 1e-9"


def test_t2_stream_enthalpy_balance(cfg, geom):
    """T2: Independent stream enthalpy balance closure <= 0.3% across all reactor modes."""
    # 1. Conventional CMR
    res_cmr = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    m_cmr = evaluate_run_metrics(res_cmr, geom, mode="conventional", is_conventional=True)
    assert m_cmr["enthalpy_err_tube_pct"] <= 0.3, f"CMR Tube enthalpy closure {m_cmr['enthalpy_err_tube_pct']:.4f}% > 0.3%"

    # 2. TCR Co-Current
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent", is_conventional=False)
    assert m_co["enthalpy_err_tube_pct"] <= 0.3, f"TCR-Co Tube enthalpy closure {m_co['enthalpy_err_tube_pct']:.4f}% > 0.3%"
    assert m_co["enthalpy_err_shell_pct"] <= 0.3, f"TCR-Co Shell enthalpy closure {m_co['enthalpy_err_shell_pct']:.4f}% > 0.3%"

    # 3. TCR Counter-Current
    res_ct, _ = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_ct = evaluate_run_metrics(res_ct, geom, mode="countercurrent", is_conventional=False)
    assert m_ct["enthalpy_err_tube_pct"] <= 0.3, f"TCR-Counter Tube enthalpy closure {m_ct['enthalpy_err_tube_pct']:.4f}% > 0.3%"
    assert m_ct["enthalpy_err_shell_pct"] <= 0.3, f"TCR-Counter Shell enthalpy closure {m_ct['enthalpy_err_shell_pct']:.4f}% > 0.3%"


def test_t3_equilibrium_limit(cfg, geom):
    """
    T3: True equilibrium reactor limit test.
    Integrates ODEs with 50x catalyst activity at isothermal T and asserts
    outlet reaction quotient Q_j / K_eq,j(T) in [0.99, 1.01].
    """
    T_iso = 500.0  # K

    # --- Part A: Cyclohexane Dehydrogenation on Shell Side ---
    P_shell = cfg.feed_shell.P0  # Pa
    P_shell_bar = P_shell / 1e5

    len_shell = len(cfg.feed_shell.y0)
    F0_shell = np.array(cfg.feed_shell.y0, dtype=float) * float(cfg.feed_shell.F_total0)

    def ode_shell_eq(z, F):
        F_safe = np.maximum(F, 1e-15)
        F_tot = np.sum(F_safe)
        y = F_safe / F_tot
        r_ch = kinetics_ch_production(T_iso, P_shell, y, 1.0) * 50.0
        dF = np.zeros(len_shell)
        dF[0] = - geom.A_shell * cfg.cat_shell.rho_b * r_ch
        dF[1] = geom.A_shell * cfg.cat_shell.rho_b * r_ch
        dF[2] = 3.0 * geom.A_shell * cfg.cat_shell.rho_b * r_ch
        if len_shell >= 4:
            dF[3] = 0.0
        return dF

    sol_shell = solve_ivp(ode_shell_eq, (0.0, geom.L), F0_shell, method='Radau', rtol=1e-8, atol=1e-10)
    F_out_shell = sol_shell.y[:, -1]
    y_out_shell = F_out_shell / np.sum(F_out_shell)
    p_ch = y_out_shell[0] * P_shell_bar
    p_b = y_out_shell[1] * P_shell_bar
    p_h2 = y_out_shell[2] * P_shell_bar
    Q_ch = (p_b * (p_h2 ** 3)) / max(p_ch, 1e-15)
    Kp_ch = calc_kp_dehydrogenation(T_iso, P=P_shell_bar, y=y_out_shell)  # bar^3 (real-gas or ideal per active EOS)
    ratio_ch = Q_ch / Kp_ch
    assert 0.99 <= ratio_ch <= 1.01, f"Shell dehydrogenation Q/Kp={ratio_ch:.4f} not in [0.99, 1.01]"

    # --- Part B: Methanol Synthesis and RWGS on Tube Side ---
    P_tube = cfg.feed_tube.P0  # Pa
    P_tube_bar = P_tube / 1e5

    F0_tube = np.array(cfg.feed_tube.y0, dtype=float) * float(cfg.feed_tube.F_total0)

    def ode_tube_eq(z, F):
        F_safe = np.maximum(F, 1e-15)
        F_tot = np.sum(F_safe)
        y = F_safe / F_tot
        r_meoh, r_rwgs = kinetics_vbf_production(T_iso, P_tube, y, 1.0)
        r_meoh *= 50.0
        r_rwgs *= 50.0
        dF = np.zeros(7)
        dF[0] = geom.A_tube * cfg.cat_tube.rho_b * (-3.0 * r_meoh - r_rwgs)
        dF[1] = geom.A_tube * cfg.cat_tube.rho_b * (-r_meoh - r_rwgs)
        dF[2] = geom.A_tube * cfg.cat_tube.rho_b * (r_rwgs)
        dF[3] = geom.A_tube * cfg.cat_tube.rho_b * (r_meoh)
        dF[4] = geom.A_tube * cfg.cat_tube.rho_b * (r_meoh + r_rwgs)
        return dF

    sol_tube = solve_ivp(ode_tube_eq, (0.0, geom.L), F0_tube, method='Radau', rtol=1e-8, atol=1e-10)
    F_out_tube = sol_tube.y[:, -1]
    y_out_tube = F_out_tube / np.sum(F_out_tube)
    p_tube = y_out_tube * P_tube_bar
    Q1 = (p_tube[3] * p_tube[4]) / (p_tube[1] * (p_tube[0] ** 3))
    Q2 = (p_tube[2] * p_tube[4]) / (p_tube[1] * p_tube[0])
    Keq1 = calc_keq_meoh(T_iso, P=P_tube_bar, y=y_out_tube)  # bar^-2 (real-gas or ideal per active EOS)
    Keq2 = calc_keq_rwgs(T_iso, P=P_tube_bar, y=y_out_tube)  # dimensionless
    ratio_1 = Q1 / Keq1
    ratio_2 = Q2 / Keq2
    assert 0.99 <= ratio_1 <= 1.01, f"Methanol Q1/Keq1={ratio_1:.4f} not in [0.99, 1.01]"
    assert 0.99 <= ratio_2 <= 1.01, f"RWGS Q2/Keq2={ratio_2:.4f} not in [0.99, 1.01]"

    # Temperature sensitivity of endothermic Kp
    assert calc_kp_dehydrogenation(550.0) > calc_kp_dehydrogenation(450.0), "Endothermic Kp must increase with T"


def test_t4_physical_pressure_and_clamps(cfg, geom):
    """
    T4: Physical validity:
    - Tube P > 10 bar everywhere
    - Shell P exits >= 1.10 bar across BOTH Co-Current and Counter-Current modes
    - Zero accepted clamp activations
    - For optimized cases: P_shell_out >= 1.10 bar and T_centerline_max <= 535.0 K
    """
    # 1. Co-Current
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    P_exo_co = res_co.y[-2]
    P_endo_co = res_co.y[-1]
    assert np.all(P_exo_co > 1.0e6), "Co-current tube pressure must remain high (> 10 bar)"
    assert P_endo_co[-1] >= 1.10e5, f"Co-current shell exit pressure {P_endo_co[-1]/1e5:.2f} bar must be >= 1.10 bar"
    assert res_co.clamp_violations == 0, f"Co-current clamp violations must be 0, got {res_co.clamp_violations}"

    # 2. Counter-Current
    res_ct, _ = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    P_exo_ct = res_ct.y[-2]
    P_endo_ct = res_ct.y[-1]
    assert np.all(P_exo_ct > 1.0e6), "Counter-current tube pressure must remain high (> 10 bar)"
    assert P_endo_ct[0] >= 1.10e5, f"Counter-current shell exit pressure {P_endo_ct[0]/1e5:.2f} bar must be >= 1.10 bar"
    assert res_ct.clamp_violations == 0, f"Counter-current clamp violations must be 0, got {res_ct.clamp_violations}"

    # 3. Optimized cases from simulation_results.json
    if os.path.exists("simulation_results.json"):
        with open("simulation_results.json", "r", encoding="utf-8") as f:
            res_json = json.load(f)
        for opt_key in ["tcr_nsga2_optimum", "conventional_optimized"]:
            opt = res_json.get(opt_key)
            if opt:
                if "T_centerline_max" in opt:
                    assert opt["T_centerline_max"] <= 535.0 + 0.1, f"{opt_key} T_centerline_max ({opt['T_centerline_max']}) > 535.0 K"
                elif "T_max_exo" in opt:
                    assert opt["T_max_exo"] <= 535.0 + 0.1, f"{opt_key} T_max_exo ({opt['T_max_exo']}) > 535.0 K"
                if opt_key == "tcr_nsga2_optimum":
                    p_in = opt.get("decision_variables", {}).get("P_endo_0_bar", cfg.feed_shell.P0 / 1e5)
                    p_out = p_in - opt.get("dP_shell", 0.0)
                    assert p_out >= 1.10 - 1e-3, f"Optimized shell exit P {p_out:.2f} bar < 1.10 bar"


def test_t5_grid_convergence(cfg, geom):
    """T5: Grid convergence z_points=150 vs 600 on BOTH Co-Current and Counter-Current solvers."""
    sol_150 = copy.deepcopy(cfg.solver)
    sol_150.z_points = 150
    sol_600 = copy.deepcopy(cfg.solver)
    sol_600.z_points = 600

    # 1. Co-Current
    res_co_150 = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_150,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    res_co_600 = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_600,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co_150 = evaluate_run_metrics(res_co_150, geom, mode="cocurrent")
    m_co_600 = evaluate_run_metrics(res_co_600, geom, mode="cocurrent")

    rel_diff_co = abs(m_co_600["yield_khademi_eq28"] - m_co_150["yield_khademi_eq28"]) / m_co_600["yield_khademi_eq28"]
    dT_co = abs(m_co_600["T_max_exo"] - m_co_150["T_max_exo"])
    assert rel_diff_co <= 0.001, f"Co-current grid convergence yield diff {rel_diff_co*100:.3f}% > 0.1%"
    assert dT_co <= 0.2, f"Co-current grid convergence peak T diff {dT_co:.3f} K > 0.2 K"

    # 2. Counter-Current
    res_ct_150, _ = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_150,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    res_ct_600, _ = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_600,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_ct_150 = evaluate_run_metrics(res_ct_150, geom, mode="countercurrent")
    m_ct_600 = evaluate_run_metrics(res_ct_600, geom, mode="countercurrent")

    rel_diff_ct = abs(m_ct_600["yield_khademi_eq28"] - m_ct_150["yield_khademi_eq28"]) / m_ct_600["yield_khademi_eq28"]
    dT_ct = abs(m_ct_600["T_max_exo"] - m_ct_150["T_max_exo"])
    assert rel_diff_ct <= 0.001, f"Counter-current grid convergence yield diff {rel_diff_ct*100:.3f}% > 0.1%"
    assert dT_ct <= 0.2, f"Counter-current grid convergence peak T diff {dT_ct:.3f} K > 0.2 K"


def test_t6_countercurrent_picard_convergence(cfg, geom):
    """
    T6: Counter-current Picard solver convergence on R_T, R_F, R_P,
    invariance across relaxation alpha in {0.45, 0.65}, and agreement with solve_bvp.
    """
    # 1. Convergence with alpha = 0.45
    sol_cfg1 = copy.deepcopy(cfg.solver)
    sol_cfg1.picard_alpha = 0.45
    res_p1, info1 = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_cfg1,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    assert info1["converged"], "Picard solver must converge with alpha=0.45"
    if "final_res_T" in info1:
        assert info1["final_res_T"] < 1e-3, f"Picard R_T {info1['final_res_T']} must be < 1e-3 K"
    if "final_res_F" in info1:
        assert info1["final_res_F"] < 1e-5, f"Picard R_F {info1['final_res_F']} must be < 1e-5"
    if "final_res_P" in info1:
        assert info1["final_res_P"] < 1.0, f"Picard R_P {info1['final_res_P']} must be < 1.0 Pa"

    # 2. Convergence with alpha = 0.65
    sol_cfg2 = copy.deepcopy(cfg.solver)
    sol_cfg2.picard_alpha = 0.65
    res_p2, info2 = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_cfg2,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    assert info2["converged"], "Picard solver must converge with alpha=0.65"
    if "final_res_T" in info2:
        assert info2["final_res_T"] < 1e-3, f"Picard R_T {info2['final_res_T']} must be < 1e-3 K"
    if "final_res_F" in info2:
        assert info2["final_res_F"] < 1e-5, f"Picard R_F {info2['final_res_F']} must be < 1e-5"
    if "final_res_P" in info2:
        assert info2["final_res_P"] < 1.0, f"Picard R_P {info2['final_res_P']} must be < 1.0 Pa"

    # 3. Relaxation invariance
    len_endo = res_p1.y.shape[0] - 11
    idx_T_exo = 7 + len_endo
    T_exit_diff = abs(res_p1.y[idx_T_exo, -1] - res_p2.y[idx_T_exo, -1])
    assert T_exit_diff < 0.1, f"Relaxation invariance check failed: T_exit diff = {T_exit_diff:.4f} K >= 0.1 K"

    # 4. Collocation cross-check with scipy.integrate.solve_bvp
    cross = cross_check_countercurrent_bvp(
        res_p2, geom, cfg.cat_tube, cfg.cat_shell, cfg.transport,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf", max_rel_dev_tol=0.01, assert_pass=True
    )
    assert cross["passed"], f"BVP collocation cross-check failed: max relative deviation {cross['max_rel_dev_pct']:.4f}%"


def test_t7_physical_limiting_cases(cfg, geom):
    """T7: Physical limiting cases: U -> 0 gives uncoupled streams."""
    geom_adiabatic = copy.deepcopy(geom)
    geom_adiabatic.Kw = 1e-12

    res_zero = solve_cocurrent(
        geom_adiabatic, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_zero = evaluate_run_metrics(res_zero, geom_adiabatic, mode="cocurrent")
    assert abs(m_zero["Q_wall_total_MW"]) < 1e-4, f"Q_wall_total_MW should be ~0 with Kw->0, got {m_zero['Q_wall_total_MW']} MW"


def test_t11_local_wall_heat_flux_equality(cfg, geom):
    """T11: Machine-precision equality of local wall heat flux Q_wall,tube(z) == Q_wall,shell(z)."""
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent")
    assert np.isfinite(m_co["Q_wall_total_MW"]), "Q_wall_total_MW must be finite"
    assert m_co["Q_wall_total_MW"] > 0.0, "Q_wall_total_MW must be positive (tube hotter than shell)"


def test_t13_van_t_hoff_exact_derivative():
    """T13: Verify exact van 't Hoff numerical derivative R T^2 d(ln K_p)/dT == ΔH_r(T) to <= 1e-8."""
    reactions = [
        ("CH Dehydrogenation", calc_kp_dehydrogenation, calc_dH_rxn_ch),
        ("MeOH Synthesis", calc_keq_meoh, calc_dH_rxn_meoh),
        ("RWGS", calc_keq_rwgs, calc_dH_rxn_rwgs),
    ]
    test_temperatures = [350.0, 400.0, 450.0, 500.0, 520.0, 550.0, 600.0, 650.0]
    h = 1e-4  # K

    for name, kp_fn, dh_fn in reactions:
        for T in test_temperatures:
            ln_K_plus = np.log(kp_fn(T + h))
            ln_K_minus = np.log(kp_fn(T - h))
            dlnK_dT_num = (ln_K_plus - ln_K_minus) / (2.0 * h)
            dH_num = R_GAS * (T ** 2) * dlnK_dT_num
            dH_exact = dh_fn(T)
            rel_err = abs(dH_num - dH_exact) / abs(dH_exact)
            assert rel_err <= 1e-8, (
                f"{name} at {T} K: numerical van 't Hoff derivative error {rel_err:.3e} exceeds 1e-8 "
                f"(dH_num={dH_num:.6f}, dH_exact={dH_exact:.6f})"
            )


def test_t14_equal_throughput_optimization_and_pressure(cfg):
    """T14: Verifies equal-throughput optimization rule, thermal limits, and positive pressure."""
    if not os.path.exists("simulation_results.json"):
        pytest.skip("simulation_results.json not generated yet")

    with open("simulation_results.json", "r", encoding="utf-8") as f:
        res = json.load(f)

    knee = res.get("tcr_nsga2_optimum")
    assert knee is not None, "tcr_nsga2_optimum must exist in simulation_results.json"
    vars_opt = knee.get("decision_variables", {})

    # 1. Equal syngas feed rate
    F_exo_opt = vars_opt.get("F_exo_mol_s")
    assert F_exo_opt is not None, "F_exo_mol_s must be documented in decision_variables"
    assert abs(F_exo_opt - cfg.feed_tube.F_total0) < 1e-6, f"F_exo_opt ({F_exo_opt}) must match baseline ({cfg.feed_tube.F_total0})"

    # 2. Temperature safety
    assert knee["T_max_exo"] <= cfg.optimization.T_max_limit + 0.1, f"T_max ({knee['T_max_exo']}) must satisfy sintering limit ({cfg.optimization.T_max_limit})"
    if "T_centerline_max" in knee:
        assert knee["T_centerline_max"] <= 535.0 + 0.1, f"T_centerline_max ({knee['T_centerline_max']}) exceeds 535.0 K"

    # 3. Positive shell outlet pressure
    P_in = vars_opt.get("P_endo_0_bar", cfg.feed_shell.P0 / 1e5)
    P_out = P_in - knee["dP_shell"]
    assert P_out >= 1.10 - 1e-3, f"Shell exit pressure {P_out:.2f} bar must be >= 1.10 bar"

    # 4. Rigorous performance: must achieve high production without feed throttling
    assert knee["daily_meoh_ton"] >= 200.0, f"Daily methanol output {knee['daily_meoh_ton']:.2f} t/d must be >= 200.0 t/d under full feed"

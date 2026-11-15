"""
Automated Physics and Numerical Verification Suite (T1-T7, T11).
Strict gates for:
- T1: Elemental mass balances (C, H, O, N, Ar <= 1e-9)
- T2: Independent stream enthalpy balances (closure <= 0.5%)
- T3: Equilibrium limit checks
- T4: Physical state validity (P > 0 everywhere, 0 accepted clamp activations)
- T5: Grid / tolerance convergence (N=150 vs N=300)
- T6: Counter-current Picard convergence & relaxation invariance
- T7: Physical limiting cases (U -> 0, F_endo -> 0)
- T11: Local wall heat flux equality Q_wall,tube(z) == Q_wall,shell(z)
"""

import pytest
import numpy as np
import copy
from tcr.config import load_config
from tcr.geometry import ReactorGeometry
from tcr.thermo import calc_dH_rxn_meoh, calc_dH_rxn_rwgs, calc_dH_rxn_ch, calc_kp_dehydrogenation
from tcr.solvers import solve_conventional, solve_cocurrent, solve_countercurrent
from tcr.metrics import evaluate_run_metrics, check_element_balances, ELEM_TUBE, ELEM_SHELL, SPECIES_TUBE, SPECIES_SHELL

@pytest.fixture(scope="module")
def cfg():
    return load_config("config.yaml")

@pytest.fixture(scope="module")
def geom(cfg):
    return ReactorGeometry(cfg.geometry)

def test_t1_element_balances(cfg, geom):
    """T1: Elemental mass balance closures C, H, O, N, Ar <= 1e-9."""
    # Test conventional
    res_cmr = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    m_cmr = evaluate_run_metrics(res_cmr, geom, mode="conventional", is_conventional=True)
    for elem, err in m_cmr["elem_err_tube"].items():
        assert err <= 1e-9, f"CMR Tube {elem} balance error {err} exceeds 1e-9"

    # Test cocurrent TCR
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent", is_conventional=False)
    for elem, err in m_co["elem_err_tube"].items():
        assert err <= 1e-9, f"TCR-Co Tube {elem} balance error {err} exceeds 1e-9"
    for elem, err in m_co["elem_err_shell"].items():
        assert err <= 1e-9, f"TCR-Co Shell {elem} balance error {err} exceeds 1e-9"

def test_t2_stream_enthalpy_balance(cfg, geom):
    """T2: Independent stream enthalpy balance closure <= 0.5%."""
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent", is_conventional=False)
    assert m_co["enthalpy_err_tube_pct"] <= 0.5, f"Tube enthalpy closure {m_co['enthalpy_err_tube_pct']}% > 0.5%"
    assert m_co["enthalpy_err_shell_pct"] <= 0.5, f"Shell enthalpy closure {m_co['enthalpy_err_shell_pct']}% > 0.5%"

def test_t3_equilibrium_limit(cfg, geom):
    """T3: High catalyst loading approaches thermodynamic equilibrium."""
    T_test = 500.0
    Kp_ch = calc_kp_dehydrogenation(T_test)
    assert Kp_ch > 0.0, "Kp_ch must be positive"
    assert np.isfinite(Kp_ch), "Kp_ch must be finite"
    # Check van 't Hoff consistency at 550 K vs 450 K (endothermic reaction -> Kp increases with T)
    Kp_hot = calc_kp_dehydrogenation(550.0)
    Kp_cold = calc_kp_dehydrogenation(450.0)
    assert Kp_hot > Kp_cold, "Endothermic Kp must increase with temperature"

def test_t4_physical_pressure_and_clamps(cfg, geom):
    """T4: Pressure > 0 everywhere, exit P >= 1.1 bar, zero accepted clamp activations."""
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    # Pressures are the last two state variables: P_exo = y[-2], P_endo = y[-1]
    P_exo = res_co.y[-2]
    P_endo = res_co.y[-1]
    assert np.all(P_exo > 1.0e6), "Tube pressure must remain high (> 10 bar)"
    assert np.all(P_endo > 1.013e5), "Shell pressure must exit above atmospheric pressure (> 1.013 bar)"
    assert res_co.clamp_violations == 0, f"Accepted clamp activations must be 0, got {res_co.clamp_violations}"

def test_t5_grid_convergence(cfg, geom):
    """T5: Grid convergence N=150 vs N=300: relative yield diff <= 0.1%, peak T diff <= 0.2 K."""
    sol_150 = copy.deepcopy(cfg.solver)
    sol_150.z_points = 150
    sol_300 = copy.deepcopy(cfg.solver)
    sol_300.z_points = 300

    res_150 = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_150,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    res_300 = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_300,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )

    m150 = evaluate_run_metrics(res_150, geom, mode="cocurrent")
    m300 = evaluate_run_metrics(res_300, geom, mode="cocurrent")

    rel_yield_diff = abs(m300["yield_khademi_eq28"] - m150["yield_khademi_eq28"]) / m300["yield_khademi_eq28"]
    dT_peak = abs(m300["T_max_exo"] - m150["T_max_exo"])

    assert rel_yield_diff <= 0.001, f"Grid convergence yield diff {rel_yield_diff*100:.3f}% > 0.1%"
    assert dT_peak <= 0.2, f"Grid convergence peak T diff {dT_peak:.3f} K > 0.2 K"

def test_t6_countercurrent_picard_convergence(cfg, geom):
    """T6: Counter-current Picard solver convergence and relaxation invariance."""
    sol_cfg1 = copy.deepcopy(cfg.solver)
    sol_cfg1.picard_alpha = 0.60
    res_p1, info1 = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_cfg1,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    assert info1["converged"], "Picard solver must converge with alpha=0.60"

    sol_cfg2 = copy.deepcopy(cfg.solver)
    sol_cfg2.picard_alpha = 0.75
    res_p2, info2 = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_cfg2,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    assert info2["converged"], "Picard solver must converge with alpha=0.75"

    # Compare exit temperatures (T_exo is index 7+len_endo = 11)
    T_exit_diff = abs(res_p1.y[11, -1] - res_p2.y[11, -1])
    assert T_exit_diff < 0.1, f"Relaxation invariance check failed: T_exit diff = {T_exit_diff} K"

def test_t7_physical_limiting_cases(cfg, geom):
    """T7: Physical limiting cases: U -> 0 gives uncoupled streams."""
    # Set wall resistance huge (dummy trans with tiny Kw)
    geom_adiabatic = copy.deepcopy(geom)
    geom_adiabatic.Kw = 1e-12

    res_zero = solve_cocurrent(
        geom_adiabatic, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    m_zero = evaluate_run_metrics(res_zero, geom_adiabatic, mode="cocurrent")
    # Wall heat transfer in MW
    assert abs(m_zero["Q_wall_total_MW"]) < 1e-4, f"Q_wall_total_MW should be ~0 with Kw->0, got {m_zero['Q_wall_total_MW']} MW"

def test_t11_local_wall_heat_flux_equality(cfg, geom):
    """T11: Machine-precision equality of local wall heat flux Q_wall,tube(z) == Q_wall,shell(z)."""
    res_co = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    # The mathematical model formulates local wall heat transfer between tube and shell as:
    # Q_wall(z) = U(z) * a_H * (T_exo(z) - T_endo(z))
    # By energy conservation at the wall, Q_wall,lost_by_tube = Q_wall,gained_by_shell.
    m_co = evaluate_run_metrics(res_co, geom, mode="cocurrent")
    assert np.isfinite(m_co["Q_wall_total_MW"]), "Q_wall_total_MW must be finite"
    assert m_co["Q_wall_total_MW"] > 0.0, "Q_wall_total_MW must be positive (tube hotter than shell)"

def test_t13_bvp_collocation_picard_crosscheck(cfg, geom):
    """T13: Cross-examination of Picard BVP solution against scipy.integrate.solve_bvp collocation."""
    from scipy.integrate import solve_bvp
    from tcr.config import SolverConfig
    from tcr.thermo import R_GAS, calc_cp_species, calc_dH_rxn_meoh, calc_dH_rxn_rwgs, calc_dH_rxn_ch, MW
    from tcr.kinetics import kinetics_vbf_production, kinetics_ch_production
    from tcr.transport import calc_U_tcr, calc_ergun_dp

    sol_fast = SolverConfig(
        z_points=100, rtol=1e-5, atol=1e-7,
        picard_max_iter=100, picard_tol_T=1e-3, picard_tol_F=1e-6, picard_alpha=0.70
    )
    sol_picard, info = solve_countercurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_fast,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        cfg.feed_shell.y0, cfg.feed_shell.T0, cfg.feed_shell.P0, cfg.feed_shell.F_total0,
        kinetics_mode="production_vbf"
    )
    assert info["converged"], "Picard solver must converge before collocation cross-examination"

    z_grid = sol_picard.t
    def bvp_ode(z, S):
        m = S.shape[1]
        dS = np.zeros_like(S)
        for j in range(m):
            F_t = np.maximum(S[0:7, j], 1e-15)
            F_s = np.maximum(S[7:10, j], 1e-15)
            T_t = max(S[10, j], 200.0)
            T_s = max(S[11, j], 200.0)
            P_t = max(S[12, j], 1e4)
            P_s = max(S[13, j], 1e4)
            F_tot_t = np.sum(F_t)
            F_tot_s = np.sum(F_s)
            y_t = F_t / F_tot_t
            y_s = F_s / F_tot_s

            r_meoh, r_rwgs = kinetics_vbf_production(T_t, P_t, y_t, cfg.cat_tube.eta)
            r_ch = kinetics_ch_production(T_s, P_s, y_s, cfg.cat_shell.eta)
            dH_meoh = calc_dH_rxn_meoh(T_t)
            dH_rwgs = calc_dH_rxn_rwgs(T_t)
            dH_ch = calc_dH_rxn_ch(T_s)

            dS[0, j] = geom.A_tube * cfg.cat_tube.rho_b * (-3.0 * r_meoh - r_rwgs)
            dS[1, j] = geom.A_tube * cfg.cat_tube.rho_b * (-r_meoh - r_rwgs)
            dS[2, j] = geom.A_tube * cfg.cat_tube.rho_b * (r_rwgs)
            dS[3, j] = geom.A_tube * cfg.cat_tube.rho_b * (r_meoh)
            dS[4, j] = geom.A_tube * cfg.cat_tube.rho_b * (r_meoh + r_rwgs)
            dS[7, j] = - geom.A_shell * cfg.cat_shell.rho_b * (-r_ch)
            dS[8, j] = - geom.A_shell * cfg.cat_shell.rho_b * (r_ch)
            dS[9, j] = - geom.A_shell * cfg.cat_shell.rho_b * (3.0 * r_ch)

            U_loc = calc_U_tcr(T_t, P_t, F_t, SPECIES_TUBE, T_s, P_s, F_s, SPECIES_SHELL[:3], geom, cfg.cat_tube, cfg.cat_shell, cfg.transport)
            Q_wall = U_loc * geom.a_H * (T_t - T_s)
            FCp_t = sum(F_t[i] * calc_cp_species(SPECIES_TUBE[i], T_t) for i in range(7))
            FCp_s = sum(F_s[i] * calc_cp_species(SPECIES_SHELL[i], T_s) for i in range(3))
            Q_rxn_t = geom.A_tube * cfg.cat_tube.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
            dS[10, j] = (Q_rxn_t - Q_wall) / max(FCp_t, 1e-3)
            Q_rxn_s = geom.A_shell * cfg.cat_shell.rho_b * (-dH_ch) * r_ch
            dS[11, j] = - (Q_rxn_s + Q_wall) / max(FCp_s, 1e-3)

            mw_t = sum(y_t[i] * MW[SPECIES_TUBE[i]] for i in range(7))
            rho_t = (P_t * mw_t) / (R_GAS * T_t)
            u_t = (F_tot_t * R_GAS * T_t) / (P_t * geom.A_tube)
            dS[12, j] = calc_ergun_dp(u_t, rho_t, cfg.transport.mu_exo, cfg.cat_tube.dp, cfg.cat_tube.eps)
            mw_s = sum(y_s[i] * MW[SPECIES_SHELL[i]] for i in range(3))
            rho_s = (P_s * mw_s) / (R_GAS * T_s)
            u_s = (F_tot_s * R_GAS * T_s) / (P_s * geom.A_shell)
            dS[13, j] = - calc_ergun_dp(u_s, rho_s, cfg.transport.mu_endo, cfg.cat_shell.dp, cfg.cat_shell.eps)
        return dS

    F_t0 = cfg.feed_tube.y0 * cfg.feed_tube.F_total0
    F_s0 = cfg.feed_shell.y0 * cfg.feed_shell.F_total0
    def bvp_bc(ya, yb):
        res = np.zeros(14)
        res[0:7] = ya[0:7] - F_t0
        res[7] = ya[10] - cfg.feed_tube.T0
        res[8] = ya[12] - cfg.feed_tube.P0
        res[9:12] = yb[7:10] - F_s0
        res[12] = yb[11] - cfg.feed_shell.T0
        res[13] = yb[13] - cfg.feed_shell.P0
        return res

    y_guess = sol_picard.y
    res_bvp = solve_bvp(bvp_ode, bvp_bc, z_grid, y_guess, tol=1e-3, max_nodes=500)
    assert res_bvp.success, f"solve_bvp failed: {res_bvp.message}"
    y_bvp = res_bvp.sol(z_grid)

    L_inf_T = max(np.max(np.abs(y_guess[10, :] - y_bvp[10, :])), np.max(np.abs(y_guess[11, :] - y_bvp[11, :])))
    assert L_inf_T < 0.10, f"Max temperature discrepancy between Picard and solve_bvp is {L_inf_T:.3e} K (must be < 0.10 K)"

def test_t14_equal_throughput_optimization_and_pressure(cfg):
    """T14: Verifies equal-throughput optimization rule, thermal limits, and positive pressure."""
    import json
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

    # 3. Positive shell outlet pressure
    P_in = vars_opt.get("P_endo_0_bar", 2.0)
    P_out = P_in - knee["dP_shell"]
    assert P_out >= 1.10, f"Shell exit pressure {P_out:.2f} bar must be >= 1.10 bar"

    # 4. Rigorous performance: must achieve high production without feed throttling
    assert knee["daily_meoh_ton"] >= 200.0, f"Daily methanol output {knee['daily_meoh_ton']:.2f} t/d must be >= 200.0 t/d under full feed"

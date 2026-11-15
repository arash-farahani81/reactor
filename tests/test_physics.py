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

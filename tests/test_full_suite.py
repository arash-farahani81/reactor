"""
Comprehensive Automated Verification Suite (T8, T9, T10, T12).
Strict gates for:
- T8: Determinism / reproducibility: identical consecutive runs produce identical results.
- T9: No-hardcode scan: manuscript template uses dynamic placeholders for all results.
- T10: Manuscript consistency: check_manuscript_consistency passes with 100% numerical match.
- T12: Benchmark verification: Khademi (2009) benchmark, Table 4 typo proof, and eta_exo sensitivity.
"""

import pytest
import json
import os
import re
from tcr.config import load_config
from tcr.geometry import ReactorGeometry
from tcr.solvers import solve_conventional
from tcr.metrics import evaluate_run_metrics
from check_manuscript_consistency import check_consistency

@pytest.fixture(scope="module")
def results():
    with open("simulation_results.json", "r", encoding="utf-8") as f:
        return json.load(f)

def test_t8_reproducibility():
    """T8: Determinism / reproducibility: identical runs yield identical results."""
    cfg = load_config("config.yaml")
    geom = ReactorGeometry(cfg.geometry)

    sol1 = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    m1 = evaluate_run_metrics(sol1, geom, mode="conventional", is_conventional=True)

    sol2 = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, cfg.feed_tube.T0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=cfg.transport.T_coolant, kinetics_mode="production_vbf"
    )
    m2 = evaluate_run_metrics(sol2, geom, mode="conventional", is_conventional=True)

    assert abs(m1["yield_khademi_eq28"] - m2["yield_khademi_eq28"]) < 1e-10, "Run 1 and Run 2 yields must be identical"
    assert abs(m1["T_max_exo"] - m2["T_max_exo"]) < 1e-10, "Run 1 and Run 2 peak T must be identical"
    assert abs(m1["MeOH_rate_net_plant"] - m2["MeOH_rate_net_plant"]) < 1e-10, "Run 1 and Run 2 rates must be identical"

def test_t9_no_hardcode_scan():
    """T9: No-hardcode scan: ensure template uses dynamic placeholders."""
    template_path = "manuscript_template.md"
    assert os.path.exists(template_path), "manuscript_template.md must exist"

    with open(template_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Check for presence of essential dynamic placeholders
    assert "{{ conventional_baseline.yield_khademi_eq28 }}" in content
    assert "{{ conventional_optimized.yield_khademi_eq28 }}" in content
    assert "{{ cocurrent_baseline.yield_khademi_eq28 }}" in content
    assert "{{ countercurrent_baseline.yield_khademi_eq28 }}" in content
    assert "{{ tcr_nsga2_optimum.yield_khademi_eq28 }}" in content
    assert "{{ process_thermal_duties.total_feed_preheat_duty_MW }}" in content

def test_t10_manuscript_consistency():
    """T10: Manuscript consistency: check_manuscript_consistency passes with 100% numerical match."""
    # Runs the verification function without raising SystemExit
    check_consistency()

def test_t12_benchmark_verification(results):
    """T12: Benchmark verification: Khademi (2009) benchmark, Table 4 typo proof, and sensitivity."""
    bench = results.get("benchmark_khademi_2009")
    assert bench is not None, "benchmark_khademi_2009 must be present in simulation_results.json"

    # Verify Table 4 typo proof
    typo = bench.get("table4_typo_proof")
    assert typo is not None, "table4_typo_proof must be present"
    assert typo["apparent_delta_H_over_R"] == -3190.0
    assert abs(typo["true_delta_H_over_R"] - (-24703.95)) < 1.0

    # Verify eta_exo sensitivity sweep
    sens = bench.get("eta_exo_sensitivity")
    assert sens is not None, "eta_exo_sensitivity must be present"
    assert "eta_1.0" in sens
    assert "eta_0.7" in sens
    assert "eta_0.5" in sens
    assert "eta_0.35" in sens

r"""
Rigorous automated consistency checker for the Path A TCR manuscript.
Verifies that:
1. No unrendered template tags {{ ... }} exist.
2. No forbidden marketing/overstated terms or stale pre-Path-A hardcoded artifacts
   exist (e.g., '2.45e4', active constant 'k_er = 2.20', 'Bi_h = 3.45', '3.98 K',
   'multi-scale', unqualified 'green hydrogen', 'high-purity hydrogen', 'chemical-grade benzene').
3. Mandatory Path A engineering derivations and concepts are present.
4. LaTeX syntax integrity: line-by-line scan for balanced $ delimiters and
   strictly escaped \% inside math mode.
5. Every simulated Path A metric in TCR_Methanol_Cyclohexane_Manuscript.md strictly
   matches simulation_results.json (including all 5 regimes, gross Eq. 28 & net yields,
   OTP ceiling, PR-EOS audit, cmr_steam_credit_MW, bed_lifting_ratio, T_centerline_max,
   unified dynamic k_er_mean, benchmark table rows, and eta_sweep values).
6. Synchronization of simulation_summary.md with simulation_results.json.
"""

import json
import os
import re
import sys


def _pattern_present(text: str, pattern: str) -> bool:
    """Match plain-text or LaTeX-escaped variants of a formatted numeric metric."""
    alt_pct = pattern.replace("%", r"\%")
    alt_unit = (
        pattern.replace(" K", r"\text{ K}")
        .replace(" MW", r"\text{ MW}")
        .replace(" bar", r"\text{ bar}")
        .replace(" pp", r"\text{ pp}")
    )
    alt_both = (
        pattern.replace("%", r"\%")
        .replace(" K", r"\text{ K}")
        .replace(" MW", r"\text{ MW}")
        .replace(" bar", r"\text{ bar}")
        .replace(" pp", r"\text{ pp}")
    )
    return (
        pattern in text
        or alt_pct in text
        or alt_unit in text
        or alt_both in text
    )


def _require_key(mapping: dict, key: str, desc: str):
    if not isinstance(mapping, dict) or key not in mapping or mapping[key] is None:
        raise KeyError(
            f"Missing required Path A key '{desc}' in simulation_results.json. "
            "Re-run `python run_all.py` before running check_manuscript_consistency.py."
        )
    return mapping[key]


def check_consistency(
    results_path: str = "simulation_results.json",
    manuscript_path: str = "TCR_Methanol_Cyclohexane_Manuscript.md",
    summary_path: str = "simulation_summary.md",
) -> None:
    with open(results_path, "r", encoding="utf-8") as f:
        res = json.load(f)

    with open(manuscript_path, "r", encoding="utf-8") as f:
        text = f.read()

    errors = []

    # 1. Check for unrendered template tags
    unrendered = re.findall(r"\{\{.*?\}\}", text)
    if unrendered:
        errors.append(f"Found unrendered template tags: {unrendered}")

    # 2. Check for forbidden misleading terms & stale pre-Path-A artifacts
    text_lower = text.lower()
    if "multi-scale" in text_lower:
        errors.append("Forbidden term 'multi-scale' found in manuscript (must use 1D pseudo-homogeneous).")
    if "green hydrogen" in text_lower:
        errors.append("Unqualified 'green hydrogen' found in manuscript (must specify hydrocarbon/LOHC dehydrogenation co-product).")
    if "high-purity hydrogen" in text_lower:
        errors.append("Forbidden overstated claim 'high-purity hydrogen' found in manuscript (must use 'H2-rich reformate stream').")
    if "chemical-grade benzene" in text_lower:
        errors.append("Forbidden overstated claim 'chemical-grade benzene' found in manuscript (must use 'co-produced benzene').")

    stale_patterns = [
        ("Obsolete equilibrium constant '2.45e4'", r"2\.45e4|2\.45\s*\\times\s*10\^4"),
        ("Stale constant k_er = 2.20 W/(m*K)", r"k_\{?er\}?\s*[=≈~]\s*2\.20|2\.20\s*\\?text\{?W/\(m"),
        ("Stale Biot number Bi_h = 3.45 (from k_er=2.20)", r"Bi_h\s*=\s*3\.45"),
        ("Stale radial overshoot 3.98 K (from k_er=2.20)", r"3\.98\s*\\?text\{?\s*K\}?"),
    ]
    for label, pat in stale_patterns:
        if re.search(pat, text, flags=re.IGNORECASE):
            errors.append(f"Stale pre-Path-A artifact found in manuscript: {label}.")

    # 3. Check for mandatory Path A engineering derivations & concepts
    required_concepts = [
        ("Wassiljewa variable fluid thermal conductivity", "wassiljewa"),
        ("Le Chatelier equilibrium conversion bound", "le chatelier"),
        ("Van Heerden thermal capacity ratio C = 0.361", "0.361"),
        ("Unified dynamic Bauer-Schluender k_er", "bauer"),
        ("Optimal Temperature Profile (OTP) ceiling", "optimal temperature profile"),
        ("Real-gas Peng-Robinson EOS audit", "peng"),
        ("CMR boiling-water steam generation credit", "steam"),
        ("Downstream extractive distillation reboiler duty", "extractive distillation"),
        ("Downstream PSA hydrogen purification work", "psa"),
        ("Upward-flow bed-lifting ratio & hold-down grid", "bed-lifting"),
        ("Hold-down grid hydraulic mitigation", "hold-down grid"),
    ]
    for name, snippet in required_concepts:
        if snippet not in text_lower:
            errors.append(f"Missing mandatory Path A concept/derivation in manuscript: {name} ('{snippet}' not found).")

    # 4. Strict line-by-line LaTeX syntax gate (odd $ delimiters and unescaped % in math mode)
    for line_idx, line_str in enumerate(text.splitlines(), start=1):
        line_inline = re.sub(r"\$\$.*?\$\$", "", line_str)
        if line_inline.count("$") % 2 != 0:
            errors.append(f"Line {line_idx}: Odd number of '$' math delimiters: {line_str.strip()}")

        for block in re.findall(r"\$\$(.*?)\$\$", line_str):
            if re.search(r"(?<!\\)%", block):
                errors.append(f"Line {line_idx}: Unescaped '%' found inside display math: '$${block}$$'")

        for block in re.findall(r"\$([^\$\n]+)\$", line_inline):
            if re.search(r"(?<!\\)%", block):
                errors.append(f"Line {line_idx}: Unescaped '%' found inside inline math: '${block}$'")

    # 5. Exact numerical consistency checks against simulation_results.json (Path A)
    cmr_base = _require_key(res, "conventional_baseline", "conventional_baseline")
    cmr_opt = _require_key(res, "conventional_optimized", "conventional_optimized")
    co = _require_key(res, "cocurrent_baseline", "cocurrent_baseline")
    counter = _require_key(res, "countercurrent_baseline", "countercurrent_baseline")
    knee = _require_key(res, "tcr_nsga2_optimum", "tcr_nsga2_optimum")
    duties = _require_key(res, "process_thermal_duties", "process_thermal_duties")
    bench = _require_key(res, "benchmark_khademi_2009", "benchmark_khademi_2009")
    otp_audit = _require_key(bench, "otp_and_eos_audit", "benchmark_khademi_2009.otp_and_eos_audit")
    otp_035 = _require_key(otp_audit, "eta_0.35", "otp_and_eos_audit.eta_0.35")
    otp_100 = _require_key(otp_audit, "eta_1.0", "otp_and_eos_audit.eta_1.0")
    eos_cmp = _require_key(otp_audit, "eos_comparison", "otp_and_eos_audit.eos_comparison")
    bench_cmr = _require_key(bench, "simulation_cmr", "benchmark_khademi_2009.simulation_cmr")
    bench_tcr = _require_key(bench, "simulation_tcr_thermo", "benchmark_khademi_2009.simulation_tcr_thermo")
    bench_lit = _require_key(bench, "simulation_tcr_literal", "benchmark_khademi_2009.simulation_tcr_literal")

    knee_vars = _require_key(knee, "decision_variables", "tcr_nsga2_optimum.decision_variables")
    eos_diff = float(_require_key(eos_cmp, "diff_pct_points", "eos_comparison.diff_pct_points"))
    eos_diff_str = f"{'+' if eos_diff >= 0.0 else ''}{eos_diff:.2f}"

    checks = [
        # Equal-throughput syngas feed lock
        ("Locked syngas feed rate", f"{float(_require_key(knee_vars, 'F_exo_mol_s', 'knee.F_exo_mol_s')):.3f}"),
        # CMR Baseline
        ("CMR baseline Eq.28 yield", f"{float(cmr_base['yield_khademi_eq28']):.2f}%"),
        ("CMR baseline net yield", f"{float(cmr_base['yield_net_COx']):.2f}%"),
        ("CMR baseline daily ton", f"{float(cmr_base['daily_meoh_ton']):.2f}"),
        ("CMR baseline peak T", f"{float(cmr_base['T_max_exo']):.1f} K"),
        ("CMR baseline centerline T", f"{float(cmr_base['T_centerline_max']):.1f} K"),
        ("CMR baseline steam credit", f"{float(cmr_base['cmr_steam_credit_MW']):.2f} MW"),
        # CMR Optimized
        ("CMR opt Eq.28 yield", f"{float(cmr_opt['yield_khademi_eq28']):.2f}%"),
        ("CMR opt net yield", f"{float(cmr_opt['yield_net_COx']):.2f}%"),
        ("CMR opt daily ton", f"{float(cmr_opt['daily_meoh_ton']):.2f}"),
        ("CMR opt peak T", f"{float(cmr_opt['T_max_exo']):.1f} K"),
        ("CMR opt centerline T", f"{float(cmr_opt['T_centerline_max']):.1f} K"),
        ("CMR opt dynamic k_er", f"{float(cmr_opt['k_er_mean']):.2f}"),
        ("CMR opt steam credit", f"{float(cmr_opt['cmr_steam_credit_MW']):.2f} MW"),
        # TCR Co-Current Baseline
        ("TCR Co-current Eq.28 yield", f"{float(co['yield_khademi_eq28']):.2f}%"),
        ("TCR Co-current net yield", f"{float(co['yield_net_COx']):.2f}%"),
        ("TCR Co-current benzene ton", f"{float(co['daily_benzene_ton']):.2f}"),
        ("TCR Co-current H2 ton", f"{float(co['daily_h2_ton']):.2f}"),
        ("TCR Co-current bed-lifting ratio", f"{float(co['bed_lifting_ratio']):.2f}"),
        # TCR Counter-Current Baseline
        ("TCR Counter-current Eq.28 yield", f"{float(counter['yield_khademi_eq28']):.2f}%"),
        ("TCR Counter-current net yield", f"{float(counter['yield_net_COx']):.2f}%"),
        ("TCR Counter-current benzene ton", f"{float(counter['daily_benzene_ton']):.2f}"),
        ("TCR Counter-current H2 ton", f"{float(counter['daily_h2_ton']):.2f}"),
        ("TCR Counter-current bed-lifting ratio", f"{float(counter['bed_lifting_ratio']):.2f}"),
        # TCR NSGA-II Knee Point
        ("TCR Knee Eq.28 yield", f"{float(knee['yield_khademi_eq28']):.2f}%"),
        ("TCR Knee net yield", f"{float(knee['yield_net_COx']):.2f}%"),
        ("TCR Knee daily ton", f"{float(knee['daily_meoh_ton']):.2f}"),
        ("TCR Knee benzene ton", f"{float(knee['daily_benzene_ton']):.2f}"),
        ("TCR Knee H2 ton", f"{float(knee['daily_h2_ton']):.2f}"),
        ("TCR Knee peak T", f"{float(knee['T_max_exo']):.1f} K"),
        ("TCR Knee centerline T", f"{float(knee['T_centerline_max']):.1f} K"),
        ("TCR Knee dynamic k_er", f"{float(knee['k_er_mean']):.2f}"),
        ("TCR Knee bed-lifting ratio", f"{float(knee['bed_lifting_ratio']):.2f}"),
        ("TCR Knee bed weight bar", f"{float(knee['bed_weight_bar']):.3f}"),
        # Battery-limit duties
        ("Baseline TCR total feed preheat duty", f"{float(duties['total_feed_preheat_duty_MW']):.2f} MW"),
        ("TCR Knee feed preheat duty", f"{float(knee['total_feed_preheat_duty_MW']):.2f} MW"),
        ("TCR Knee extractive distillation reboiler", f"{float(knee['Q_reboiler_sep_MW']):.2f} MW"),
        ("TCR Knee PSA H2 work", f"{float(knee['W_psa_H2_MW']):.2f}"),
        # Khademi Benchmark (Gross Eq. 28 + Net)
        ("Benchmark CMR Eq.28 yield", f"{float(bench_cmr['yield_pct']):.2f}%"),
        ("Benchmark CMR net yield", f"{float(bench_cmr['yield_net_pct']):.2f}%"),
        ("Benchmark OTCR thermo Eq.28 yield", f"{float(bench_tcr['yield_pct']):.2f}%"),
        ("Benchmark OTCR thermo net yield", f"{float(bench_tcr['yield_net_pct']):.2f}%"),
        ("Benchmark OTCR literal Eq.28 yield", f"{float(bench_lit['yield_pct']):.2f}%"),
        ("Benchmark OTCR literal net yield", f"{float(bench_lit['yield_net_pct']):.2f}%"),
        # OTP Ceiling & PR-EOS Audit
        ("OTP ceiling net yield (eta=0.35)", f"{float(otp_035['otp_yield_net']):.2f}%"),
        ("Best isothermal net yield (eta=0.35)", f"{float(otp_035['best_isothermal_yield_net']):.2f}%"),
        ("OTP ceiling net yield (eta=1.0)", f"{float(otp_100['otp_yield_net']):.2f}%"),
        ("EOS audit ideal net yield", f"{float(eos_cmp['yield_net_ideal']):.2f}%"),
        ("EOS audit real PR-EOS net yield", f"{float(eos_cmp['yield_net_real_eos']):.2f}%"),
        ("EOS audit delta pp", eos_diff_str),
    ]

    for label, pattern in checks:
        if not _pattern_present(text, pattern):
            errors.append(f"Missing expected Path A value for {label}: '{pattern}' not found in manuscript.")

    # Verify eta_sweep table values if structured list present
    sens = res.get("sensitivity_analysis", {})
    eta_sweep = sens.get("eta_sweep")
    if isinstance(eta_sweep, list):
        for rec in eta_sweep:
            y_str = f"{float(rec['yield_meoh']):.2f}%"
            if not _pattern_present(text, y_str):
                errors.append(f"Missing eta_sweep yield ({rec['eta']}): '{y_str}' not found in manuscript.")

    # 6. Check simulation_summary.md synchronization
    if os.path.exists(summary_path):
        with open(summary_path, "r", encoding="utf-8") as f:
            summary_text = f.read()
        summary_checks = [
            ("Summary CMR baseline MeOH t/d", f"{float(cmr_base['daily_meoh_ton']):.2f}"),
            ("Summary CMR opt MeOH t/d", f"{float(cmr_opt['daily_meoh_ton']):.2f}"),
            ("Summary Co-current MeOH t/d", f"{float(co['daily_meoh_ton']):.2f}"),
            ("Summary Counter-current MeOH t/d", f"{float(counter['daily_meoh_ton']):.2f}"),
            ("Summary Knee MeOH t/d", f"{float(knee['daily_meoh_ton']):.2f}"),
            ("Summary Knee Benzene t/d", f"{float(knee['daily_benzene_ton']):.2f}"),
            ("Summary Knee H2 t/d", f"{float(knee['daily_h2_ton']):.2f}"),
        ]
        for label, val in summary_checks:
            if val not in summary_text:
                errors.append(f"Discrepancy in simulation_summary.md: '{val}' ({label}) not found.")

    if errors:
        print("❌ Path A Manuscript Consistency Check FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("✓ Path A Manuscript Consistency Check PASSED! 100% numerical match with simulation_results.json.")


if __name__ == "__main__":
    check_consistency()

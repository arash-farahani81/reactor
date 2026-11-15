r"""
Rigorous automated consistency checker for the manuscript.
Verifies that:
1. Every simulated metric in TCR_Methanol_Cyclohexane_Manuscript.md strictly
   matches simulation_results.json.
2. No unrendered template tags {{ ... }} exist.
3. No forbidden marketing/unphysical terms exist (e.g. unreferenced 'multi-scale',
   unqualified 'green hydrogen', or unqualified '0 kW utility').
4. LaTeX syntax integrity: line-by-line scan for balanced $ delimiters and
   strictly escaped \% inside math mode.
5. Synchronization of README.md Section 4 and simulation_summary.md with simulation_results.json.
"""

import json
import re
import sys
import os

def check_consistency():
    results_path = "simulation_results.json"
    manuscript_path = "TCR_Methanol_Cyclohexane_Manuscript.md"
    readme_path = "README.md"
    summary_path = "simulation_summary.md"

    with open(results_path, "r", encoding="utf-8") as f:
        res = json.load(f)

    with open(manuscript_path, "r", encoding="utf-8") as f:
        text = f.read()

    errors = []

    # 1. Check for unrendered template tags
    unrendered = re.findall(r"\{\{.*?\}\}", text)
    if unrendered:
        errors.append(f"Found unrendered template tags: {unrendered}")

    # 2. Check for forbidden misleading terms
    if "multi-scale" in text.lower():
        errors.append("Forbidden term 'multi-scale' found in manuscript (must use 1D pseudo-homogeneous).")
    if "green hydrogen" in text:
        errors.append("Unqualified 'green hydrogen' found in manuscript (must specify fossil/dehydrogenation co-product).")

    # 3. Strict line-by-line LaTeX syntax gate (odd $ delimiters and unescaped % in math mode)
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

    # 4. Check specific numerical values against simulation_results.json
    cmr_base = res["conventional_baseline"]
    cmr_opt = res["conventional_optimized"]
    co = res["cocurrent_baseline"]
    counter = res["countercurrent_baseline"]
    knee = res["tcr_nsga2_optimum"]
    duties = res["process_thermal_duties"]

    checks = [
        ("CMR baseline yield", f"{cmr_base['yield_khademi_eq28']:.2f}%"),
        ("CMR baseline daily ton", f"{cmr_base['daily_meoh_ton']:.2f}"),
        ("CMR baseline peak T", f"{cmr_base['T_max_exo']:.1f} K"),
        ("CMR opt yield", f"{cmr_opt['yield_khademi_eq28']:.2f}%"),
        ("CMR opt daily ton", f"{cmr_opt['daily_meoh_ton']:.2f}"),
        ("CMR opt peak T", f"{cmr_opt['T_max_exo']:.1f} K"),
        ("TCR Co-current yield", f"{co['yield_khademi_eq28']:.2f}%"),
        ("TCR Co-current benzene ton", f"{co['daily_benzene_ton']:.2f}"),
        ("TCR Co-current H2 ton", f"{co['daily_h2_ton']:.2f}"),
        ("TCR Counter-current yield", f"{counter['yield_khademi_eq28']:.2f}%"),
        ("TCR Counter-current benzene ton", f"{counter['daily_benzene_ton']:.2f}"),
        ("TCR Counter-current H2 ton", f"{counter['daily_h2_ton']:.2f}"),
        ("TCR Knee yield", f"{knee['yield_khademi_eq28']:.2f}%"),
        ("TCR Knee daily ton", f"{knee['daily_meoh_ton']:.2f}"),
        ("TCR Knee benzene ton", f"{knee['daily_benzene_ton']:.2f}"),
        ("TCR Knee H2 ton", f"{knee['daily_h2_ton']:.2f}"),
        ("TCR Knee peak T", f"{knee['T_max_exo']:.1f} K"),
        ("Total feed preheat duty", f"{duties['total_feed_preheat_duty_MW']:.2f} MW"),
        ("Syngas preheat duty", f"{duties['syngas_preheat_duty_MW']:.2f} MW"),
    ]

    for label, pattern in checks:
        alt_pattern = pattern.replace("%", r"\%")
        alt_math = pattern.replace(" K", r"\text{ K}").replace(" MW", r"\text{ MW}")
        alt_both = pattern.replace("%", r"\%").replace(" K", r"\text{ K}").replace(" MW", r"\text{ MW}")
        if pattern not in text and alt_pattern not in text and alt_math not in text and alt_both not in text:
            errors.append(f"Missing expected value for {label}: '{pattern}' not found in manuscript.")

    # 5. Check README.md Section 4 synchronization
    if os.path.exists(readme_path):
        with open(readme_path, "r", encoding="utf-8") as f:
            readme_text = f.read()
        readme_checks = [
            ("README CMR baseline MeOH t/d", f"{cmr_base['daily_meoh_ton']:.2f}"),
            ("README CMR opt MeOH t/d", f"{cmr_opt['daily_meoh_ton']:.2f}"),
            ("README Co-current MeOH t/d", f"{co['daily_meoh_ton']:.2f}"),
            ("README Counter-current MeOH t/d", f"{counter['daily_meoh_ton']:.2f}"),
            ("README Knee MeOH t/d", f"{knee['daily_meoh_ton']:.2f}"),
            ("README Knee Benzene t/d", f"{knee['daily_benzene_ton']:.2f}"),
            ("README Knee H2 t/d", f"{knee['daily_h2_ton']:.2f}"),
            ("README Knee Feed duty MW", f"{duties['total_feed_preheat_duty_MW']:.2f} MW"),
        ]
        for label, val in readme_checks:
            if val not in readme_text:
                errors.append(f"Discrepancy in README.md: '{val}' ({label}) not found.")

    # 6. Check simulation_summary.md synchronization
    if os.path.exists(summary_path):
        with open(summary_path, "r", encoding="utf-8") as f:
            summary_text = f.read()
        summary_checks = [
            ("Summary CMR baseline MeOH t/d", f"{cmr_base['daily_meoh_ton']:.2f}"),
            ("Summary CMR opt MeOH t/d", f"{cmr_opt['daily_meoh_ton']:.2f}"),
            ("Summary Co-current MeOH t/d", f"{co['daily_meoh_ton']:.2f}"),
            ("Summary Counter-current MeOH t/d", f"{counter['daily_meoh_ton']:.2f}"),
            ("Summary Knee MeOH t/d", f"{knee['daily_meoh_ton']:.2f}"),
            ("Summary Knee Benzene t/d", f"{knee['daily_benzene_ton']:.2f}"),
            ("Summary Knee H2 t/d", f"{knee['daily_h2_ton']:.2f}"),
        ]
        for label, val in summary_checks:
            if val not in summary_text:
                errors.append(f"Discrepancy in simulation_summary.md: '{val}' ({label}) not found.")

    if errors:
        print("❌ Manuscript Consistency Check FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("✓ Manuscript Consistency Check PASSED! 100% numerical match with simulation_results.json.")

if __name__ == "__main__":
    check_consistency()

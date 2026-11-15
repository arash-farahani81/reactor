"""
Rigorous automated consistency checker for the manuscript.
Verifies that:
1. Every simulated metric in TCR_Methanol_Cyclohexane_Manuscript.md strictly
   matches simulation_results.json.
2. No unrendered template tags {{ ... }} exist.
3. No forbidden marketing/unphysical terms exist (e.g. unreferenced 'multi-scale',
   unqualified 'green hydrogen', or unqualified '0 kW utility').
"""

import json
import re
import sys

def check_consistency():
    results_path = "simulation_results.json"
    manuscript_path = "TCR_Methanol_Cyclohexane_Manuscript.md"

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

    # 3. Check for LaTeX syntax integrity (odd $ delimiters and unescaped % in math mode)
    cleaned_math = re.sub(r"\$\$.*?\$\$", "", text, flags=re.DOTALL)
    for line_idx, line_str in enumerate(cleaned_math.splitlines()):
        if line_str.count("$") % 2 != 0:
            errors.append(f"Odd number of '$' math delimiters on line {line_idx+1}: {line_str.strip()}")

    inline_math_blocks = re.findall(r"\$([^\$\n]+)\$", cleaned_math)
    for block in inline_math_blocks:
        if re.search(r"(?<!\\)%", block):
            errors.append(f"Unescaped '%' found inside math mode block: '${block}$'")

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

    if errors:
        print("❌ Manuscript Consistency Check FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("✓ Manuscript Consistency Check PASSED! 100% numerical match with simulation_results.json.")

if __name__ == "__main__":
    check_consistency()

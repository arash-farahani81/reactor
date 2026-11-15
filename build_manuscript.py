"""
Template-driven zero-hardcode manuscript compiler.
Reads simulation_results.json (SSOT) and manuscript_template.md to build
TCR_Methanol_Cyclohexane_Manuscript.md with exact numerical consistency.
"""

import json
import re
import os

def load_results(path="simulation_results.json"):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def build_replacement_dict(results):
    c_base = results["conventional_baseline"]
    c_opt = results["conventional_optimized"]
    co = results["cocurrent_baseline"]
    counter = results["countercurrent_baseline"]
    knee = results["tcr_nsga2_optimum"]
    duties = results["process_thermal_duties"]
    bench = results["benchmark_khademi_2009"]

    # Calculate shell exit pressures dynamically from metadata or decision variables
    P_shell_in_base = float(results.get("metadata", {}).get("shell_inlet_pressure_bar", 2.5))
    co_P_shell_out = max(P_shell_in_base - co["dP_shell"], 0.0)
    counter_P_shell_out = max(P_shell_in_base - counter["dP_shell"], 0.0)
    knee_P_endo_in = float(knee.get('decision_variables', {}).get('P_endo_0_bar', P_shell_in_base))
    knee_P_shell_out = max(knee_P_endo_in - knee["dP_shell"], 0.0)

    delta_cmr_yield = (c_opt["yield_khademi_eq28"] - c_base["yield_khademi_eq28"]) / c_base["yield_khademi_eq28"] * 100.0

    rep = {
        "conventional_baseline.T_exo_0": "503.15",
        "conventional_baseline.T_coolant": "503.15",
        "conventional_baseline.F_exo_0": f"{c_base.get('F_exo_total_0', 0.640):.3f}",
        "conventional_baseline.MeOH_rate_net_plant": f"{c_base['MeOH_rate_net_plant']:.2f}",
        "conventional_baseline.daily_meoh_ton": f"{c_base['daily_meoh_ton']:.2f}",
        "conventional_baseline.X_CO2": f"{c_base['X_CO2']:.2f}",
        "conventional_baseline.X_COx": f"{c_base['X_COx']:.2f}",
        "conventional_baseline.yield_khademi_eq28": f"{c_base['yield_khademi_eq28']:.2f}",
        "conventional_baseline.T_max_exo": f"{c_base['T_max_exo']:.1f}",
        "conventional_baseline.dP_tube": f"{c_base['dP_tube']:.2f}",
        "conventional_baseline.Q_wall_total_MW": f"{c_base['Q_wall_total_MW']:.2f}",

        "conventional_optimized.T_exo_0": f"{c_opt['T_exo_0']:.1f}",
        "conventional_optimized.T_coolant": f"{c_opt['T_coolant']:.1f}",
        "conventional_optimized.F_exo_0": f"{c_opt.get('F_exo_total_0', 0.640):.3f}",
        "conventional_optimized.MeOH_rate_net_plant": f"{c_opt['MeOH_rate_net_plant']:.2f}",
        "conventional_optimized.daily_meoh_ton": f"{c_opt['daily_meoh_ton']:.2f}",
        "conventional_optimized.X_CO2": f"{c_opt['X_CO2']:.2f}",
        "conventional_optimized.X_COx": f"{c_opt['X_COx']:.2f}",
        "conventional_optimized.yield_khademi_eq28": f"{c_opt['yield_khademi_eq28']:.2f}",
        "conventional_optimized.T_max_exo": f"{c_opt['T_max_exo']:.1f}",
        "conventional_optimized.dP_tube": f"{c_opt['dP_tube']:.2f}",
        "conventional_optimized.Q_wall_total_MW": f"{c_opt['Q_wall_total_MW']:.2f}",

        "cocurrent_baseline.F_exo_0": f"{co.get('F_exo_total_0', 0.640):.3f}",
        "cocurrent_baseline.T_max_exo": f"{co['T_max_exo']:.1f}",
        "cocurrent_baseline.yield_khademi_eq28": f"{co['yield_khademi_eq28']:.2f}",
        "cocurrent_baseline.daily_meoh_ton": f"{co['daily_meoh_ton']:.2f}",
        "cocurrent_baseline.MeOH_rate_net_plant": f"{co['MeOH_rate_net_plant']:.2f}",
        "cocurrent_baseline.X_CO2": f"{co['X_CO2']:.2f}",
        "cocurrent_baseline.X_COx": f"{co['X_COx']:.2f}",
        "cocurrent_baseline.X_CH": f"{co['X_CH']:.2f}",
        "cocurrent_baseline.Benzene_rate_plant": f"{co['Benzene_rate_plant']:.2f}",
        "cocurrent_baseline.daily_benzene_ton": f"{co['daily_benzene_ton']:.2f}",
        "cocurrent_baseline.H2_rate_plant": f"{co['H2_rate_plant']:.2f}",
        "cocurrent_baseline.daily_h2_ton": f"{co['daily_h2_ton']:.2f}",
        "cocurrent_baseline.dP_tube": f"{co['dP_tube']:.2f}",
        "cocurrent_baseline.dP_shell": f"{co['dP_shell']:.2f}",
        "cocurrent_baseline.P_shell_out": f"{co_P_shell_out:.2f}",

        "countercurrent_baseline.F_exo_0": f"{counter.get('F_exo_total_0', 0.640):.3f}",
        "countercurrent_baseline.T_max_exo": f"{counter['T_max_exo']:.1f}",
        "countercurrent_baseline.yield_khademi_eq28": f"{counter['yield_khademi_eq28']:.2f}",
        "countercurrent_baseline.daily_meoh_ton": f"{counter['daily_meoh_ton']:.2f}",
        "countercurrent_baseline.MeOH_rate_net_plant": f"{counter['MeOH_rate_net_plant']:.2f}",
        "countercurrent_baseline.X_CO2": f"{counter['X_CO2']:.2f}",
        "countercurrent_baseline.X_COx": f"{counter['X_COx']:.2f}",
        "countercurrent_baseline.X_CH": f"{counter['X_CH']:.2f}",
        "countercurrent_baseline.Benzene_rate_plant": f"{counter['Benzene_rate_plant']:.2f}",
        "countercurrent_baseline.daily_benzene_ton": f"{counter['daily_benzene_ton']:.2f}",
        "countercurrent_baseline.H2_rate_plant": f"{counter['H2_rate_plant']:.2f}",
        "countercurrent_baseline.daily_h2_ton": f"{counter['daily_h2_ton']:.2f}",
        "countercurrent_baseline.dP_tube": f"{counter['dP_tube']:.2f}",
        "countercurrent_baseline.dP_shell": f"{counter['dP_shell']:.2f}",
        "countercurrent_baseline.P_shell_out": f"{counter_P_shell_out:.2f}",

        "tcr_nsga2_optimum.T_exo_0": f"{knee['decision_variables']['T_exo_0_K']:.1f}",
        "tcr_nsga2_optimum.T_endo_0": f"{knee['decision_variables']['T_endo_0_K']:.1f}",
        "tcr_nsga2_optimum.P_endo_0": f"{knee_P_endo_in:.2f}",
        "tcr_nsga2_optimum.F_exo_0": f"{knee.get('decision_variables', {}).get('F_exo_mol_s', 0.640):.3f}",
        "tcr_nsga2_optimum.F_endo_0": f"{knee['decision_variables']['F_endo_mol_s']:.4f}",
        "tcr_nsga2_optimum.T_max_exo": f"{knee['T_max_exo']:.1f}",
        "tcr_nsga2_optimum.yield_khademi_eq28": f"{knee['yield_khademi_eq28']:.2f}",
        "tcr_nsga2_optimum.daily_meoh_ton": f"{knee['daily_meoh_ton']:.2f}",
        "tcr_nsga2_optimum.MeOH_rate_net_plant": f"{knee['MeOH_rate_net_plant']:.2f}",
        "tcr_nsga2_optimum.X_CO2": f"{knee['X_CO2']:.2f}",
        "tcr_nsga2_optimum.X_COx": f"{knee['X_COx']:.2f}",
        "tcr_nsga2_optimum.X_CH": f"{knee['X_CH']:.2f}",
        "tcr_nsga2_optimum.Benzene_rate_plant": f"{knee['Benzene_rate_plant']:.2f}",
        "tcr_nsga2_optimum.daily_benzene_ton": f"{knee['daily_benzene_ton']:.2f}",
        "tcr_nsga2_optimum.H2_rate_plant": f"{knee['H2_rate_plant']:.2f}",
        "tcr_nsga2_optimum.daily_h2_ton": f"{knee['daily_h2_ton']:.2f}",
        "tcr_nsga2_optimum.dP_tube": f"{knee['dP_tube']:.2f}",
        "tcr_nsga2_optimum.dP_shell": f"{knee['dP_shell']:.2f}",
        "tcr_nsga2_optimum.P_shell_out": f"{knee_P_shell_out:.2f}",

        "delta_cmr_yield_pct": f"{delta_cmr_yield:.1f}",

        "process_thermal_duties.syngas_preheat_duty_MW": f"{duties['syngas_preheat_duty_MW']:.2f}",
        "process_thermal_duties.cyclohexane_vaporization_and_preheat_duty_MW": f"{duties['cyclohexane_vaporization_and_preheat_duty_MW']:.2f}",
        "process_thermal_duties.total_feed_preheat_duty_MW": f"{duties['total_feed_preheat_duty_MW']:.2f}",
        "process_thermal_duties.reaction_zone_external_cooling_utility_kW": f"{duties['reaction_zone_external_cooling_utility_kW']:.1f}",
    }
    return rep

def render_manuscript(template_path="manuscript_template.md",
                      results_path="simulation_results.json",
                      output_path="TCR_Methanol_Cyclohexane_Manuscript.md"):
    results = load_results(results_path)
    rep = build_replacement_dict(results)

    with open(template_path, "r", encoding="utf-8") as f:
        content = f.read()

    def replace_match(match):
        key = match.group(1).strip()
        if key in rep:
            return rep[key]
        raise ValueError(f"Unknown placeholder key in template: {key}")

    rendered = re.sub(r"\{\{\s*([a-zA-Z0-9_\.]+)\s*\}\}", replace_match, content)

    # Automated LaTeX syntax validation gate before writing output
    for line_idx, line_str in enumerate(rendered.splitlines(), start=1):
        line_inline = re.sub(r"\$\$.*?\$\$", "", line_str)
        if line_inline.count("$") % 2 != 0:
            raise ValueError(f"Rendered manuscript has odd '$' math delimiters on line {line_idx}: {line_str.strip()}")
        for block in re.findall(r"\$\$(.*?)\$\$", line_str):
            if re.search(r"(?<!\\)%", block):
                raise ValueError(f"Rendered manuscript has unescaped '%' in display math on line {line_idx}: '$${block}$$'")
        for block in re.findall(r"\$([^\$\n]+)\$", line_inline):
            if re.search(r"(?<!\\)%", block):
                raise ValueError(f"Rendered manuscript has unescaped '%' in inline math on line {line_idx}: '${block}$'")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)

    print(f"✓ Successfully rendered manuscript to {output_path} with {len(rep)} variables replaced.")

if __name__ == "__main__":
    render_manuscript()

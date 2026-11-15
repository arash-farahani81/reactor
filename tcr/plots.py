"""
Publication-quality plotting module generating all 5 master project figures:
1. comparison_coupled_vs_conventional.png
2. reactor_results.png
3. sensitivity_analysis.png
4. tcr_multiphysics_profiles.png
5. nsga2_pareto_front.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Publication typography & styling
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
mpl.rcParams['axes.edgecolor'] = '#333333'
mpl.rcParams['axes.linewidth'] = 1.1
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 4.5
mpl.rcParams['ytick.major.size'] = 4.5

def plot_multiphysics_profiles(sol_co, sol_ct, sol_cv, geom, out_path="tcr_multiphysics_profiles.png"):
    """Figure 1: Comprehensive 6-panel multiphysics comparison."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=300)
    fig.subplots_adjust(hspace=0.28, wspace=0.28)

    z = sol_co.t
    c_exo_co = '#D62828'
    c_endo_co = '#0077B6'
    c_exo_ct = '#F77F00'
    c_endo_ct = '#023E8A'
    c_conv = '#2B9348'

    # (a) Temperature profiles
    ax = axes[0, 0]
    ax.plot(z, sol_co.y[10, :], '-', color=c_exo_co, lw=2.2, label='Tube (Co-Current)')
    ax.plot(z, sol_co.y[11, :], '--', color=c_endo_co, lw=2.2, label='Shell (Co-Current)')
    ax.plot(z, sol_ct.y[10, :], '-.', color=c_exo_ct, lw=2.2, label='Tube (Counter-Current)')
    ax.plot(z, sol_ct.y[11, :], ':', color=c_endo_ct, lw=2.2, label='Shell (Counter-Current)')
    ax.plot(sol_cv.t, sol_cv.y[7, :], '-', color=c_conv, lw=1.8, label='Conventional (Jacketed)')
    ax.axhline(535.0, color='gray', linestyle=':', lw=1.0, label='Cu Sintering (535 K)')
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Temperature (K)', fontsize=11, fontweight='bold')
    ax.set_title('(a) Thermal Profiles & Hot-Spot Moderation', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8, loc='best', framealpha=0.9)

    # (b) CO2 conversion
    ax = axes[0, 1]
    co2_in = sol_co.y[1, 0]
    X_co = (co2_in - sol_co.y[1, :]) / co2_in * 100.0
    X_ct = (co2_in - sol_ct.y[1, :]) / co2_in * 100.0
    X_cv = (co2_in - sol_cv.y[1, :]) / co2_in * 100.0
    ax.plot(z, X_co, '-', color=c_exo_co, lw=2.2, label='Coupled Co-Current')
    ax.plot(z, X_ct, '-.', color=c_exo_ct, lw=2.2, label='Coupled Counter-Current')
    ax.plot(sol_cv.t, X_cv, '-', color=c_conv, lw=1.8, label='Conventional Reactor')
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('CO₂ Conversion (%)', fontsize=11, fontweight='bold')
    ax.set_title('(b) CO₂ Conversion along Catalyst Bed', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8.5, loc='lower right')

    # (c) Cyclohexane conversion
    ax = axes[0, 2]
    X_ch_co = (sol_co.y[7, 0] - sol_co.y[7, :]) / sol_co.y[7, 0] * 100.0
    X_ch_ct = (sol_ct.y[7, -1] - sol_ct.y[7, :]) / sol_ct.y[7, -1] * 100.0
    ax.plot(z, X_ch_co, '-', color=c_endo_co, lw=2.2, label='Co-Current (Shell)')
    ax.plot(z, X_ch_ct, '-.', color=c_endo_ct, lw=2.2, label='Counter-Current (Shell)')
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('C₆H₁₂ Conversion (%)', fontsize=11, fontweight='bold')
    ax.set_title('(c) Endothermic Dehydrogenation Conversion', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8.5, loc='best')

    # (d) Product formation
    ax = axes[1, 0]
    ax.plot(z, sol_co.y[3, :] * 1000, '-', color='#007F5F', lw=2.2, label='CH₃OH (Co-Current)')
    ax.plot(z, sol_ct.y[3, :] * 1000, '-.', color='#55A630', lw=2.2, label='CH₃OH (Counter-Current)')
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('CH₃OH Flow Rate (mmol/s·tube)', fontsize=11, fontweight='bold', color='#007F5F')
    ax2 = ax.twinx()
    ax2.plot(z, sol_co.y[8, :] * 1000, '--', color='#6A040F', lw=2.0, label='C₆H₆ (Co-Current)')
    ax2.plot(z, sol_ct.y[8, :] * 1000, ':', color='#9D0208', lw=2.0, label='C₆H₆ (Counter-Current)')
    ax2.set_ylabel('C₆H₆ Flow Rate (mmol/s·tube)', fontsize=11, fontweight='bold', color='#6A040F')
    ax.set_title('(d) Product Formation Profiles', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)

    # (e) Pressure profiles
    ax = axes[1, 1]
    ax.plot(z, sol_co.y[12, :] / 1e5, '-', color='#D62828', lw=2.2, label='Tube Side (Methanol)')
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Tube Pressure (bar)', fontsize=11, fontweight='bold', color='#D62828')
    ax3 = ax.twinx()
    ax3.plot(z, sol_co.y[13, :] / 1e5, '--', color='#0077B6', lw=2.2, label='Shell Side (Co-Current)')
    ax3.plot(z, sol_ct.y[13, :] / 1e5, ':', color='#023E8A', lw=2.2, label='Shell Side (Counter-Current)')
    ax3.set_ylabel('Shell Pressure (bar)', fontsize=11, fontweight='bold', color='#0077B6')
    ax.set_title('(e) Packed-Bed Ergun Pressure Profiles', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)

    # (f) Heat flux distribution
    ax = axes[1, 2]
    q_co = sol_co.Q_wall_arr / (geom.a_H * 1000.0) if hasattr(sol_co, 'Q_wall_arr') else np.zeros_like(z)
    q_ct = sol_ct.Q_wall_arr / (geom.a_H * 1000.0) if hasattr(sol_ct, 'Q_wall_arr') else np.zeros_like(z)
    ax.plot(z, q_co, '-', color='#E85D04', lw=2.2, label='Co-Current Heat Flux')
    ax.plot(z, q_ct, '-.', color='#8338EC', lw=2.2, label='Counter-Current Heat Flux')
    ax.axhline(0, color='gray', linestyle='--', lw=1.0)
    ax.set_xlabel('Axial Position, z (m)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Transferred Wall Heat Flux (kW/m²)', fontsize=11, fontweight='bold')
    ax.set_title('(f) Inter-Stream Heat Exchange Flux', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8.5, loc='best')

    plt.suptitle('Thermally Coupled Shell-and-Tube Reactor (TCR) Multiphysics Performance Benchmarks\n'
                 'Exothermic Methanol Synthesis (Tube) Coupled with Cyclohexane Dehydrogenation (Shell)',
                 fontsize=14, fontweight='bold', y=0.99)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_pareto_frontier(opt_res: dict, out_path="nsga2_pareto_front.png"):
    """Figure 2: NSGA-II 2D and 3D Pareto optimal trade-offs."""
    F_pareto = opt_res['F_pareto']
    X_pareto = opt_res['X_pareto']
    knee_idx = opt_res['knee_idx']

    meoh_yields = -F_pareto[:, 0]
    benzene_rates = -F_pareto[:, 1]
    hotspot_rises = F_pareto[:, 2]

    fig = plt.figure(figsize=(15, 6), dpi=300)

    # Subplot 1: 2D Pareto Front
    ax1 = fig.add_subplot(1, 2, 1)
    scatter = ax1.scatter(
        benzene_rates, meoh_yields,
        c=X_pareto[:, 0] - 273.15, cmap='plasma', s=65, edgecolors='k', lw=0.6, alpha=0.9
    )
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Tube Inlet Temperature (°C)', fontsize=10, fontweight='bold')
    ax1.set_xlabel('Benzene Production Rate (mol/s)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Methanol Yield (%)', fontsize=11, fontweight='bold')
    ax1.set_title('(a) Multi-Objective Pareto Frontier (NSGA-II)', fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.4)

    ax1.plot(
        benzene_rates[knee_idx], meoh_yields[knee_idx],
        marker='*', color='#D62828', markersize=16, markeredgecolor='black',
        label=f'Optimal Knee Point\n(Yield: {meoh_yields[knee_idx]:.2f}%, BZ: {benzene_rates[knee_idx]:.2f} mol/s)'
    )
    ax1.legend(fontsize=9, loc='lower left')

    # Subplot 2: 3D Objective Space
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    p3d = ax2.scatter(
        benzene_rates, meoh_yields, hotspot_rises,
        c=X_pareto[:, 3], cmap='viridis', s=60, edgecolors='k', lw=0.5
    )
    cbar2 = plt.colorbar(p3d, ax=ax2, pad=0.1)
    cbar2.set_label('Tube Feed F_exo_0 (mol/s)', fontsize=10, fontweight='bold')
    ax2.set_xlabel('Benzene (mol/s)', fontsize=9, fontweight='bold')
    ax2.set_ylabel('Methanol Yield (%)', fontsize=9, fontweight='bold')
    ax2.set_zlabel('Hotspot Rise (K)', fontsize=9, fontweight='bold')
    ax2.set_title('(b) 3D Objective Space', fontsize=12, fontweight='bold')

    plt.suptitle('Multi-Objective Pareto Optimization of Thermally Coupled Reactor using NSGA-II',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_reactor_results(sol_co, geom, out_path="reactor_results.png"):
    """Figure 3: Co-current reactor standalone performance."""
    z = sol_co.t
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=300)
    fig.suptitle('Thermally Coupled Shell-and-Tube Reactor — Co-Current Baseline Profiles',
                 fontsize=14, fontweight='bold')

    # (a) Temperature
    ax = axes[0, 0]
    ax.plot(z, sol_co.y[10, :], '-', color='#E63946', lw=2.5, label='Tube (Methanol Synthesis)')
    ax.plot(z, sol_co.y[11, :], '-', color='#457B9D', lw=2.5, label='Shell (Dehydrogenation)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Temperature (K)', fontsize=11)
    ax.set_title('(a) Temperature Profiles', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Conversions
    ax = axes[0, 1]
    co2_in = sol_co.y[1, 0]
    X_co2 = (co2_in - sol_co.y[1, :]) / co2_in * 100.0
    ch_in = sol_co.y[7, 0]
    X_ch = (ch_in - sol_co.y[7, :]) / ch_in * 100.0
    ax.plot(z, X_co2, '-', color='#264653', lw=2.5, label='CO₂ Conversion')
    ax.plot(z, X_ch, '--', color='#457B9D', lw=2.5, label='C₆H₁₂ Conversion')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(b) Conversion Profiles', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (c) Molar flow rates
    ax = axes[0, 2]
    ax.plot(z, sol_co.y[3, :] * 1000, '-', color='#2A9D8F', lw=2.5, label='CH₃OH (Tube)')
    ax2 = ax.twinx()
    ax2.plot(z, sol_co.y[8, :] * 1000, '--', color='#457B9D', lw=2.5, label='C₆H₆ (Shell)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('F_CH₃OH (mmol/s·tube)', fontsize=11, color='#2A9D8F')
    ax2.set_ylabel('F_C₆H₆ (mmol/s·tube)', fontsize=11, color='#457B9D')
    ax.set_title('(c) Product Flow Rates', fontsize=12)
    ax.grid(True, alpha=0.3)

    # (d) Pressure
    ax = axes[1, 0]
    ax.plot(z, sol_co.y[12, :] / 1e5, '-', color='#E63946', lw=2.5, label='P_tube')
    ax3 = ax.twinx()
    ax3.plot(z, sol_co.y[13, :] / 1e5, '--', color='#457B9D', lw=2.5, label='P_shell')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('P_tube (bar)', fontsize=11, color='#E63946')
    ax3.set_ylabel('P_shell (bar)', fontsize=11, color='#457B9D')
    ax.set_title('(d) Pressure Profiles', fontsize=12)
    ax.grid(True, alpha=0.3)

    # (e) Syngas flows
    ax = axes[1, 1]
    labels_plot = ['H₂', 'CO₂', 'CO', 'CH₃OH', 'H₂O']
    for i in range(5):
        ax.plot(z, sol_co.y[i, :] * 1000, lw=2.0, label=labels_plot[i])
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Molar Flow (mmol/s·tube)', fontsize=11)
    ax.set_title('(e) Tube-Side Molar Flows', fontsize=12)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    # (f) H2 coproduction
    ax = axes[1, 2]
    h2_prod = (sol_co.y[9, :] - sol_co.y[9, 0]) * 1000
    ax.plot(z, h2_prod, '-', color='#F4A261', lw=2.5)
    ax.fill_between(z, 0, h2_prod, alpha=0.2, color='#F4A261')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Coproduced H₂ (mmol/s·tube)', fontsize=11)
    ax.set_title('(f) Shell-Side H₂ Coproduction', fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_comparison_coupled_vs_conventional(sol_co, sol_cv, out_path="comparison_coupled_vs_conventional.png"):
    """Figure 4: Direct comparison between TCR and Conventional reactor."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=300)
    fig.suptitle('Direct Reactor Comparison: Thermally Coupled (Co-Current) vs. Conventional (Jacketed)',
                 fontsize=14, fontweight='bold')

    # Temperature
    ax = axes[0]
    ax.plot(sol_co.t, sol_co.y[10, :], '-', color='#E63946', lw=2.5, label='Coupled (Tube)')
    ax.plot(sol_co.t, sol_co.y[11, :], '--', color='#457B9D', lw=2.5, label='Coupled (Shell)')
    ax.plot(sol_cv.t, sol_cv.y[7, :], '-.', color='#2A9D8F', lw=2.5, label='Conventional (Jacketed)')
    ax.axhline(503.15, color='gray', linestyle=':', label='Boiling Coolant (503.15 K)')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('Temperature (K)', fontsize=11)
    ax.set_title('(a) Temperature Profiles', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # CO2 Conversion
    ax = axes[1]
    co2_in = sol_co.y[1, 0]
    X_co = (co2_in - sol_co.y[1, :]) / co2_in * 100.0
    X_cv = (co2_in - sol_cv.y[1, :]) / co2_in * 100.0
    ax.plot(sol_co.t, X_co, '-', color='#E63946', lw=2.5, label='Coupled Reactor')
    ax.plot(sol_cv.t, X_cv, '--', color='#2A9D8F', lw=2.5, label='Conventional Reactor')
    ax.set_xlabel('Reactor Length (m)', fontsize=11)
    ax.set_ylabel('CO₂ Conversion (%)', fontsize=11)
    ax.set_title('(b) CO₂ Conversion Comparison', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_sensitivity_analysis(sens_results: dict, out_path="sensitivity_analysis.png"):
    """Figure 5: Parametric sensitivity sweeps."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=300)
    fig.suptitle('Parametric Sensitivity Analysis — Thermally Coupled Reactor', fontsize=14, fontweight='bold')

    # (a) Tube inlet temperature sweep
    ax = axes[0]
    t_sweep = sens_results.get('T_exo_sweep', {})
    if t_sweep:
        ax.plot(t_sweep['T_vals'], t_sweep['X_CO2'], 'o-', color='#E63946', lw=2, label='CO₂ Conv.')
        ax.plot(t_sweep['T_vals'], t_sweep['X_CH'], 's--', color='#457B9D', lw=2, label='C₆H₁₂ Conv.')
    ax.set_xlabel('Tube Inlet Temperature (K)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(a) Effect of Tube Inlet Temperature', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Shell feed flow rate sweep
    ax = axes[1]
    f_sweep = sens_results.get('F_endo_sweep', {})
    if f_sweep:
        ax.plot(np.array(f_sweep['F_vals']) * 1000, f_sweep['X_CO2'], 'o-', color='#E63946', lw=2, label='CO₂ Conv.')
        ax.plot(np.array(f_sweep['F_vals']) * 1000, f_sweep['X_CH'], 's--', color='#457B9D', lw=2, label='C₆H₁₂ Conv.')
    ax.set_xlabel('Shell Feed Flow (mmol/s·tube)', fontsize=11)
    ax.set_ylabel('Conversion (%)', fontsize=11)
    ax.set_title('(b) Effect of Shell Feed Flow Rate', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (c) Catalyst effectiveness factor sweep
    ax = axes[2]
    eta_sweep = sens_results.get('eta_sweep', {})
    if eta_sweep:
        ax.plot(eta_sweep['eta_vals'], eta_sweep['X_CO2'], 'o-', color='#E63946', lw=2, label='CO₂ Conv.')
        ax.plot(eta_sweep['eta_vals'], eta_sweep['yield_meoh'], '^-', color='#2A9D8F', lw=2, label='MeOH Yield')
    ax.set_xlabel('Pellet Effectiveness Factor η_exo (-)', fontsize=11)
    ax.set_ylabel('Metric (%)', fontsize=11)
    ax.set_title('(c) Effect of Catalyst Effectiveness Factor', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

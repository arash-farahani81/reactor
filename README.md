# Thermally Coupled Reactor (TCR) Simulation & Optimization

An industrially-sized, rigorous simulation and multi-objective optimization framework for a **Thermally Coupled Shell-and-Tube Reactor (TCR)** directly coupling:
- **Tube Side (Exothermic):** Methanol Synthesis from syngas ($\text{CO}/\text{CO}_2/\text{H}_2$) over commercial $\text{CuO/ZnO/Al}_2\text{O}_3$ catalyst ($\eta_{\text{exo}} = 0.35$).
- **Shell Side (Endothermic):** Cyclohexane Dehydrogenation to benzene and hydrogen over $\text{Pt/Al}_2\text{O}_3$ catalyst ($\eta_{\text{endo}} = 1.00$).

The repository features single-source-of-truth architecture, zero hardcoded manuscript results, rigorous mass and enthalpy balance closures, under-relaxed Picard iteration for counter-current boundary value problems, and NSGA-II Pareto optimization.

---

## 1. Quick Start

### Installation
Ensure Python 3.10+ is installed:
```bash
pip install -r requirements.txt
```

### Run Full Simulation Pipeline
Executes all reactor models (conventional baseline, conventional optimized, co-current TCR, counter-current TCR), NSGA-II optimization, Khademi et al. (2009) benchmark reproduction, generates all 5 publication-quality figures, updates `simulation_results.json`, and compiles the manuscript:
```bash
python run_all.py
```

### Run Automated Physics Test Suite
Executes all 12 rigorous physics and software tests (mass balances $\le 10^{-9}$, enthalpy closures $\le 0.5\%$, grid convergence, Picard convergence, manuscript consistency):
```bash
pytest tests/ -v
```

### Compile & Verify Manuscript
```bash
python build_manuscript.py
python check_manuscript_consistency.py
```

---

## 2. Repository Architecture

```
reactor/
├── tcr/                           # Core modular simulation package
│   ├── config.py                  # Typed configuration dataclasses
│   ├── thermo.py                  # Cp(T), Hi(T), dHr(T), Kp(T) functions
│   ├── geometry.py                # Exact tube and annular flow geometry
│   ├── transport.py               # Dixon/Li & Finlayson heat transfer & Ergun dP
│   ├── kinetics.py                # VBF 1996 and Khademi/Itoh kinetics
│   ├── solvers.py                 # Radau IIA IVP and Picard BVP solvers
│   ├── metrics.py                 # Element and independent enthalpy balances
│   ├── optimize.py                # CMR optimization and TCR NSGA-II
│   ├── benchmark.py               # Khademi 2009 benchmark & typo proof
│   ├── plots.py                   # Generates all 5 publication figures
│   └── io.py                      # SSOT JSON, summary, and thermal duties
├── tests/                         # Automated verification suite
│   ├── test_physics.py            # Physics gates T1-T7, T11
│   └── test_full_suite.py         # Full verification gates T8-T10, T12
├── legacy/                        # Superseded legacy scripts (quarantined)
│   ├── README.md
│   ├── tcr_advanced.py
│   ├── thermally_coupled_reactor.py
│   └── tcr_simulation_ssot.py
├── baseline/                      # Frozen pre-refactor benchmark outputs
├── figures/                       # High-resolution generated PNG figures
├── config.yaml                    # Single source of truth configuration
├── run_all.py                     # Single entry-point simulation runner
├── build_manuscript.py            # Template-driven manuscript generator
├── check_manuscript_consistency.py# Automated consistency assertion script
├── manuscript_template.md         # Dynamic manuscript template
├── TCR_Methanol_Cyclohexane_Manuscript.md # Compiled research manuscript
├── simulation_results.json        # SSOT simulation outputs
├── simulation_summary.md          # Markdown executive summary
├── REVIEW_LOG.md                  # Detailed phase-by-phase review log
├── FINAL_AUDIT.md                 # Complete D1-D10 defect audit report
└── OPEN_ISSUES.md                 # Physical issue tracking ledger
```

---

## 3. Modeling & Physics Highlights

1. **Reaction Kinetics:**
   - Methanol synthesis modeled via authentic Vanden Bussche & Froment (1996) LHHW equations with temperature-dependent adsorption terms.
   - Cyclohexane dehydrogenation modeled via reversible LHHW kinetics with thermodynamically consistent van 't Hoff equilibrium constant ($\ln K_p = -205400/(RT) + 43.30$).
2. **Heat and Momentum Transport:**
   - Dixon and Li & Finlayson packed-bed wall Nusselt correlations yielding realistic bed-to-wall heat transfer coefficients ($U \approx 200\text{–}350\text{ W/(m}^2\cdot\text{K)}$).
   - Ergun equation for packed-bed pressure drops across tube and shell.
   - Boiling-water film coefficient ($h_{\text{cool}} = 6000\text{ W/(m}^2\cdot\text{K)}$) for conventional reactor cooling jacket.
3. **Rigorous Conservation Closures:**
   - Elemental mass balance closures (C, H, O, N, Ar) $\le 10^{-15}$.
   - Independent stream enthalpy balance closures $\le 0.005\%$.
   - Zero state-clamp activations on accepted solution trajectories.
4. **Counter-Current Picard Solver:**
   - Integrates both tube and shell streams strictly forward along their physical flow paths.
   - Employs under-relaxation ($\alpha = 0.70$) with continuous interpolation to achieve dual convergence on temperature ($\Delta T < 10^{-3}\text{ K}$) and molar flows ($\Delta F_i/F_i < 10^{-3}$).
5. **Fair Multi-Objective Optimization:**
   - Baseline CMR optimized over $(T_{\text{exo},0}, T_{\text{coolant}})$ to establish a rigorous, fair comparison.
   - NSGA-II Pareto optimization mapping trade-offs between methanol yield, benzene/hydrogen coproduction, and catalyst thermal protection ($T_{\text{max}} \le 535\text{ K}$).

---

## 4. Key Results Summary

| Performance Metric | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |
|---|---|---|---|---|---|
| **Methanol Production (t/d)** | 209.62 | **226.87** | 211.04 | 207.56 | **164.08** |
| **Methanol Yield (%)** | 28.53% | **30.88%** | 28.72% | 28.25% | **25.50%** |
| **Benzene Coproduct (t/d)** | 0.00 | 0.00 | **100.71** | **112.37** | **66.56** |
| **Hydrogen Coproduct (t/d)** | 0.00 | 0.00 | **7.80** | **8.70** | **5.15** |
| **Peak Bed Temp (K)** | 524.4 | 535.0 | 539.6 | 550.5 | **525.7** |
| **Shell Exit Pressure (bar)** | — | — | 1.02 | 0.95 | **1.13** |
| **Reaction Zone Cooling Utility**| 5.57 MW | 5.64 MW | **0.0 kW** | **0.0 kW** | **0.0 kW** |
| **Feed Preheating Duty** | 10.51 MW | 10.51 MW | 17.17 MW | 17.17 MW | 17.17 MW |

---

## 5. Audit & Compliance

All phases (Phase 0 through Phase 6) have been completed and audited. Detailed technical logs are available in:
- `FINAL_AUDIT.md`: Mapping of defects D1–D10 before and after refactoring.
- `REVIEW_LOG.md`: Step-by-step logs with execution evidence and bilingual summaries.
- `OPEN_ISSUES.md`: Detailed resolution history of physical modeling issues.

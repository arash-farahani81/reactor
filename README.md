# A Critical Re-Assessment and Thermodynamic Ceiling Analysis of Thermally Coupled Methanol–Cyclohexane Reactors

A rigorous, single-source-of-truth (SSOT) multiphysics simulation, forensic benchmarking, and multi-objective optimization framework evaluating a **Thermally Coupled Shell-and-Tube Reactor (TCR)** that pairs:

- **Tube Side (Exothermic):** Pressurized methanol synthesis ($\text{CO}/\text{CO}_2$ hydrogenation + RWGS over commercial $\text{Cu/ZnO/Al}_2\text{O}_3$, $P_0 = 76.98\text{ bar}$, $L = 7.022\text{ m}$, $N_{\text{tubes}} = 2962$, locked equal-throughput syngas feed $F_{\text{exo},0} = 0.640\text{ mol/s}$ per tube, baseline $\eta_{\text{exo}} = 0.35$).
- **Shell Side (Endothermic):** Cyclohexane dehydrogenation ($\text{C}_6\text{H}_{12} \rightleftharpoons \text{C}_6\text{H}_6 + 3\text{H}_2$ over $\text{Pt/Al}_2\text{O}_3$, $5\text{ mol}\%$ $\text{H}_2$ production co-feed, pressurized shell inlet to guarantee positive discharge pressure).

All numerical outputs are generated dynamically from `config.yaml` via `run_all.py` into `simulation_results.json` (SSOT), then compiled automatically into `TCR_Methanol_Cyclohexane_Manuscript.md` and `simulation_summary.md` with zero hardcoded simulation results and strict automated consistency & LaTeX verification.

---

## 1. Five Core Scientific Contributions (Path A)

1. **Khademi et al. (2009) Forensic Benchmark Reconciliation:**
   - Reproduces the literature benchmark using exact Table 7/8 shell geometry ($D_{\text{shell}} = 0.060\text{ m}$, $A_{\text{shell}} = 1.44199 \times 10^{-3}\text{ m}^2$, $\varepsilon_{\text{shell}} = 0.39$) and authentic 3-reaction Graaf et al. (1988) kinetics via exact 3-to-2 stoichiometric projection ($r_{\text{meoh}} = r_1 + r_3$, $r_{\text{rwgs}} = r_2 - r_1$).
   - Documents two severe typographical errors in Khademi Table 1 (CO hydrogenation $E_a$ misprinted as $63{,}000$ vs. $113{,}000\text{ J/mol}$; RWGS pre-exponential misprinted as $9.64 \times 10^6$ vs. $9.64 \times 10^{11}$) and the Table 4 cyclohexane equilibrium typo ($B = +3190\text{ K}$ vs. exact NIST-anchored Kirchhoff van 't Hoff thermodynamics, $-\Delta H_r^\circ / R \approx -24{,}704\text{ K}$).
   - Reports both **gross Eq. (28) methanol yield** (which includes the $+3.57\text{ percentage-point}$ recycle methanol offset) and **net $\text{CO}_x$ methanol yield**, explaining the cold-shell-inlet ($423.0\text{ K}$) sensible-heating and tube over-cooling mechanism.
2. **Optimal Temperature Profile (OTP) Thermodynamic Yield Ceiling:**
   - Computes the theoretical upper bound on single-pass methanol yield via a 5-node L-BFGS-B Optimal Temperature Profile (bounded by $T_{\max} = 535.0\text{ K}$) alongside the best isothermal tube benchmark at both $\eta_{\text{exo}} = 0.35$ and $\eta_{\text{exo}} = 1.00$.
   - Demonstrates that a fairly optimized conventional boiling-water methanol reactor (CMR) at identical throughput already operates close to the theoretical OTP ceiling, refuting yield-driven claims of TCR superiority.
3. **First-Law Battery-Limit Energy & Steam-Credit Accounting:**
   - Strictly separates First-Law thermodynamic energy qualities—**external heating & steam-credit duties ($\text{MW}_{\text{th}}$)**, **cooling-water rejection ($\text{MW}_{\text{cw}}$)**, and **electrical shaft work ($\text{MW}_e$)**—and never sums across distinct energy grades.
   - Credits the conventional boiling-water reactor with its exportable medium-pressure steam generation while debiting the coupled reactor for cyclohexane vaporization/preheating, extractive distillation of the benzene–cyclohexane azeotrope ($\Delta T_{\text{bp}} = 0.6\text{ K}$), flash condensation, and hydrogen PSA purification work.
4. **Unified Dynamic Bauer–Schlünder Conduction, Radial Overshoot Warning & Shell Hydraulics:**
   - Replaces static radial conductivity approximations with unified dynamic Bauer–Schlünder $k_{er}(z)$ coupled to Wassiljewa–Mason–Saxena mixture gas conductivity and Dixon / Li & Finlayson wall heat transfer correlations.
   - Quantifies radial centerline temperature overshoot ($\Delta T_{\text{radial}}$) under 1D bulk temperature constraints ($T_{\text{bulk,max}} \le 535.0\text{ K}$), demonstrating why industrial 1D design optimization must incorporate explicit radial centerline safety derating to prevent copper catalyst sintering.
   - Audits upward-flow shell pressure drop ($\Delta P_{\text{shell}}$ via Ergun + Eisfeld–Schnitzlein wall corrections) against static catalyst bed weight ($P_{\text{bed}} = \rho_{b,\text{endo}} g L$) and flags mechanical hold-down grid requirements whenever $\Delta P_{\text{shell}} / P_{\text{bed}} \ge 0.85$.
5. **High-Pressure Real-Gas Peng–Robinson EOS Audit:**
   - Compares ideal-gas partial pressures ($\phi_i = 1$) against full Peng–Robinson EOS fugacities ($f_i = \phi_i y_i P$) at $76.98\text{ bar}$ evaluated at the exact optimized CMR operating point, capturing how non-ideal fugacities shift RWGS water formation and intensify kinetic water inhibition on $\text{Cu/ZnO/Al}_2\text{O}_3$.

---

## 2. Quick Start & Reproduction Commands

Ensure Python 3.10+ is installed, then install scientific dependencies and run the full SSOT pipeline:

```bash
# 1. Install Python scientific dependencies
pip install -r tcr/requirements.txt

# 2. Run the full numerical pipeline (CMR, TCR co/counter-current, NSGA-II,
#    Khademi benchmark, OTP ceiling, PR-EOS audit, sensitivity sweeps, and figures)
python3 run_all.py

# 3. Compile the zero-hardcode manuscript from simulation_results.json (SSOT)
python3 build_manuscript.py

# 4. Run automated numerical, semantic, and LaTeX consistency verification
python3 check_manuscript_consistency.py

# 5. Build the unified review & replication dossier (optional)
python3 build_ai_dossier.py

# 6. Execute the full automated verification & physics test suite (T1-T15)
pytest tests/ -q
```

---

## 3. Repository Architecture (SSOT)

```text
.
├── config.yaml                            # Single source of truth physical & solver parameters
├── run_all.py                             # Master pipeline entry point -> simulation_results.json
├── build_manuscript.py                    # Zero-hardcode template compiler + LaTeX syntax gate
├── check_manuscript_consistency.py        # Automated SSOT numerical & semantic consistency checker
├── build_ai_dossier.py                    # Review & replication dossier generator
├── manuscript_template.md                 # Dynamic Path A manuscript template ({{ ... }} tags)
├── TCR_Methanol_Cyclohexane_Manuscript.md # Compiled manuscript (SSOT output)
├── simulation_results.json                # Master SSOT numerical database
├── simulation_summary.md                  # Auto-generated executive numerical summary (SSOT)
├── tcr/                                   # Core modular multiphysics package
│   ├── requirements.txt                   # Python scientific dependency manifest
│   ├── config.py                          # Typed configuration dataclasses & YAML loader
│   ├── thermo.py                          # Cp(T), Hi(T), NIST-anchored Kirchhoff Kp(T), PR-EOS
│   ├── geometry.py                        # Tube and annular unit-cell flow geometry
│   ├── transport.py                       # Wassiljewa k_g, dynamic Bauer-Schluender k_er, Dixon/Ergun
│   ├── kinetics.py                        # Production VBF (1996), authentic Graaf (1988), Itoh (1987)
│   ├── solvers.py                         # Stiff Radau IIA IVP, Picard BVP & solve_bvp cross-check
│   ├── metrics.py                         # Mass/enthalpy closures, steam credit, bed-lifting & k_er
│   ├── optimize.py                        # Equal-throughput CMR optimization & TCR NSGA-II
│   ├── benchmark.py                       # Khademi (2009) forensic audit, OTP ceiling & PR-EOS audit
│   ├── plots.py                           # Publication-quality figures (Figures 1-5)
│   └── io.py                              # SSOT JSON serialization, summary & battery-limit duties
├── tests/                                 # Automated physics & SSOT verification suite (T1-T15)
│   ├── test_physics.py                    # Conservation, Kirchhoff consistency, PR-EOS, OTP, hydraulics
│   └── test_full_suite.py                 # Grid convergence, Picard/BVP agreement, SSOT no-hardcode gate
├── figures/                               # High-resolution generated PNG figures
├── baseline/                              # Archived pre-refactor benchmark outputs
└── legacy/                                # Quarantined pre-modular scripts
```

---

## 4. Single Source of Truth for Numerical Results

To prevent stale documentation drift, **this `README.md` intentionally contains no hardcoded simulation outputs**. After running `python3 run_all.py` and `python3 build_manuscript.py`, consult the auto-generated SSOT artifacts for exact, verified numerical results:

- **`simulation_summary.md`** — Concise executive tables covering all five equal-throughput reactor regimes, OTP ceiling, PR-EOS audit, and First-Law battery-limit energy accounting.
- **`TCR_Methanol_Cyclohexane_Manuscript.md`** — Complete publication-ready Path A manuscript with all dynamically injected tables, equations, and embedded figure references.
- **`simulation_results.json`** — Machine-readable SSOT dataset backing every reported number and checked by `check_manuscript_consistency.py`.

---

## 5. Verification & Test Suite Coverage (`T1–T15`)

Running `pytest tests/ -q` validates:
- **Conservation & Thermodynamics (`test_physics.py`):** C/H/O/N/Ar elemental closure ($\le 10^{-12}$), independent stream enthalpy balance ($\le 0.5\%$), exact Kirchhoff van 't Hoff consistency ($R T^2 d\ln K_p / dT \equiv \Delta H_r^\circ(T)$), Graaf 3-to-2 stoichiometric projection & Hess's Law equivalence, Peng–Robinson EOS fugacity audit, dynamic Bauer–Schlünder $k_{er}$ & centerline temperature calculation, CMR steam credit, and upward-flow shell bed-lifting ratio.
- **Solvers, Optimization & SSOT Integrity (`test_full_suite.py`):** Axial grid convergence, counter-current Picard vs. `scipy.integrate.solve_bvp` cross-check, equal-throughput $F_{\text{exo},0} = 0.640\text{ mol/s}$ lock, OTP ceiling monotonicity, zero-hardcode template rendering, and strict manuscript consistency checks.

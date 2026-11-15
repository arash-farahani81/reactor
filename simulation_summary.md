# Thermally Coupled Reactor (TCR) Simulation Performance Summary

Generated from Single Source of Truth architecture (`tcr/` package).

| Performance Indicator | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Methanol Single-Pass Yield (%)** | 28.55% | 30.89% | 28.36% | 27.76% | **29.17%** |
| **Methanol Production (mol/s)** | 75.78 | 81.99 | 75.25 | 73.68 | **77.42** |
| **Methanol Production (ton/day)** | 209.79 | 226.98 | 208.33 | 203.97 | **214.33** |
| **CO₂ Conversion (%)** | 19.08% | 21.43% | 22.91% | 21.93% | **21.59%** |
| **Total COx Conversion (%)** | 28.55% | 30.89% | 28.36% | 27.76% | **29.17%** |
| **Cyclohexane Conversion (%)** | N/A | N/A | 11.06% | 12.22% | **7.89%** |
| **Benzene Production (mol/s)** | 0.00 | 0.00 | 13.11 | 14.48 | **15.64** |
| **Benzene Production (ton/day)** | 0.00 | 0.00 | 88.45 | 97.75 | **105.56** |
| **Coproduced H₂ Rate (mol/s)** | 0.00 | 0.00 | 39.32 | 43.45 | **46.93** |
| **Coproduced H₂ Rate (ton/day)** | 0.00 | 0.00 | 6.85 | 7.57 | **8.17** |
| **Peak Bed Temperature (K)** | 524.37 | 535.00 | 540.88 | 552.50 | **527.51** |
| **Tube Pressure Drop (bar)** | 0.45 | 0.45 | 0.46 | 0.47 | **0.45** |
| **Shell Pressure Drop (bar)** | N/A | N/A | 0.68 | 0.72 | **1.79** |
| **Reaction Cooling Utility Required** | Boiling Water (5.58 MW) | Boiling Water (5.50 MW) | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** |

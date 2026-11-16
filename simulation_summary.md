# Thermally Coupled Reactor (TCR) Simulation Performance Summary

Generated from Single Source of Truth architecture (`tcr/` package).

| Performance Indicator | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Methanol Single-Pass Yield (%)** | 33.44% | 36.16% | 32.00% | 28.74% | **33.62%** |
| **Methanol Production (mol/s)** | 88.76 | 95.97 | 84.94 | 76.29 | **89.21** |
| **Methanol Production (ton/day)** | 245.72 | 265.68 | 235.14 | 211.19 | **246.98** |
| **CO₂ Conversion (%)** | 21.32% | 23.67% | 25.86% | 24.32% | **25.52%** |
| **Total COx Conversion (%)** | 33.44% | 36.16% | 32.00% | 28.74% | **33.62%** |
| **Cyclohexane Conversion (%)** | N/A | N/A | 13.14% | 12.02% | **4.59%** |
| **Benzene Production (mol/s)** | 0.00 | 0.00 | 14.79 | 13.53 | **9.20** |
| **Benzene Production (ton/day)** | 0.00 | 0.00 | 99.80 | 91.28 | **62.12** |
| **Coproduced H₂ Rate (mol/s)** | 0.00 | 0.00 | 44.37 | 40.58 | **27.61** |
| **Coproduced H₂ Rate (ton/day)** | 0.00 | 0.00 | 7.73 | 7.07 | **4.81** |
| **Peak Bed Temperature (K)** | 514.27 | 519.42 | 545.73 | 544.53 | **531.00** |
| **Tube Pressure Drop (bar)** | 0.45 | 0.45 | 0.47 | 0.48 | **0.46** |
| **Shell Pressure Drop (bar)** | N/A | N/A | 0.69 | 0.70 | **1.84** |
| **Shell Exit Pressure (bar)** | N/A | N/A | 1.81 | 1.80 | **1.16** |
| **Syngas Feed Rate (mol/s)** | 0.640 | 0.640 | 0.640 | 0.640 | **0.640** |
| **Reaction Cooling Utility Required** | Boiling Water (7.04 MW) | Boiling Water (6.50 MW) | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** |

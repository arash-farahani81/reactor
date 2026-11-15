# Thermally Coupled Reactor (TCR) Simulation Performance Summary

Generated from Single Source of Truth architecture (`tcr/` package).

| Performance Indicator | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Methanol Single-Pass Yield (%)** | 28.53% | 30.88% | 28.72% | 28.25% | **25.50%** |
| **Methanol Production (mol/s)** | 75.72 | 81.95 | 76.23 | 74.97 | **59.27** |
| **Methanol Production (ton/day)** | 209.62 | 226.87 | 211.04 | 207.56 | **164.08** |
| **CO₂ Conversion (%)** | 19.06% | 21.34% | 22.74% | 21.87% | **19.87%** |
| **Total COx Conversion (%)** | 28.53% | 30.88% | 28.72% | 28.25% | **25.50%** |
| **Cyclohexane Conversion (%)** | N/A | N/A | 12.60% | 14.05% | **8.25%** |
| **Benzene Production (mol/s)** | 0.00 | 0.00 | 14.92 | 16.65 | **9.86** |
| **Benzene Production (ton/day)** | 0.00 | 0.00 | 100.71 | 112.37 | **66.56** |
| **Coproduced H₂ Rate (mol/s)** | 0.00 | 0.00 | 44.77 | 49.95 | **29.59** |
| **Coproduced H₂ Rate (ton/day)** | 0.00 | 0.00 | 7.80 | 8.70 | **5.15** |
| **Peak Bed Temperature (K)** | 524.36 | 535.00 | 539.59 | 550.48 | **525.74** |
| **Tube Pressure Drop (bar)** | 0.45 | 0.45 | 0.46 | 0.47 | **0.34** |
| **Shell Pressure Drop (bar)** | N/A | N/A | 0.98 | 1.05 | **0.87** |
| **Reaction Cooling Utility Required** | Boiling Water (5.57 MW) | Boiling Water (5.64 MW) | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** |

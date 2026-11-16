# Thermally Coupled Reactor (TCR) Simulation Performance Summary

Generated from Single Source of Truth architecture (`tcr/` package).

| Performance Indicator | CMR Baseline | CMR Optimized | TCR Co-Current | TCR Counter-Current | TCR NSGA-II Optimum |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Methanol Single-Pass Yield (%)** | 20.94% | 30.05% | 29.07% | 29.10% | **29.93%** |
| **Methanol Production (mol/s)** | 55.58 | 79.75 | 77.14 | 77.23 | **79.44** |
| **Methanol Production (ton/day)** | 153.88 | 220.77 | 213.55 | 213.81 | **219.92** |
| **CO₂ Conversion (%)** | 15.19% | 20.96% | 22.70% | 21.11% | **21.94%** |
| **Total COx Conversion (%)** | 20.94% | 30.05% | 29.07% | 29.10% | **29.93%** |
| **Cyclohexane Conversion (%)** | N/A | N/A | 12.55% | 17.54% | **10.15%** |
| **Benzene Production (mol/s)** | 0.00 | 0.00 | 14.13 | 19.74 | **18.88** |
| **Benzene Production (ton/day)** | 0.00 | 0.00 | 95.35 | 133.23 | **127.43** |
| **Coproduced H₂ Rate (mol/s)** | 0.00 | 0.00 | 42.38 | 59.23 | **56.64** |
| **Coproduced H₂ Rate (ton/day)** | 0.00 | 0.00 | 7.38 | 10.32 | **9.87** |
| **Peak Bed Temperature (K)** | 507.34 | 525.23 | 535.01 | 554.05 | **529.47** |
| **Tube Pressure Drop (bar)** | 0.44 | 0.45 | 0.46 | 0.47 | **0.45** |
| **Shell Pressure Drop (bar)** | N/A | N/A | 0.66 | 0.75 | **1.75** |
| **Shell Exit Pressure (bar)** | N/A | N/A | 1.84 | 1.75 | **1.21** |
| **Syngas Feed Rate (mol/s)** | 0.640 | 0.640 | 0.640 | 0.640 | **0.640** |
| **Reaction Cooling Utility Required** | Boiling Water (4.23 MW) | Boiling Water (5.12 MW) | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** | **0.0 kW (Autothermal)** |

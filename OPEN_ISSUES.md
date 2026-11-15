# OPEN ISSUES — Thermally Coupled Reactor (TCR) Project

Tracked scientific, physical, and implementation open questions and benchmark discrepancies.

| Issue ID | Category | Description | Status | Resolution Path |
|---|---|---|---|---|
| **ISSUE-01** | Benchmark Kinetics | Khademi et al. (2009) Table 4 typo: printed $K_p = 4.89 \times 10^{35} \exp(-3190/T)\text{ Pa}^3$ has apparent activation temperature $-3190\text{ K}$, whereas $\Delta H^\circ / R \approx 205400 / 8.314 = 24704\text{ K}$. | Open | Implement both `khademi_literal` and `khademi_thermo` in benchmark mode; document the exact mathematical discrepancy. |
| **ISSUE-02** | Wall Heat Transfer | The Chilton-Colburn packed bed correlation in Khademi Table 5 is for fluid-particle mass/heat transfer ($j_D, j_H$), not bed-to-wall heat transfer ($h_w$). | Open | Implement Dixon / Li & Finlayson packed-bed-to-wall correlation for gas beds and high boiling water film coefficient ($h_{\text{cool}} = 6000\text{ W/(m}^2\cdot\text{K)}$) for CMR. |
| **ISSUE-03** | Thermodynamic Properties | Constant $C_p$ and constant $\Delta H_{r,298}^\circ$ cause stream enthalpy balance errors of $13\text{–}16\%$. | Open | Implement NASA-7 / NIST polynomial $C_p(T)$ and temperature-dependent $\Delta H_r(T) = \Delta H_{r,298}^\circ + \int_{298.15}^T \Delta C_p(T')\,dT'$. |
| **ISSUE-04** | Shell Geometry | The $6.5\times$ area multiplier for shell flow area lacked rigorous geometric derivation. | Open | Define shell flow area and hydraulic diameter strictly from tube pitch bundle geometry ($p_{\text{pitch}} = 1.33 D_o$). |
| **ISSUE-05** | Fair Optimization | Unoptimized conventional reactor vs. NSGA-II optimized TCR. | Open | Formulate and solve CMR single-objective optimization over $(T_{\text{exo},0}, T_{\text{coolant}})$ to provide fair benchmark. |

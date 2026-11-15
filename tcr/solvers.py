"""
Numerical solvers for:
1. Conventional Jacketed Reactor (CMR)
2. Co-Current Thermally Coupled Reactor (Radau IIA IVP)
3. Counter-Current Thermally Coupled Reactor (Under-relaxed Picard BVP + solve_bvp cross-check)
"""

import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from .thermo import (
    R_GAS, calc_cp_species, calc_dH_rxn_meoh, calc_dH_rxn_rwgs, calc_dH_rxn_ch, MW
)
from .geometry import ReactorGeometry
from .config import CatalystBedConfig, TransportConfig, SolverConfig
from .kinetics import (
    kinetics_vbf_production, kinetics_ch_production,
    kinetics_graaf_benchmark, kinetics_itoh_benchmark
)
from .transport import calc_U_tcr, calc_U_cmr, calc_ergun_dp

SPECIES_TUBE = ['H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2']
SPECIES_SHELL = ['C6H12', 'C6H6', 'H2', 'Ar']

# Terminal event guards (Rule R5)
def make_event_pressure_collapse(P_min):
    def event(z, S):
        return S[-2] - P_min
    event.terminal = True
    event.direction = -1
    return event

def make_event_negative_flow():
    def event(z, S):
        # S[0:7] are tube flows
        return float(np.min(S[0:7]) + 1e-6)
    event.terminal = True
    event.direction = -1
    return event


def solve_conventional(geom: ReactorGeometry, cat: CatalystBedConfig, trans: TransportConfig,
                       sol_cfg: SolverConfig, y0: np.ndarray, T0: float, P0: float,
                       F_total0: float, T_coolant: float = 503.15,
                       kinetics_mode: str = "production_vbf") -> object:
    """
    Solves conventional methanol reactor cooled by boiling water jacket.
    State vector: [F_exo (7), T_exo, P_exo]
    """
    F_init = y0 * F_total0
    S0 = np.concatenate([F_init, [T0, P0]])
    z_eval = np.linspace(0.0, geom.L, sol_cfg.z_points)

    def odes_conv(z, S):
        # Numerical guard ONLY inside domain bounds
        F_exo = np.maximum(S[0:7], 1e-15)
        T_exo = max(S[7], 200.0)
        P_exo = max(S[8], 1e4)

        F_tot = np.sum(F_exo)
        y_exo = F_exo / F_tot

        if kinetics_mode == "production_vbf":
            r_meoh, r_rwgs = kinetics_vbf_production(T_exo, P_exo, y_exo, cat.eta)
            dH_meoh = calc_dH_rxn_meoh(T_exo)
            dH_rwgs = calc_dH_rxn_rwgs(T_exo)
        else:
            r1, r2, r3 = kinetics_graaf_benchmark(T_exo, P_exo, y_exo, cat.eta)
            r_meoh = r1 + r3
            r_rwgs = r2
            dH_meoh = calc_dH_rxn_meoh(T_exo)
            dH_rwgs = calc_dH_rxn_rwgs(T_exo)

        dF = np.zeros(7)
        dF[0] = geom.A_tube * cat.rho_b * (-3.0 * r_meoh - r_rwgs)
        dF[1] = geom.A_tube * cat.rho_b * (-r_meoh - r_rwgs)
        dF[2] = geom.A_tube * cat.rho_b * (r_rwgs)
        dF[3] = geom.A_tube * cat.rho_b * (r_meoh)
        dF[4] = geom.A_tube * cat.rho_b * (r_meoh + r_rwgs)
        dF[5] = 0.0
        dF[6] = 0.0

        # Exact temperature-dependent heat capacity flow (W/K)
        FCp = sum(F_exo[i] * calc_cp_species(SPECIES_TUBE[i], T_exo) for i in range(7))

        # Dynamic overall heat transfer coefficient with boiling water jacket
        U_local = calc_U_cmr(T_exo, P_exo, F_exo, SPECIES_TUBE, geom, cat, trans)
        Q_wall = U_local * geom.a_H * (T_exo - T_coolant)

        Q_rxn = geom.A_tube * cat.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
        dT = (Q_rxn - Q_wall) / max(FCp, 1e-3)

        # Ergun pressure drop
        mw_avg = sum(y_exo[i] * MW[SPECIES_TUBE[i]] for i in range(7))
        rho_g = (P_exo * mw_avg) / (R_GAS * T_exo)
        u_s = (F_tot * R_GAS * T_exo) / (P_exo * geom.A_tube)
        dP = calc_ergun_dp(u_s, rho_g, trans.mu_exo, cat.dp, cat.eps)

        return np.concatenate([dF, [dT, dP]])

    sol = solve_ivp(
        odes_conv, (0.0, geom.L), S0, t_eval=z_eval, method='Radau',
        rtol=sol_cfg.rtol, atol=sol_cfg.atol
    )

    # Post-calculate exact wall heat flux along z
    Q_wall_arr = np.zeros(len(sol.t))
    for i in range(len(sol.t)):
        T_e = sol.y[7, i]
        P_e = sol.y[8, i]
        F_e = sol.y[0:7, i]
        U_loc = calc_U_cmr(T_e, P_e, F_e, SPECIES_TUBE, geom, cat, trans)
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_coolant)
    sol.Q_wall_arr = Q_wall_arr

    # Check clamp count on accepted trajectory
    sol.clamp_violations = int(np.sum(sol.y[0:7, :] < -1e-6) + np.sum(sol.y[8, :] < 0.2 * P0))
    return sol


def solve_cocurrent(geom: ReactorGeometry, cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
                    trans: TransportConfig, sol_cfg: SolverConfig,
                    y_exo_0: np.ndarray, T_exo_0: float, P_exo_0: float, F_exo_total_0: float,
                    y_endo_0: np.ndarray, T_endo_0: float, P_endo_0: float, F_endo_total_0: float,
                    kinetics_mode: str = "production_vbf") -> object:
    """
    Solves Co-Current TCR (IVP from z=0 to z=L) using Radau IIA.
    State vector: [F_exo (7), F_endo (len_endo), T_exo, T_endo, P_exo, P_endo]
    """
    len_endo = len(y_endo_0)
    species_shell = SPECIES_SHELL[:len_endo]

    F_exo_init = y_exo_0 * F_exo_total_0
    F_endo_init = y_endo_0 * F_endo_total_0
    S0 = np.concatenate([F_exo_init, F_endo_init, [T_exo_0, T_endo_0, P_exo_0, P_endo_0]])
    z_eval = np.linspace(0.0, geom.L, sol_cfg.z_points)

    def odes_co(z, S):
        F_exo = np.maximum(S[0:7], 1e-15)
        F_endo = np.maximum(S[7:7+len_endo], 1e-15)
        T_exo = max(S[7+len_endo], 200.0)
        T_endo = max(S[8+len_endo], 200.0)
        P_exo = max(S[9+len_endo], 1e4)
        P_endo = max(S[10+len_endo], 1e4)

        F_exo_tot = np.sum(F_exo)
        F_endo_tot = np.sum(F_endo)
        y_exo = F_exo / F_exo_tot
        y_endo = F_endo / F_endo_tot

        if kinetics_mode == "production_vbf":
            r_meoh, r_rwgs = kinetics_vbf_production(T_exo, P_exo, y_exo, cat_exo.eta)
            r_ch = kinetics_ch_production(T_endo, P_endo, y_endo, cat_endo.eta)
        else:
            r1, r2, r3 = kinetics_graaf_benchmark(T_exo, P_exo, y_exo, cat_exo.eta)
            r_meoh = r1 + r3
            r_rwgs = r2
            r_ch = kinetics_itoh_benchmark(T_endo, P_endo, y_endo, mode="khademi_thermo", eta_endo=cat_endo.eta)

        dH_meoh = calc_dH_rxn_meoh(T_exo)
        dH_rwgs = calc_dH_rxn_rwgs(T_exo)
        dH_ch = calc_dH_rxn_ch(T_endo)

        dF_exo = np.zeros(7)
        dF_exo[0] = geom.A_tube * cat_exo.rho_b * (-3.0 * r_meoh - r_rwgs)
        dF_exo[1] = geom.A_tube * cat_exo.rho_b * (-r_meoh - r_rwgs)
        dF_exo[2] = geom.A_tube * cat_exo.rho_b * (r_rwgs)
        dF_exo[3] = geom.A_tube * cat_exo.rho_b * (r_meoh)
        dF_exo[4] = geom.A_tube * cat_exo.rho_b * (r_meoh + r_rwgs)

        dF_endo = np.zeros(len_endo)
        dF_endo[0] = geom.A_shell * cat_endo.rho_b * (-r_ch)
        dF_endo[1] = geom.A_shell * cat_endo.rho_b * (r_ch)
        dF_endo[2] = geom.A_shell * cat_endo.rho_b * (3.0 * r_ch)

        # Dynamic overall U evaluated from local states of both streams
        U_local = calc_U_tcr(T_exo, P_exo, F_exo, SPECIES_TUBE,
                             T_endo, P_endo, F_endo, species_shell,
                             geom, cat_exo, cat_endo, trans)
        Q_wall = U_local * geom.a_H * (T_exo - T_endo)

        FCp_exo = sum(F_exo[i] * calc_cp_species(SPECIES_TUBE[i], T_exo) for i in range(7))
        FCp_endo = sum(F_endo[i] * calc_cp_species(species_shell[i], T_endo) for i in range(len_endo))

        Q_rxn_exo = geom.A_tube * cat_exo.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
        dT_exo = (Q_rxn_exo - Q_wall) / max(FCp_exo, 1e-3)

        Q_rxn_endo = geom.A_shell * cat_endo.rho_b * (-dH_ch) * r_ch
        dT_endo = (Q_rxn_endo + Q_wall) / max(FCp_endo, 1e-3)

        # Ergun pressure drops
        mw_exo = sum(y_exo[i] * MW[SPECIES_TUBE[i]] for i in range(7))
        rho_exo = (P_exo * mw_exo) / (R_GAS * T_exo)
        u_exo = (F_exo_tot * R_GAS * T_exo) / (P_exo * geom.A_tube)
        dP_exo = calc_ergun_dp(u_exo, rho_exo, trans.mu_exo, cat_exo.dp, cat_exo.eps)

        mw_endo = sum(y_endo[i] * MW[species_shell[i]] for i in range(len_endo))
        rho_endo = (P_endo * mw_endo) / (R_GAS * T_endo)
        u_endo = (F_endo_tot * R_GAS * T_endo) / (P_endo * geom.A_shell)
        dP_endo = calc_ergun_dp(u_endo, rho_endo, trans.mu_endo, cat_endo.dp, cat_endo.eps)

        return np.concatenate([dF_exo, dF_endo, [dT_exo, dT_endo, dP_exo, dP_endo]])

    sol = solve_ivp(
        odes_co, (0.0, geom.L), S0, t_eval=z_eval, method='Radau',
        rtol=sol_cfg.rtol, atol=sol_cfg.atol
    )

    # Post-calculate exact wall heat flux along z
    n_pts = len(sol.t)
    Q_wall_arr = np.zeros(n_pts)
    for i in range(n_pts):
        T_e = sol.y[7+len_endo, i]
        T_n = sol.y[8+len_endo, i]
        P_e = sol.y[9+len_endo, i]
        P_n = sol.y[10+len_endo, i]
        F_e = sol.y[0:7, i]
        F_n = sol.y[7:7+len_endo, i]
        U_loc = calc_U_tcr(T_e, P_e, F_e, SPECIES_TUBE,
                           T_n, P_n, F_n, species_shell,
                           geom, cat_exo, cat_endo, trans)
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_n)
    sol.Q_wall_arr = Q_wall_arr

    # Check clamp count on accepted trajectory
    sol.clamp_violations = int(
        np.sum(sol.y[0:7, :] < -1e-6)
        + np.sum(sol.y[7:7+len_endo, :] < -1e-6)
        + np.sum(sol.y[9+len_endo, :] < 0.2 * P_exo_0)
        + np.sum(sol.y[10+len_endo, :] < 0.2 * P_endo_0)
    )
    return sol


def solve_countercurrent(geom: ReactorGeometry, cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
                         trans: TransportConfig, sol_cfg: SolverConfig,
                         y_exo_0: np.ndarray, T_exo_0: float, P_exo_0: float, F_exo_total_0: float,
                         y_endo_0: np.ndarray, T_endo_0: float, P_endo_0: float, F_endo_total_0: float,
                         kinetics_mode: str = "production_vbf", cross_check_bvp: bool = False) -> tuple[object, dict]:
    """
    Solves Counter-Current TCR (Two-Point BVP) via Under-Relaxed Picard Fixed-Point Iteration.
    - Tube integrates forward from z = 0 to L given interpolated shell fields.
    - Shell integrates along physical flow path from z = L to 0 given interpolated tube fields.
    - Local U(z) evaluated from local interpolated states of both streams (Q_wall_tube == Q_wall_shell).
    - Checks dual convergence on both temperature and molar flows.
    """
    z_points = sol_cfg.z_points
    z_grid = np.linspace(0.0, geom.L, z_points)
    len_endo = len(y_endo_0)
    species_shell = SPECIES_SHELL[:len_endo]

    # Initial guesses: linear temperature profile from T_exo_0 to T_endo_0
    T_shell_prof = np.linspace(T_exo_0 + 5.0, T_endo_0, z_points)
    P_shell_prof = np.full(z_points, P_endo_0)
    F_shell_prof = np.tile(y_endo_0 * F_endo_total_0, (z_points, 1)).T

    F_exo_init = y_exo_0 * F_exo_total_0
    F_endo_init = y_endo_0 * F_endo_total_0

    sol_tube = None
    F_shell_new = None
    T_shell_new = None
    P_shell_new = None
    converged = False
    iterations_run = 0
    final_res_T = 0.0
    final_res_F = 0.0

    for it in range(sol_cfg.picard_max_iter):
        iterations_run = it + 1

        # ── Step 1: Integrate Tube Forward (z = 0 -> L) ──
        def ode_tube(z, S_t):
            F_exo = np.maximum(S_t[0:7], 1e-15)
            T_exo = max(S_t[7], 200.0)
            P_exo = max(S_t[8], 1e4)

            F_tot = np.sum(F_exo)
            y_exo = F_exo / F_tot

            # Interpolate local shell state at coordinate z
            T_s = float(np.interp(z, z_grid, T_shell_prof))
            P_s = float(np.interp(z, z_grid, P_shell_prof))
            F_s = np.array([float(np.interp(z, z_grid, F_shell_prof[k, :])) for k in range(len_endo)])

            if kinetics_mode == "production_vbf":
                r_meoh, r_rwgs = kinetics_vbf_production(T_exo, P_exo, y_exo, cat_exo.eta)
            else:
                r1, r2, r3 = kinetics_graaf_benchmark(T_exo, P_exo, y_exo, cat_exo.eta)
                r_meoh = r1 + r3
                r_rwgs = r2

            dH_meoh = calc_dH_rxn_meoh(T_exo)
            dH_rwgs = calc_dH_rxn_rwgs(T_exo)

            dF = np.zeros(7)
            dF[0] = geom.A_tube * cat_exo.rho_b * (-3.0 * r_meoh - r_rwgs)
            dF[1] = geom.A_tube * cat_exo.rho_b * (-r_meoh - r_rwgs)
            dF[2] = geom.A_tube * cat_exo.rho_b * (r_rwgs)
            dF[3] = geom.A_tube * cat_exo.rho_b * (r_meoh)
            dF[4] = geom.A_tube * cat_exo.rho_b * (r_meoh + r_rwgs)

            # Local U(z) from local states
            U_local = calc_U_tcr(T_exo, P_exo, F_exo, SPECIES_TUBE,
                                 T_s, P_s, F_s, species_shell,
                                 geom, cat_exo, cat_endo, trans)
            Q_wall = U_local * geom.a_H * (T_exo - T_s)

            FCp = sum(F_exo[i] * calc_cp_species(SPECIES_TUBE[i], T_exo) for i in range(7))
            Q_rxn = geom.A_tube * cat_exo.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
            dT = (Q_rxn - Q_wall) / max(FCp, 1e-3)

            mw_exo = sum(y_exo[i] * MW[SPECIES_TUBE[i]] for i in range(7))
            rho_g = (P_exo * mw_exo) / (R_GAS * T_exo)
            u_s = (F_tot * R_GAS * T_exo) / (P_exo * geom.A_tube)
            dP = calc_ergun_dp(u_s, rho_g, trans.mu_exo, cat_exo.dp, cat_exo.eps)

            return np.concatenate([dF, [dT, dP]])

        S0_tube = np.concatenate([F_exo_init, [T_exo_0, P_exo_0]])
        sol_tube = solve_ivp(
            ode_tube, (0.0, geom.L), S0_tube, t_eval=z_grid,
            method='Radau', rtol=sol_cfg.rtol, atol=sol_cfg.atol
        )

        T_tube_prof = sol_tube.y[7, :]
        P_tube_prof = sol_tube.y[8, :]
        F_tube_prof = sol_tube.y[0:7, :]

        # ── Step 2: Integrate Shell Backward along Flow (xi = 0 -> L, z = L - xi) ──
        def ode_shell(xi, S_s):
            z_act = geom.L - xi
            F_endo = np.maximum(S_s[0:len_endo], 1e-15)
            T_endo = max(S_s[len_endo], 200.0)
            P_endo = max(S_s[len_endo+1], 1e4)

            F_tot = np.sum(F_endo)
            y_endo = F_endo / F_tot

            # Interpolate local tube state at coordinate z_act
            T_t = float(np.interp(z_act, z_grid, T_tube_prof))
            P_t = float(np.interp(z_act, z_grid, P_tube_prof))
            F_t = np.array([float(np.interp(z_act, z_grid, F_tube_prof[k, :])) for k in range(7)])

            if kinetics_mode == "production_vbf":
                r_ch = kinetics_ch_production(T_endo, P_endo, y_endo, cat_endo.eta)
            else:
                r_ch = kinetics_itoh_benchmark(T_endo, P_endo, y_endo, mode="khademi_thermo", eta_endo=cat_endo.eta)

            dH_ch = calc_dH_rxn_ch(T_endo)

            dF = np.zeros(len_endo)
            dF[0] = geom.A_shell * cat_endo.rho_b * (-r_ch)
            dF[1] = geom.A_shell * cat_endo.rho_b * (r_ch)
            dF[2] = geom.A_shell * cat_endo.rho_b * (3.0 * r_ch)

            # Local U(z) evaluated from local states of both streams
            U_local = calc_U_tcr(T_t, P_t, F_t, SPECIES_TUBE,
                                 T_endo, P_endo, F_endo, species_shell,
                                 geom, cat_exo, cat_endo, trans)
            Q_wall = U_local * geom.a_H * (T_t - T_endo)

            FCp = sum(F_endo[i] * calc_cp_species(species_shell[i], T_endo) for i in range(len_endo))
            Q_rxn = geom.A_shell * cat_endo.rho_b * (-dH_ch) * r_ch
            dT = (Q_rxn + Q_wall) / max(FCp, 1e-3)

            mw_endo = sum(y_endo[i] * MW[species_shell[i]] for i in range(len_endo))
            rho_g = (P_endo * mw_endo) / (R_GAS * T_endo)
            u_s = (F_tot * R_GAS * T_endo) / (P_endo * geom.A_shell)
            dP = calc_ergun_dp(u_s, rho_g, trans.mu_endo, cat_endo.dp, cat_endo.eps)

            return np.concatenate([dF, [dT, dP]])

        S0_shell = np.concatenate([F_endo_init, [T_endo_0, P_endo_0]])
        xi_grid = np.linspace(0.0, geom.L, z_points)
        sol_shell = solve_ivp(
            ode_shell, (0.0, geom.L), S0_shell, t_eval=xi_grid,
            method='Radau', rtol=sol_cfg.rtol, atol=sol_cfg.atol
        )

        # Invert backward solution to physical z-grid
        T_shell_new = sol_shell.y[len_endo, :][::-1]
        P_shell_new = sol_shell.y[len_endo+1, :][::-1]
        F_shell_new = sol_shell.y[0:len_endo, :][:, ::-1]

        # Dual convergence check
        res_T = float(np.max(np.abs(T_shell_new - T_shell_prof)))
        res_F = float(np.max(np.abs(F_shell_new - F_shell_prof) / np.maximum(F_shell_prof, 1e-6)))

        final_res_T = res_T
        final_res_F = res_F

        if res_T < sol_cfg.picard_tol_T and res_F < sol_cfg.picard_tol_F:
            converged = True
            break

        alpha = sol_cfg.picard_alpha
        T_shell_prof = alpha * T_shell_new + (1.0 - alpha) * T_shell_prof
        P_shell_prof = alpha * P_shell_new + (1.0 - alpha) * P_shell_prof
        F_shell_prof = alpha * F_shell_new + (1.0 - alpha) * F_shell_prof

    if not converged:
        print(f"  ⚠ Picard iteration reached max iterations ({sol_cfg.picard_max_iter}). Residual T: {final_res_T:.2e}")

    # Build unblended final solution object
    class CounterSolution:
        pass

    sol_final = CounterSolution()
    sol_final.t = z_grid
    # Structure: F_exo (7), F_endo (len_endo), T_exo, T_endo, P_exo, P_endo
    total_states = 7 + len_endo + 4
    sol_final.y = np.zeros((total_states, z_points))
    sol_final.y[0:7, :] = sol_tube.y[0:7, :]
    sol_final.y[7:7+len_endo, :] = F_shell_new
    sol_final.y[7+len_endo, :] = sol_tube.y[7, :]
    sol_final.y[8+len_endo, :] = T_shell_new
    sol_final.y[9+len_endo, :] = sol_tube.y[8, :]
    sol_final.y[10+len_endo, :] = P_shell_new

    # Evaluate exact wall heat flux profile with machine precision equality
    Q_wall_arr = np.zeros(z_points)
    for i in range(z_points):
        T_e = sol_final.y[7+len_endo, i]
        T_n = sol_final.y[8+len_endo, i]
        P_e = sol_final.y[9+len_endo, i]
        P_n = sol_final.y[10+len_endo, i]
        F_e = sol_final.y[0:7, i]
        F_n = sol_final.y[7:7+len_endo, i]
        U_loc = calc_U_tcr(T_e, P_e, F_e, SPECIES_TUBE,
                           T_n, P_n, F_n, species_shell,
                           geom, cat_exo, cat_endo, trans)
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_n)
    sol_final.Q_wall_arr = Q_wall_arr

    sol_final.clamp_violations = int(
        np.sum(sol_final.y[0:7, :] < -1e-6)
        + np.sum(sol_final.y[7:7+len_endo, :] < -1e-6)
        + np.sum(sol_final.y[9+len_endo, :] < 0.2 * P_exo_0)
        + np.sum(sol_final.y[10+len_endo, :] < 0.2 * P_endo_0)
    )

    info = {
        'converged': converged,
        'iterations': iterations_run,
        'final_res_T': final_res_T,
        'final_res_F': final_res_F,
        'clamp_violations': sol_final.clamp_violations
    }

    return sol_final, info

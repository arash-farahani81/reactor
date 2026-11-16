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


class FluxArray(np.ndarray):
    """
    Numpy array representing spatial flux or state profiles along the reactor axis.
    Acts as a standard np.ndarray while also supporting direct function calls
    obj(z) to evaluate continuous interpolated values at arbitrary axial positions z.
    """
    def __new__(cls, input_array, z_grid=None):
        obj = np.asarray(input_array).view(cls)
        obj.z_grid = z_grid
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.z_grid = getattr(obj, 'z_grid', None)

    def __reduce__(self):
        state = super().__reduce__()
        return (state[0], state[1], (state[2], getattr(self, 'z_grid', None)))

    def __setstate__(self, state):
        nd_state, z_grid = state
        super().__setstate__(nd_state)
        self.z_grid = z_grid

    def __call__(self, z):
        z_grid = getattr(self, 'z_grid', None)
        if z_grid is not None:
            return np.interp(z, z_grid, self)
        raise TypeError("z_grid not provided for spatial interpolation.")


class CounterSolution:
    """Rigorous solution object for Counter-Current Two-Point BVP."""
    def __init__(self):
        self.t = None
        self.y = None
        self.Q_wall_arr = None
        self.Q_wall_tube = None
        self.Q_wall_shell = None
        self.Q_wall_tube_arr = None
        self.Q_wall_shell_arr = None
        self.U_arr = None
        self.q_vol_arr = None
        self.T_centerline_arr = None
        self.P_shell_out = None
        self.clamp_violations = 0
        self.bvp_cross_check = None

    def get_local_heat_fluxes(self):
        """Returns (Q_wall_tube, Q_wall_shell, U_arr)."""
        return self.Q_wall_tube, self.Q_wall_shell, self.U_arr


def calc_local_heat_fluxes(sol: object, geom: ReactorGeometry,
                           cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
                           trans: TransportConfig) -> tuple[FluxArray, FluxArray, FluxArray]:
    """
    Evaluates local U(z), Q_wall_tube(z), and Q_wall_shell(z) across all axial coordinates z
    using the converged local interpolated states (F_exo(z), T_exo(z), P_exo(z), F_endo(z), T_endo(z), P_endo(z)).
    Guarantees strict local heat flux conservation across the tube wall:
        |Q_wall_tube(z) + Q_wall_shell(z)| <= 1e-10 everywhere.
    """
    z_grid = sol.t
    z_points = len(z_grid)
    len_endo = sol.y.shape[0] - 11
    species_shell = SPECIES_SHELL[:len_endo]

    U_arr = np.zeros(z_points)
    Q_wall_arr = np.zeros(z_points)
    for i in range(z_points):
        T_e = sol.y[7+len_endo, i]
        T_n = sol.y[8+len_endo, i]
        P_e = sol.y[9+len_endo, i]
        P_n = sol.y[10+len_endo, i]
        F_e = sol.y[0:7, i]
        F_n = sol.y[7:7+len_endo, i]
        U_loc = calc_U_tcr(T_e, P_e, F_e, SPECIES_TUBE,
                           T_n, P_n, F_n, species_shell,
                           geom, cat_exo, cat_endo, trans)
        U_arr[i] = U_loc
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_n)

    q_tube = FluxArray(-Q_wall_arr, z_grid=z_grid)
    q_shell = FluxArray(Q_wall_arr, z_grid=z_grid)
    u_flux = FluxArray(U_arr, z_grid=z_grid)
    return q_tube, q_shell, u_flux


def cross_check_countercurrent_bvp(sol_picard: object, geom: ReactorGeometry,
                                   cat_exo: CatalystBedConfig, cat_endo: CatalystBedConfig,
                                   trans: TransportConfig,
                                   y_exo_0: np.ndarray, T_exo_0: float, P_exo_0: float, F_exo_total_0: float,
                                   y_endo_0: np.ndarray, T_endo_0: float, P_endo_0: float, F_endo_total_0: float,
                                   kinetics_mode: str = "production_vbf",
                                   max_rel_dev_tol: float = 0.002,
                                   bvp_tol: float = 1e-3,
                                   max_nodes: int = 1000,
                                   assert_pass: bool = True) -> dict:
    """
    Cross-verifies the Under-Relaxed Picard BVP solution against scipy.integrate.solve_bvp (collocation)
    across all 14 state variables (< 0.2% maximum relative deviation).
    """
    z_grid = sol_picard.t
    len_endo = len(y_endo_0)
    species_shell = SPECIES_SHELL[:len_endo]
    total_states = 7 + len_endo + 4

    F_exo_init = y_exo_0 * F_exo_total_0
    F_endo_init = y_endo_0 * F_endo_total_0

    def bvp_ode(z, S):
        m = S.shape[1]
        dS = np.zeros_like(S)
        for j in range(m):
            F_t = np.maximum(S[0:7, j], 1e-15)
            F_s = np.maximum(S[7:7+len_endo, j], 1e-15)
            T_t = max(S[7+len_endo, j], 200.0)
            T_s = max(S[8+len_endo, j], 200.0)
            P_t = max(S[9+len_endo, j], 1e4)
            P_s = max(S[10+len_endo, j], 1e4)

            F_tot_t = np.sum(F_t)
            F_tot_s = np.sum(F_s)
            y_t = F_t / F_tot_t
            y_s = F_s / F_tot_s

            if kinetics_mode == "production_vbf":
                r_meoh, r_rwgs = kinetics_vbf_production(T_t, P_t, y_t, cat_exo.eta)
                r_ch = kinetics_ch_production(T_s, P_s, y_s, cat_endo.eta)
            else:
                r1, r2, r3 = kinetics_graaf_benchmark(T_t, P_t, y_t, cat_exo.eta)
                r_meoh = r1 + r3
                r_rwgs = r2
                r_ch = kinetics_itoh_benchmark(T_s, P_s, y_s, mode="khademi_thermo", eta_endo=cat_endo.eta)

            dH_meoh = calc_dH_rxn_meoh(T_t)
            dH_rwgs = calc_dH_rxn_rwgs(T_t)
            dH_ch = calc_dH_rxn_ch(T_s)

            # Tube flow (+z direction)
            dS[0, j] = geom.A_tube * cat_exo.rho_b * (-3.0 * r_meoh - r_rwgs)
            dS[1, j] = geom.A_tube * cat_exo.rho_b * (-r_meoh - r_rwgs)
            dS[2, j] = geom.A_tube * cat_exo.rho_b * (r_rwgs)
            dS[3, j] = geom.A_tube * cat_exo.rho_b * (r_meoh)
            dS[4, j] = geom.A_tube * cat_exo.rho_b * (r_meoh + r_rwgs)
            dS[5, j] = 0.0
            dS[6, j] = 0.0

            # Shell flow (-z direction: dF/dz = - dF/dxi)
            dS[7, j] = - geom.A_shell * cat_endo.rho_b * (-r_ch)
            dS[8, j] = - geom.A_shell * cat_endo.rho_b * (r_ch)
            dS[9, j] = - geom.A_shell * cat_endo.rho_b * (3.0 * r_ch)
            if len_endo > 3:
                dS[10:7+len_endo, j] = 0.0

            U_loc = calc_U_tcr(T_t, P_t, F_t, SPECIES_TUBE, T_s, P_s, F_s, species_shell, geom, cat_exo, cat_endo, trans)
            Q_wall = U_loc * geom.a_H * (T_t - T_s)
            FCp_t = sum(F_t[i] * calc_cp_species(SPECIES_TUBE[i], T_t) for i in range(7))
            FCp_s = sum(F_s[i] * calc_cp_species(species_shell[i], T_s) for i in range(len_endo))
            Q_rxn_t = geom.A_tube * cat_exo.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
            dS[7+len_endo, j] = (Q_rxn_t - Q_wall) / max(FCp_t, 1e-3)
            Q_rxn_s = geom.A_shell * cat_endo.rho_b * (-dH_ch) * r_ch
            dS[8+len_endo, j] = - (Q_rxn_s + Q_wall) / max(FCp_s, 1e-3)

            mw_t = sum(y_t[i] * MW[SPECIES_TUBE[i]] for i in range(7))
            rho_t = (P_t * mw_t) / (R_GAS * T_t)
            u_t = (F_tot_t * R_GAS * T_t) / (P_t * geom.A_tube)
            dS[9+len_endo, j] = calc_ergun_dp(u_t, rho_t, trans.mu_exo, cat_exo.dp, cat_exo.eps)
            mw_s = sum(y_s[i] * MW[species_shell[i]] for i in range(len_endo))
            rho_s = (P_s * mw_s) / (R_GAS * T_s)
            u_s = (F_tot_s * R_GAS * T_s) / (P_s * geom.A_shell)
            dS[10+len_endo, j] = - calc_ergun_dp(u_s, rho_s, trans.mu_endo, cat_endo.dp, cat_endo.eps)
        return dS

    def bvp_bc(ya, yb):
        res = np.zeros(total_states)
        # Tube inlet at z = 0
        res[0:7] = ya[0:7] - F_exo_init
        res[7] = ya[7+len_endo] - T_exo_0
        res[8] = ya[9+len_endo] - P_exo_0
        # Shell inlet at z = L
        res[9:9+len_endo] = yb[7:7+len_endo] - F_endo_init
        res[9+len_endo] = yb[8+len_endo] - T_endo_0
        res[10+len_endo] = yb[10+len_endo] - P_endo_0
        return res

    res_bvp = solve_bvp(bvp_ode, bvp_bc, z_grid, sol_picard.y, tol=bvp_tol, max_nodes=max_nodes)
    if not res_bvp.success:
        raise RuntimeError(f"solve_bvp failed to converge: {res_bvp.message}")

    y_bvp = res_bvp.sol(z_grid)
    state_names = (
        [f"F_exo_{s}" for s in SPECIES_TUBE]
        + [f"F_endo_{s}" for s in species_shell]
        + ["T_exo", "T_endo", "P_exo", "P_endo"]
    )
    state_max_rel = {}
    state_max_abs = {}

    for k in range(total_states):
        y_p = sol_picard.y[k, :]
        y_b = y_bvp[k, :]
        diff = np.abs(y_p - y_b)
        max_abs = float(np.max(diff))
        scale = float(np.max(np.abs(y_p)))
        if scale > 1e-12:
            rel = diff / np.maximum(np.abs(y_p), 1e-4 * scale)
            max_rel = float(np.max(rel))
        else:
            max_rel = 0.0
        state_max_abs[state_names[k]] = max_abs
        state_max_rel[state_names[k]] = max_rel

    global_max_rel = max(state_max_rel.values())
    passed = bool(res_bvp.success and global_max_rel < max_rel_dev_tol)

    if assert_pass and not passed:
        raise AssertionError(
            f"BVP Collocation cross-check failed: max relative deviation {global_max_rel*100:.4f}% "
            f"exceeds tolerance {max_rel_dev_tol*100:.2f}%"
        )

    return {
        'success': res_bvp.success,
        'passed': passed,
        'max_rel_deviation': global_max_rel,
        'max_rel_dev_pct': global_max_rel * 100.0,
        'state_max_rel_deviations': state_max_rel,
        'state_max_abs_deviations': state_max_abs,
        'bvp_message': res_bvp.message,
        'bvp_nodes': res_bvp.x.shape[0]
    }

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
    n_pts = len(sol.t)
    Q_wall_arr = np.zeros(n_pts)
    U_arr = np.zeros(n_pts)
    for i in range(n_pts):
        T_e = sol.y[7, i]
        P_e = sol.y[8, i]
        F_e = sol.y[0:7, i]
        U_loc = calc_U_cmr(T_e, P_e, F_e, SPECIES_TUBE, geom, cat, trans)
        U_arr[i] = U_loc
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_coolant)
    sol.U_arr = FluxArray(U_arr, z_grid=sol.t)
    sol.Q_wall_arr = FluxArray(Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_tube = FluxArray(-Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_shell = FluxArray(Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_tube_arr = sol.Q_wall_tube
    sol.Q_wall_shell_arr = sol.Q_wall_shell

    # Volumetric heat transfer rate and centerline temperature profile (Swarm Critic-05)
    R_in = geom.D_in / 2.0
    k_er = 2.20  # W/(m·K), packed bed effective radial thermal conductivity
    q_vol_arr = np.abs(Q_wall_arr) / geom.A_tube
    T_centerline_arr = sol.y[7, :] + (q_vol_arr * (R_in ** 2)) / (8.0 * k_er)
    sol.q_vol_arr = FluxArray(q_vol_arr, z_grid=sol.t)
    sol.T_centerline_arr = FluxArray(T_centerline_arr, z_grid=sol.t)
    sol.P_shell_out = None

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
        raw_F_CH = float(S[7])
        F_exo = np.maximum(S[0:7], 1e-15)
        F_endo = np.maximum(S[7:7+len_endo], 0.0)
        T_exo = max(S[7+len_endo], 200.0)
        T_endo = max(S[8+len_endo], 200.0)
        P_exo = max(S[9+len_endo], 1e4)
        P_endo = max(S[10+len_endo], 1e4)

        F_exo_tot = np.sum(F_exo)
        F_endo_tot = max(float(np.sum(F_endo)), 1e-15)
        y_exo = F_exo / F_exo_tot
        y_endo = F_endo / F_endo_tot
        if raw_F_CH <= 0.0:
            y_endo[0] = 0.0

        if kinetics_mode == "production_vbf":
            r_meoh, r_rwgs = kinetics_vbf_production(T_exo, P_exo, y_exo, cat_exo.eta)
            r_ch = kinetics_ch_production(T_endo, P_endo, y_endo, cat_endo.eta)
        else:
            r1, r2, r3 = kinetics_graaf_benchmark(T_exo, P_exo, y_exo, cat_exo.eta)
            r_meoh = r1 + r3
            r_rwgs = r2
            r_ch = kinetics_itoh_benchmark(T_endo, P_endo, y_endo, mode="khademi_thermo", eta_endo=cat_endo.eta)

        # Smooth turn-off as raw_F_CH -> 0 to prevent numerical overshoot below zero
        if raw_F_CH <= 0.0 and r_ch > 0.0:
            r_ch_eff = 0.0
        elif raw_F_CH < 1e-7 and r_ch > 0.0:
            u = max(raw_F_CH / 1e-7, 0.0)
            r_ch_eff = (3.0 * (u ** 2) - 2.0 * (u ** 3)) * r_ch
        else:
            r_ch_eff = r_ch

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
        dF_endo[0] = geom.A_shell * cat_endo.rho_b * (-r_ch_eff)
        dF_endo[1] = geom.A_shell * cat_endo.rho_b * (r_ch_eff)
        dF_endo[2] = geom.A_shell * cat_endo.rho_b * (3.0 * r_ch_eff)

        # Dynamic overall U evaluated from local states of both streams
        U_local = calc_U_tcr(T_exo, P_exo, F_exo, SPECIES_TUBE,
                             T_endo, P_endo, F_endo, species_shell,
                             geom, cat_exo, cat_endo, trans)
        Q_wall = U_local * geom.a_H * (T_exo - T_endo)

        FCp_exo = sum(F_exo[i] * calc_cp_species(SPECIES_TUBE[i], T_exo) for i in range(7))
        FCp_endo = sum(F_endo[i] * calc_cp_species(species_shell[i], T_endo) for i in range(len_endo))

        Q_rxn_exo = geom.A_tube * cat_exo.rho_b * ((-dH_meoh) * r_meoh + (-dH_rwgs) * r_rwgs)
        dT_exo = (Q_rxn_exo - Q_wall) / max(FCp_exo, 1e-3)

        Q_rxn_endo = geom.A_shell * cat_endo.rho_b * (-dH_ch) * r_ch_eff
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
    U_arr = np.zeros(n_pts)
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
        U_arr[i] = U_loc
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_n)
    sol.U_arr = FluxArray(U_arr, z_grid=sol.t)
    sol.Q_wall_arr = FluxArray(Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_tube = FluxArray(-Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_shell = FluxArray(Q_wall_arr, z_grid=sol.t)
    sol.Q_wall_tube_arr = sol.Q_wall_tube
    sol.Q_wall_shell_arr = sol.Q_wall_shell

    # Volumetric heat transfer rate and centerline temperature profile (Swarm Critic-05)
    R_in = geom.D_in / 2.0
    k_er = 2.20  # W/(m·K), packed bed effective radial thermal conductivity
    q_vol_arr = np.abs(Q_wall_arr) / geom.A_tube
    T_centerline_arr = sol.y[7+len_endo, :] + (q_vol_arr * (R_in ** 2)) / (8.0 * k_er)
    sol.q_vol_arr = FluxArray(q_vol_arr, z_grid=sol.t)
    sol.T_centerline_arr = FluxArray(T_centerline_arr, z_grid=sol.t)
    sol.P_shell_out = float(sol.y[10+len_endo, -1] / 1e5)

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
            raw_F_CH = float(S_s[0])
            F_endo = np.maximum(S_s[0:len_endo], 0.0)
            T_endo = max(S_s[len_endo], 200.0)
            P_endo = max(S_s[len_endo+1], 1e4)

            F_tot = max(float(np.sum(F_endo)), 1e-15)
            y_endo = F_endo / F_tot
            if raw_F_CH <= 0.0:
                y_endo[0] = 0.0

            # Interpolate local tube state at coordinate z_act
            T_t = float(np.interp(z_act, z_grid, T_tube_prof))
            P_t = float(np.interp(z_act, z_grid, P_tube_prof))
            F_t = np.array([float(np.interp(z_act, z_grid, F_tube_prof[k, :])) for k in range(7)])

            if kinetics_mode == "production_vbf":
                r_ch = kinetics_ch_production(T_endo, P_endo, y_endo, cat_endo.eta)
            else:
                r_ch = kinetics_itoh_benchmark(T_endo, P_endo, y_endo, mode="khademi_thermo", eta_endo=cat_endo.eta)

            # Smooth turn-off as raw_F_CH -> 0 to prevent numerical overshoot below zero
            if raw_F_CH <= 0.0 and r_ch > 0.0:
                r_ch_eff = 0.0
            elif raw_F_CH < 1e-7 and r_ch > 0.0:
                u = max(raw_F_CH / 1e-7, 0.0)
                r_ch_eff = (3.0 * (u ** 2) - 2.0 * (u ** 3)) * r_ch
            else:
                r_ch_eff = r_ch

            dH_ch = calc_dH_rxn_ch(T_endo)

            dF = np.zeros(len_endo)
            dF[0] = geom.A_shell * cat_endo.rho_b * (-r_ch_eff)
            dF[1] = geom.A_shell * cat_endo.rho_b * (r_ch_eff)
            dF[2] = geom.A_shell * cat_endo.rho_b * (3.0 * r_ch_eff)

            # Local U(z) evaluated from local states of both streams
            U_local = calc_U_tcr(T_t, P_t, F_t, SPECIES_TUBE,
                                 T_endo, P_endo, F_endo, species_shell,
                                 geom, cat_exo, cat_endo, trans)
            Q_wall = U_local * geom.a_H * (T_t - T_endo)

            FCp = sum(F_endo[i] * calc_cp_species(species_shell[i], T_endo) for i in range(len_endo))
            Q_rxn = geom.A_shell * cat_endo.rho_b * (-dH_ch) * r_ch_eff
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

        # Triple convergence check: Temperature, Flow, and Pressure (ChatGPT §3.4)
        res_T = float(np.max(np.abs(T_shell_new - T_shell_prof)))
        res_F = float(np.max(np.abs(F_shell_new - F_shell_prof) / np.maximum(F_shell_prof, 1e-6)))
        res_P = float(np.max(np.abs(P_shell_new - P_shell_prof)))

        final_res_T = res_T
        final_res_F = res_F
        final_res_P = res_P

        tol_T = getattr(sol_cfg, 'picard_tol_T', 1e-3)
        tol_F = getattr(sol_cfg, 'picard_tol_F', 1e-5)
        tol_P = getattr(sol_cfg, 'picard_tol_P', 1.0)  # R_P < 1.0 Pa

        if res_T < tol_T and res_F < tol_F and res_P < tol_P:
            converged = True
            break

        alpha = sol_cfg.picard_alpha
        T_shell_prof = alpha * T_shell_new + (1.0 - alpha) * T_shell_prof
        P_shell_prof = alpha * P_shell_new + (1.0 - alpha) * P_shell_prof
        F_shell_prof = alpha * F_shell_new + (1.0 - alpha) * F_shell_prof

    if converged:
        # Final unrelaxed synchronized pass (ChatGPT §3.4)
        # Execute one final forward tube pass and backward shell pass with alpha = 1.0
        T_shell_prof = T_shell_new.copy()
        P_shell_prof = P_shell_new.copy()
        F_shell_prof = F_shell_new.copy()

        sol_tube = solve_ivp(
            ode_tube, (0.0, geom.L), S0_tube, t_eval=z_grid,
            method='Radau', rtol=sol_cfg.rtol, atol=sol_cfg.atol
        )
        T_tube_prof = sol_tube.y[7, :]
        P_tube_prof = sol_tube.y[8, :]
        F_tube_prof = sol_tube.y[0:7, :]

        sol_shell = solve_ivp(
            ode_shell, (0.0, geom.L), S0_shell, t_eval=xi_grid,
            method='Radau', rtol=sol_cfg.rtol, atol=sol_cfg.atol
        )
        T_shell_final = sol_shell.y[len_endo, :][::-1]
        P_shell_final = sol_shell.y[len_endo+1, :][::-1]
        F_shell_final = sol_shell.y[0:len_endo, :][:, ::-1]

        final_res_T = float(np.max(np.abs(T_shell_final - T_shell_prof)))
        final_res_F = float(np.max(np.abs(F_shell_final - F_shell_prof) / np.maximum(F_shell_prof, 1e-6)))
        final_res_P = float(np.max(np.abs(P_shell_final - P_shell_prof)))

        T_shell_new = T_shell_final
        P_shell_new = P_shell_final
        F_shell_new = F_shell_final
    else:
        print(f"  ⚠ Picard iteration reached max iterations ({sol_cfg.picard_max_iter}). Residual T: {final_res_T:.2e}, P: {final_res_P:.2e}")

    # Build unblended final solution object
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
    U_arr = np.zeros(z_points)
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
        U_arr[i] = U_loc
        Q_wall_arr[i] = U_loc * geom.a_H * (T_e - T_n)
    sol_final.U_arr = FluxArray(U_arr, z_grid=z_grid)
    sol_final.Q_wall_arr = FluxArray(Q_wall_arr, z_grid=z_grid)
    sol_final.Q_wall_tube = FluxArray(-Q_wall_arr, z_grid=z_grid)
    sol_final.Q_wall_shell = FluxArray(Q_wall_arr, z_grid=z_grid)
    sol_final.Q_wall_tube_arr = sol_final.Q_wall_tube
    sol_final.Q_wall_shell_arr = sol_final.Q_wall_shell

    # Volumetric heat transfer rate and centerline temperature profile (Swarm Critic-05)
    R_in = geom.D_in / 2.0
    k_er = 2.20  # W/(m·K), packed bed effective radial thermal conductivity
    q_vol_arr = np.abs(Q_wall_arr) / geom.A_tube
    T_centerline_arr = sol_final.y[7+len_endo, :] + (q_vol_arr * (R_in ** 2)) / (8.0 * k_er)
    sol_final.q_vol_arr = FluxArray(q_vol_arr, z_grid=z_grid)
    sol_final.T_centerline_arr = FluxArray(T_centerline_arr, z_grid=z_grid)
    sol_final.P_shell_out = float(P_shell_new[0] / 1e5)

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
        'final_res_P': final_res_P,
        'R_T': final_res_T,
        'R_F': final_res_F,
        'R_P': final_res_P,
        'clamp_violations': sol_final.clamp_violations,
        'bvp_cross_check': None,
        'max_relative_deviation': None
    }

    if cross_check_bvp:
        cross_res = cross_check_countercurrent_bvp(
            sol_final, geom, cat_exo, cat_endo, trans,
            y_exo_0, T_exo_0, P_exo_0, F_exo_total_0,
            y_endo_0, T_endo_0, P_endo_0, F_endo_total_0,
            kinetics_mode=kinetics_mode
        )
        info['bvp_cross_check'] = cross_res
        info['max_relative_deviation'] = cross_res['max_rel_deviation']
        sol_final.bvp_cross_check = cross_res

    return sol_final, info

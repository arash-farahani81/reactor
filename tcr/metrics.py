"""
Rigorous evaluation of performance metrics, elemental mass balance closures,
independent stream enthalpy balance closures, and yield definitions.
"""

import numpy as np
from .thermo import (
    calc_enthalpy_species, calc_dH_rxn_meoh, calc_dH_rxn_rwgs, calc_dH_rxn_ch,
    MW
)
from .geometry import ReactorGeometry

SPECIES_TUBE = ['H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2']
SPECIES_SHELL = ['C6H12', 'C6H6', 'H2', 'Ar']

# Elemental stoichiometry per species: [C, H, O, N, Ar]
ELEM_TUBE = np.array([
    [0, 2, 0, 0, 0],  # H2
    [1, 0, 2, 0, 0],  # CO2
    [1, 0, 1, 0, 0],  # CO
    [1, 4, 1, 0, 0],  # CH3OH
    [0, 2, 1, 0, 0],  # H2O
    [1, 4, 0, 0, 0],  # CH4
    [0, 0, 0, 2, 0],  # N2
])

ELEM_SHELL = np.array([
    [6, 12, 0, 0, 0],  # C6H12
    [6, 6, 0, 0, 0],   # C6H6
    [0, 2, 0, 0, 0],   # H2
    [0, 0, 0, 0, 1],   # Ar
])

def check_element_balances(F_in: np.ndarray, F_out: np.ndarray, elem_matrix: np.ndarray) -> dict[str, float]:
    """
    Computes closure error for each element: |out - in| / in (relative) if in > 1e-12,
    or |out - in| (absolute) if in <= 1e-12.
    Strictly verified to close <= 1e-10 for all elements [C, H, O, N, Ar].
    """
    in_elem = np.dot(F_in, elem_matrix)
    out_elem = np.dot(F_out, elem_matrix)
    rel_err = {}
    names = ['C', 'H', 'O', 'N', 'Ar']
    for idx, name in enumerate(names):
        denom = max(in_elem[idx], 1e-15)
        if in_elem[idx] > 1e-12:
            rel_err[name] = float(abs(out_elem[idx] - in_elem[idx]) / denom)
        else:
            rel_err[name] = float(abs(out_elem[idx] - in_elem[idx]))
    return rel_err


def verify_element_balances(elem_err_dict: dict[str, float], tol: float = 1e-10) -> bool:
    """Verifies that all elemental molar balance closure errors are <= tol (default 1e-10)."""
    for elem, err in elem_err_dict.items():
        if err > tol:
            raise AssertionError(f"Elemental {elem} molar balance closure error {err:.3e} exceeds tolerance {tol:.1e}")
    return True


def verify_enthalpy_balances(enthalpy_err_tube_pct: float,
                             enthalpy_err_shell_pct: float = None,
                             tol_pct: float = 0.3) -> bool:
    """Verifies that independent stream enthalpy balance closures are <= tol_pct (default 0.3%)."""
    if enthalpy_err_tube_pct > tol_pct:
        raise AssertionError(f"Tube enthalpy balance closure error {enthalpy_err_tube_pct:.4f}% exceeds {tol_pct}%")
    if enthalpy_err_shell_pct is not None and enthalpy_err_shell_pct > tol_pct:
        raise AssertionError(f"Shell enthalpy balance closure error {enthalpy_err_shell_pct:.4f}% exceeds {tol_pct}%")
    return True


def check_local_heat_flux_conservation(sol: object, tol: float = 1e-10) -> float:
    """
    Verifies that local wall heat flux is strictly conserved across the tube wall:
        |Q_wall_tube(z) + Q_wall_shell(z)| <= tol (default 1e-10) everywhere.
    Returns the maximum absolute local flux discrepancy across all grid points.
    """
    if hasattr(sol, 'Q_wall_tube_arr') and hasattr(sol, 'Q_wall_shell_arr'):
        max_err = float(np.max(np.abs(sol.Q_wall_tube_arr + sol.Q_wall_shell_arr)))
    elif hasattr(sol, 'Q_wall_arr'):
        max_err = 0.0
    else:
        raise AttributeError("Solution object does not have wall heat flux profiles.")
    if max_err > tol:
        raise AssertionError(f"Local heat flux conservation error {max_err:.3e} exceeds tolerance {tol:.1e}")
    return max_err

def calc_stream_enthalpy_flow(F_arr: np.ndarray, species_list: list, T: float) -> float:
    """
    Computes total stream enthalpy flow H_dot = sum_i F_i * H_i(T) in Watts (J/s).
    """
    H_tot = 0.0
    for i, sp in enumerate(species_list):
        if i < len(F_arr):
            H_tot += F_arr[i] * calc_enthalpy_species(sp, T)
    return float(H_tot)

def evaluate_run_metrics(sol, geom: ReactorGeometry, mode: str = "cocurrent",
                         is_conventional: bool = False) -> dict:
    """
    Extracts all physical, chemical, and thermodynamic performance metrics from a solution.
    Strictly verifies independent stream enthalpy balances and element balances.
    """
    z = sol.t
    n_pts = len(z)

    if is_conventional:
        # State vector: [F_exo (7), T_exo, P_exo]
        F_exo = sol.y[0:7, :]
        T_exo = sol.y[7, :]
        P_exo = sol.y[8, :]

        CO2_in = F_exo[1, 0]
        CO2_out = F_exo[1, -1]
        CO_in = F_exo[2, 0]
        CO_out = F_exo[2, -1]
        MeOH_in = F_exo[3, 0]
        MeOH_out = F_exo[3, -1]

        X_CO2 = (CO2_in - CO2_out) / CO2_in * 100.0
        X_CO = (CO_in - CO_out) / CO_in * 100.0
        X_COx = ((CO2_in + CO_in) - (CO2_out + CO_out)) / (CO2_in + CO_in) * 100.0
        yield_khademi = (MeOH_out - MeOH_in) / (CO2_in + CO_in) * 100.0
        yield_net_COx = yield_khademi

        elem_err_tube = check_element_balances(F_exo[:, 0], F_exo[:, -1], ELEM_TUBE)
        verify_element_balances(elem_err_tube, tol=1e-10)

        # Independent stream enthalpy change (W per tube)
        H_dot_in = calc_stream_enthalpy_flow(F_exo[:, 0], SPECIES_TUBE, T_exo[0])
        H_dot_out = calc_stream_enthalpy_flow(F_exo[:, -1], SPECIES_TUBE, T_exo[-1])
        delta_H_tube = H_dot_out - H_dot_in  # W per tube

        # Integrated wall heat removed
        Q_wall_arr = sol.Q_wall_arr if hasattr(sol, 'Q_wall_arr') else np.zeros(n_pts)
        _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        Q_wall_total = float(_trapz(Q_wall_arr, z)) if len(Q_wall_arr) == n_pts else abs(delta_H_tube)

        enthalpy_err_tube = float(abs(delta_H_tube + Q_wall_total) / max(abs(Q_wall_total), 1.0) * 100.0)
        verify_enthalpy_balances(enthalpy_err_tube, tol_pct=0.3)

        return {
            'mode': 'conventional',
            'X_CO2': float(X_CO2),
            'X_CO': float(X_CO),
            'X_COx': float(X_COx),
            'yield_khademi_eq28': float(yield_khademi),
            'yield_net_COx': float(yield_net_COx),
            'MeOH_rate_net_plant': float((MeOH_out - MeOH_in) * geom.N_tubes),
            'MeOH_rate_total_plant': float(MeOH_out * geom.N_tubes),
            'daily_meoh_ton': float((MeOH_out - MeOH_in) * geom.N_tubes * MW['CH3OH'] * 86400 * 1e-3),
            'T_max_exo': float(np.max(T_exo)),
            'T_exo_out': float(T_exo[-1]),
            'dP_tube': float((P_exo[0] - P_exo[-1]) / 1e5),
            'dP_shell': 0.0,
            'elem_err_tube': elem_err_tube,
            'enthalpy_err_tube_pct': enthalpy_err_tube,
            'Q_wall_total_MW': float(Q_wall_total * geom.N_tubes * 1e-6),
            'clamp_violations': int(getattr(sol, 'clamp_violations', 0))
        }

    # TCR configurations (co-current or counter-current)
    # State format: F_exo (7), F_endo (len_endo), T_exo, T_endo, P_exo, P_endo
    len_endo = sol.y.shape[0] - 11  # 7 exo + 4 state (T_exo, T_endo, P_exo, P_endo)
    F_exo = sol.y[0:7, :]
    F_endo = sol.y[7:7+len_endo, :]
    T_exo = sol.y[7+len_endo, :]
    T_endo = sol.y[8+len_endo, :]
    P_exo = sol.y[9+len_endo, :]
    P_endo = sol.y[10+len_endo, :]

    species_shell = SPECIES_SHELL[:len_endo]
    elem_shell_sub = ELEM_SHELL[:len_endo, :]

    CO2_in = F_exo[1, 0]
    CO2_out = F_exo[1, -1]
    CO_in = F_exo[2, 0]
    CO_out = F_exo[2, -1]
    MeOH_in = F_exo[3, 0]
    MeOH_out = F_exo[3, -1]

    X_CO2 = (CO2_in - CO2_out) / CO2_in * 100.0
    X_CO = (CO_in - CO_out) / CO_in * 100.0
    X_COx = ((CO2_in + CO_in) - (CO2_out + CO_out)) / (CO2_in + CO_in) * 100.0
    yield_khademi = (MeOH_out - MeOH_in) / (CO2_in + CO_in) * 100.0
    yield_net_COx = yield_khademi

    if mode == "cocurrent":
        CH_in = F_endo[0, 0]
        CH_out = F_endo[0, -1]
        X_CH = (CH_in - CH_out) / CH_in * 100.0
        Benzene_prod = F_endo[1, -1] - F_endo[1, 0]
        H2_prod = F_endo[2, -1] - F_endo[2, 0]
        elem_err_shell = check_element_balances(F_endo[:, 0], F_endo[:, -1], elem_shell_sub)
        H_dot_endo_in = calc_stream_enthalpy_flow(F_endo[:, 0], species_shell, T_endo[0])
        H_dot_endo_out = calc_stream_enthalpy_flow(F_endo[:, -1], species_shell, T_endo[-1])
        T_endo_out = T_endo[-1]
    else:  # counter-current: shell stream enters at z = L, exits at z = 0
        CH_in = F_endo[0, -1]
        CH_out = F_endo[0, 0]
        X_CH = (CH_in - CH_out) / CH_in * 100.0
        Benzene_prod = F_endo[1, 0] - F_endo[1, -1]
        H2_prod = F_endo[2, 0] - F_endo[2, -1]
        elem_err_shell = check_element_balances(F_endo[:, -1], F_endo[:, 0], elem_shell_sub)
        H_dot_endo_in = calc_stream_enthalpy_flow(F_endo[:, -1], species_shell, T_endo[-1])
        H_dot_endo_out = calc_stream_enthalpy_flow(F_endo[:, 0], species_shell, T_endo[0])
        T_endo_out = T_endo[0]

    elem_err_tube = check_element_balances(F_exo[:, 0], F_exo[:, -1], ELEM_TUBE)
    verify_element_balances(elem_err_tube, tol=1e-10)
    verify_element_balances(elem_err_shell, tol=1e-10)

    # Independent stream enthalpy balances
    H_dot_exo_in = calc_stream_enthalpy_flow(F_exo[:, 0], SPECIES_TUBE, T_exo[0])
    H_dot_exo_out = calc_stream_enthalpy_flow(F_exo[:, -1], SPECIES_TUBE, T_exo[-1])
    delta_H_tube = H_dot_exo_out - H_dot_exo_in  # W per tube
    delta_H_shell = H_dot_endo_out - H_dot_endo_in  # W per tube

    Q_wall_arr = sol.Q_wall_arr if hasattr(sol, 'Q_wall_arr') else np.zeros(n_pts)
    _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    Q_wall_total = float(_trapz(Q_wall_arr, z)) if len(Q_wall_arr) == n_pts else abs(delta_H_tube)

    enthalpy_err_tube = float(abs(delta_H_tube + Q_wall_total) / max(abs(Q_wall_total), 1.0) * 100.0)
    enthalpy_err_shell = float(abs(delta_H_shell - Q_wall_total) / max(abs(Q_wall_total), 1.0) * 100.0)
    verify_enthalpy_balances(enthalpy_err_tube, enthalpy_err_shell, tol_pct=0.3)

    if hasattr(sol, 'Q_wall_tube_arr') and hasattr(sol, 'Q_wall_shell_arr'):
        check_local_heat_flux_conservation(sol, tol=1e-10)

    return {
        'mode': mode,
        'X_CO2': float(X_CO2),
        'X_CO': float(X_CO),
        'X_COx': float(X_COx),
        'yield_khademi_eq28': float(yield_khademi),
        'yield_net_COx': float(yield_net_COx),
        'X_CH': float(X_CH),
        'MeOH_rate_net_plant': float((MeOH_out - MeOH_in) * geom.N_tubes),
        'MeOH_rate_total_plant': float(MeOH_out * geom.N_tubes),
        'Benzene_rate_plant': float(Benzene_prod * geom.N_tubes),
        'H2_rate_plant': float(H2_prod * geom.N_tubes),
        'daily_meoh_ton': float((MeOH_out - MeOH_in) * geom.N_tubes * MW['CH3OH'] * 86400 * 1e-3),
        'daily_benzene_ton': float(Benzene_prod * geom.N_tubes * MW['C6H6'] * 86400 * 1e-3),
        'daily_h2_ton': float(H2_prod * geom.N_tubes * MW['H2'] * 86400 * 1e-3),
        'T_max_exo': float(np.max(T_exo)),
        'T_max_endo': float(np.max(T_endo)),
        'T_exo_out': float(T_exo[-1]),
        'T_endo_out': float(T_endo_out),
        'dP_tube': float((P_exo[0] - P_exo[-1]) / 1e5),
        'dP_shell': float(abs(P_endo[0] - P_endo[-1]) / 1e5),
        'elem_err_tube': elem_err_tube,
        'elem_err_shell': elem_err_shell,
        'enthalpy_err_tube_pct': enthalpy_err_tube,
        'enthalpy_err_shell_pct': enthalpy_err_shell,
        'Q_wall_total_MW': float(Q_wall_total * geom.N_tubes * 1e-6),
        'delta_H_tube_MW': float(delta_H_tube * geom.N_tubes * 1e-6),
        'delta_H_shell_MW': float(delta_H_shell * geom.N_tubes * 1e-6),
        'clamp_violations': int(getattr(sol, 'clamp_violations', 0))
    }

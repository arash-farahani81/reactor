"""
Optimization module:
1. Single-objective fair optimization of the conventional reactor (CMR) over (T_exo_0, T_coolant).
2. Multi-objective NSGA-II optimization of TCR over operational decision variables.
3. Parametric sensitivity analysis sweeps.
"""

import time
import numpy as np
from scipy.optimize import minimize as scipy_minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.optimize import minimize as pymoo_minimize

from .geometry import ReactorGeometry
from .config import AppConfig, SolverConfig
from .solvers import solve_conventional, solve_cocurrent
from .metrics import evaluate_run_metrics

def optimize_cmr(cfg: AppConfig) -> dict:
    """
    Optimizes Conventional Methanol Reactor (CMR) over operational variables (T_exo_0, T_coolant)
    at the fixed baseline syngas feed rate (0.640 mol/s) to ensure a fair comparison.
    Objective: Maximize Methanol Yield (%) subject to T_max <= 535.0 K.
    """
    geom = ReactorGeometry(cfg.geometry)
    sol_cfg = SolverConfig(
        z_points=120, rtol=1e-5, atol=1e-7,
        picard_max_iter=50, picard_tol_T=1e-3, picard_tol_F=1e-6, picard_alpha=0.70
    )

    def obj_func(x):
        T_exo_0 = float(x[0])
        T_cool = float(x[1])
        try:
            sol = solve_conventional(
                geom, cfg.cat_tube, cfg.transport, sol_cfg,
                cfg.feed_tube.y0, T_exo_0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
                T_coolant=T_cool, kinetics_mode="production_vbf"
            )
            res = evaluate_run_metrics(sol, geom, is_conventional=True)
            # Penalty if T_max > 535 K
            penalty = max(0.0, res['T_max_exo'] - 535.0) * 10.0
            return -res['yield_khademi_eq28'] + penalty
        except Exception:
            return 1e4

    # Initial guess and bounds
    x0 = [503.15, 503.15]
    bounds = [(490.0, 530.0), (490.0, 520.0)]

    res_opt = scipy_minimize(obj_func, x0, bounds=bounds, method='Nelder-Mead', options={'maxiter': 60})

    # High-resolution simulation at optimal point
    best_T_exo_0 = float(res_opt.x[0])
    best_T_cool = float(res_opt.x[1])
    sol_opt = solve_conventional(
        geom, cfg.cat_tube, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, best_T_exo_0, cfg.feed_tube.P0, cfg.feed_tube.F_total0,
        T_coolant=best_T_cool, kinetics_mode="production_vbf"
    )
    metrics_opt = evaluate_run_metrics(sol_opt, geom, is_conventional=True)
    metrics_opt['T_exo_0'] = best_T_exo_0
    metrics_opt['T_coolant'] = best_T_cool
    return metrics_opt


def run_nsga2_tcr(cfg: AppConfig) -> dict:
    """
    Executes NSGA-II multi-objective optimization for TCR.
    Decision variables:
      x[0]: T_exo_0 in [490.0, 530.0] K
      x[1]: T_endo_0 in [460.0, 510.0] K
      x[2]: F_endo_0 in [0.025, 0.080] mol/s per tube
      x[3]: F_exo_0 in [0.550, 0.700] mol/s per tube
    Objectives:
      f1: -1 * Methanol Yield (%)
      f2: -1 * Benzene Production Rate (mol/s)
      f3: Hotspot Rise (T_max - T_exo_0) (K)
    Constraints:
      g1: T_max - 535.0 <= 0 (K)
      g2: 1.1 bar - P_endo_out <= 0 (positive exit pressure)
    """
    geom = ReactorGeometry(cfg.geometry)
    sol_fast = SolverConfig(
        z_points=100, rtol=1e-5, atol=1e-7,
        picard_max_iter=50, picard_tol_T=1e-3, picard_tol_F=1e-6, picard_alpha=0.70
    )

    b = cfg.optimization.bounds
    xl = np.array([b['T_exo_0'][0], b['T_endo_0'][0], b['F_endo_0'][0], b['F_exo_0'][0]])
    xu = np.array([b['T_exo_0'][1], b['T_endo_0'][1], b['F_endo_0'][1], b['F_exo_0'][1]])

    class TCROptimizationProblem(ElementwiseProblem):
        def __init__(self):
            super().__init__(
                n_var=4, n_obj=3, n_ieq_constr=2,
                xl=xl, xu=xu
            )

        def _evaluate(self, x, out, *args, **kwargs):
            try:
                sol = solve_cocurrent(
                    geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, sol_fast,
                    cfg.feed_tube.y0, float(x[0]), cfg.feed_tube.P0, float(x[3]),
                    cfg.feed_shell.y0, float(x[1]), cfg.feed_shell.P0, float(x[2]),
                    kinetics_mode="production_vbf"
                )
                if not sol.success:
                    out["F"] = [1e4, 1e4, 1e4]
                    out["G"] = [1e4, 1e4]
                    return

                res = evaluate_run_metrics(sol, geom, mode="cocurrent")
                f1 = -res['yield_khademi_eq28']
                f2 = -res['Benzene_rate_plant']
                f3 = res['T_max_exo'] - float(x[0])

                # Constraints: T_max <= 535 K, and P_endo_out >= 1.1 bar
                P_endo_out = (cfg.feed_shell.P0 / 1e5) - res['dP_shell']
                g1 = res['T_max_exo'] - cfg.optimization.T_max_limit
                g2 = 1.1 - P_endo_out

                out["F"] = [f1, f2, f3]
                out["G"] = [g1, g2]
            except Exception:
                out["F"] = [1e4, 1e4, 1e4]
                out["G"] = [1e4, 1e4]

    problem = TCROptimizationProblem()
    algorithm = NSGA2(
        pop_size=cfg.optimization.pop_size,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True
    )

    t0 = time.time()
    res_opt = pymoo_minimize(
        problem, algorithm, ('n_gen', cfg.optimization.n_gen),
        seed=cfg.optimization.seed, verbose=False
    )
    opt_time = time.time() - t0

    if res_opt.F is not None and len(res_opt.F) > 0:
        F_pareto = np.atleast_2d(res_opt.F)
        X_pareto = np.atleast_2d(res_opt.X)
    else:
        F_pareto = np.atleast_2d(res_opt.pop.get("F"))
        X_pareto = np.atleast_2d(res_opt.pop.get("X"))

    meoh_yields = -F_pareto[:, 0]
    benzene_rates = -F_pareto[:, 1]
    hotspot_rises = F_pareto[:, 2]

    # Select Knee Point using 3D Euclidean distance to utopia
    norm_y = (meoh_yields - np.min(meoh_yields)) / (np.max(meoh_yields) - np.min(meoh_yields) + 1e-10)
    norm_b = (benzene_rates - np.min(benzene_rates)) / (np.max(benzene_rates) - np.min(benzene_rates) + 1e-10)
    norm_h = (hotspot_rises - np.min(hotspot_rises)) / (np.max(hotspot_rises) - np.min(hotspot_rises) + 1e-10)
    dist_ideal = np.sqrt((1.0 - norm_y) ** 2 + (1.0 - norm_b) ** 2 + (norm_h) ** 2)
    knee_idx = int(np.argmin(dist_ideal))

    # High-resolution simulation at knee point
    knee_x = X_pareto[knee_idx]
    sol_knee = solve_cocurrent(
        geom, cfg.cat_tube, cfg.cat_shell, cfg.transport, cfg.solver,
        cfg.feed_tube.y0, float(knee_x[0]), cfg.feed_tube.P0, float(knee_x[3]),
        cfg.feed_shell.y0, float(knee_x[1]), cfg.feed_shell.P0, float(knee_x[2]),
        kinetics_mode="production_vbf"
    )
    res_knee = evaluate_run_metrics(sol_knee, geom, mode="cocurrent")
    res_knee['decision_variables'] = {
        'T_exo_0_K': float(knee_x[0]),
        'T_endo_0_K': float(knee_x[1]),
        'F_endo_mol_s': float(knee_x[2]),
        'F_exo_mol_s': float(knee_x[3])
    }

    return {
        'X_pareto': X_pareto,
        'F_pareto': F_pareto,
        'knee_idx': knee_idx,
        'opt_time_s': opt_time,
        'knee_results': res_knee
    }

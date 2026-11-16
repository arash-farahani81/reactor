"""
Configuration loader and typed dataclasses for TCR simulation.
Reads parameters from config.yaml with full source tracking.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import yaml
import os
import numpy as np
from .thermo import set_use_real_gas_eos

@dataclass
class Parameter:
    value: Any
    unit: str
    source: str

@dataclass
class GeometryConfig:
    L: float
    D_in: float
    D_out: float
    D_shell_eff: float
    N_tubes: int
    Kw: float
    p_t: Optional[float] = None
    # Derived quantities
    A_tube: float = field(init=False)
    A_shell: float = field(init=False)
    a_H: float = field(init=False)

    def __post_init__(self):
        if self.p_t is None:
            self.p_t = float(np.sqrt(np.pi / (2.0 * np.sqrt(3.0))) * self.D_shell_eff)
        self.A_tube = (np.pi / 4.0) * (self.D_in ** 2)
        self.A_shell = (np.pi / 4.0) * (self.D_shell_eff ** 2 - self.D_out ** 2)
        self.a_H = np.pi * self.D_in

@dataclass
class CatalystBedConfig:
    dp: float
    rho_b: float
    eps: float
    eta: float

@dataclass
class StreamFeedConfig:
    y0: np.ndarray
    T0: float
    P0: float
    F_total0: float
    components: Optional[List[str]] = None

@dataclass
class TransportConfig:
    h_cool: float
    T_coolant: float
    Kg_exo: float
    Kg_endo: float
    mu_exo: float
    mu_endo: float

@dataclass
class SolverConfig:
    z_points: int
    rtol: float
    atol: float
    picard_max_iter: int
    picard_tol_T: float
    picard_tol_F: float
    picard_alpha: float
    picard_tol_P: float = 1.0

@dataclass
class OptimizationConfig:
    seed: int
    pop_size: int
    n_gen: int
    bounds: Dict[str, List[float]]
    T_max_limit: float
    seeds: List[int] = field(default_factory=lambda: [42, 101, 777])
    pop_size_counter: int = 16
    n_gen_counter: int = 10
    P_shell_out_min: float = 1.10
    P_shell_out_knee: float = 1.15

@dataclass
class AppConfig:
    geometry: GeometryConfig
    cat_tube: CatalystBedConfig
    cat_shell: CatalystBedConfig
    feed_tube: StreamFeedConfig
    feed_shell: StreamFeedConfig
    feed_tube_bench: StreamFeedConfig
    feed_shell_bench: StreamFeedConfig
    transport: TransportConfig
    solver: SolverConfig
    optimization: OptimizationConfig
    raw_yaml: Dict[str, Any]
    use_real_gas_eos: bool = True
    use_dynamic_thiele: bool = True

def load_config(config_path: str = "config.yaml") -> AppConfig:
    """Loads configuration from YAML file into typed dataclasses."""
    if not os.path.exists(config_path):
        # Fallback to search in root directory
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        config_path = os.path.join(base_dir, "config.yaml")

    with open(config_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    geo_data = data["geometry"]
    p_t_val = float(geo_data["p_t"]["value"]) if "p_t" in geo_data and "value" in geo_data["p_t"] else None
    geom = GeometryConfig(
        L=float(geo_data["L"]["value"]),
        D_in=float(geo_data["D_in"]["value"]),
        D_out=float(geo_data["D_out"]["value"]),
        D_shell_eff=float(geo_data["D_shell_eff"]["value"]),
        N_tubes=int(geo_data["N_tubes"]["value"]),
        Kw=float(geo_data["Kw"]["value"]),
        p_t=p_t_val
    )

    cat_data = data["catalyst"]
    cat_tube = CatalystBedConfig(
        dp=float(cat_data["tube"]["dp"]["value"]),
        rho_b=float(cat_data["tube"]["rho_b"]["value"]),
        eps=float(cat_data["tube"]["eps"]["value"]),
        eta=float(cat_data["tube"]["eta"]["value"])
    )
    cat_shell = CatalystBedConfig(
        dp=float(cat_data["shell"]["dp"]["value"]),
        rho_b=float(cat_data["shell"]["rho_b"]["value"]),
        eps=float(cat_data["shell"]["eps"]["value"]),
        eta=float(cat_data["shell"]["eta"]["value"])
    )

    feed_prod = data["feed"]["production"]
    feed_tube = StreamFeedConfig(
        y0=np.array(feed_prod["tube"]["y0"], dtype=float),
        T0=float(feed_prod["tube"]["T0"]["value"]),
        P0=float(feed_prod["tube"]["P0"]["value"]),
        F_total0=float(feed_prod["tube"]["F_total0"]["value"]),
        components=feed_prod["tube"].get("components", ['H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2'])
    )
    feed_shell = StreamFeedConfig(
        y0=np.array(feed_prod["shell"]["y0"], dtype=float),
        T0=float(feed_prod["shell"]["T0"]["value"]),
        P0=float(feed_prod["shell"]["P0"]["value"]),
        F_total0=float(feed_prod["shell"]["F_total0"]["value"]),
        components=feed_prod["shell"].get("components", ['C6H12', 'C6H6', 'H2', 'Ar'][:len(feed_prod["shell"]["y0"])])
    )

    feed_bench = data["feed"]["benchmark_khademi"]
    feed_tube_bench = StreamFeedConfig(
        y0=np.array(feed_bench["tube"]["y0"], dtype=float),
        T0=float(feed_bench["tube"]["T0"]["value"]),
        P0=float(feed_bench["tube"]["P0"]["value"]),
        F_total0=float(feed_bench["tube"]["F_total0"]["value"]),
        components=feed_bench["tube"].get("components", ['H2', 'CO2', 'CO', 'CH3OH', 'H2O', 'CH4', 'N2'])
    )
    feed_shell_bench = StreamFeedConfig(
        y0=np.array(feed_bench["shell"]["y0"], dtype=float),
        T0=float(feed_bench["shell"]["T0"]["value"]),
        P0=float(feed_bench["shell"]["P0"]["value"]),
        F_total0=float(feed_bench["shell"]["F_total0"]["value"]),
        components=feed_bench["shell"].get("components", ['C6H12', 'C6H6', 'H2', 'Ar'])
    )

    trans_data = data["transport"]
    transport = TransportConfig(
        h_cool=float(trans_data["conventional"]["h_cool"]["value"]),
        T_coolant=float(trans_data["conventional"]["T_coolant"]["value"]),
        Kg_exo=float(trans_data["wall_heat_transfer"]["Kg_exo"]),
        Kg_endo=float(trans_data["wall_heat_transfer"]["Kg_endo"]),
        mu_exo=float(trans_data["wall_heat_transfer"]["mu_exo"]),
        mu_endo=float(trans_data["wall_heat_transfer"]["mu_endo"])
    )

    sol_data = data["solver"]
    solver = SolverConfig(
        z_points=int(sol_data["z_points"]),
        rtol=float(sol_data["rtol"]),
        atol=float(sol_data["atol"]),
        picard_max_iter=int(sol_data["picard_max_iter"]),
        picard_tol_T=float(sol_data["picard_tol_T"]),
        picard_tol_F=float(sol_data["picard_tol_F"]),
        picard_alpha=float(sol_data["picard_alpha"]),
        picard_tol_P=float(sol_data.get("picard_tol_P", 1.0))
    )

    use_real_gas_eos = bool(data.get("use_real_gas_eos", True))
    set_use_real_gas_eos(use_real_gas_eos)
    use_dynamic_thiele = bool(data.get("use_dynamic_thiele", True))

    opt_data = data["optimization"]
    optimization = OptimizationConfig(
        seed=int(opt_data["seed"]),
        pop_size=int(opt_data["pop_size"]),
        n_gen=int(opt_data["n_gen"]),
        bounds=opt_data["bounds"],
        T_max_limit=float(opt_data["T_max_limit"]),
        seeds=[int(s) for s in opt_data.get("seeds", [42, 101, 777])],
        pop_size_counter=int(opt_data.get("pop_size_counter", 16)),
        n_gen_counter=int(opt_data.get("n_gen_counter", 10)),
        P_shell_out_min=float(opt_data.get("P_shell_out_min", 1.10)),
        P_shell_out_knee=float(opt_data.get("P_shell_out_knee", 1.15))
    )

    return AppConfig(
        geometry=geom,
        cat_tube=cat_tube,
        cat_shell=cat_shell,
        feed_tube=feed_tube,
        feed_shell=feed_shell,
        feed_tube_bench=feed_tube_bench,
        feed_shell_bench=feed_shell_bench,
        transport=transport,
        solver=solver,
        optimization=optimization,
        raw_yaml=data,
        use_real_gas_eos=use_real_gas_eos,
        use_dynamic_thiele=use_dynamic_thiele
    )

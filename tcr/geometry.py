"""
Reactor geometry module: unified tube and shell dimensions with zero undocumented multipliers.
Based on industrial shell-and-tube specifications in Khademi et al. (2009), Table 6.
Supports tube and shell geometry for multi-component reactions (including 4-component shell mixtures).
"""

import numpy as np
from .config import GeometryConfig

class ReactorGeometry:
    def __init__(self, cfg: GeometryConfig):
        self.L = float(cfg.L)
        self.D_in = float(cfg.D_in)
        self.D_out = float(cfg.D_out)
        self.D_shell_eff = float(cfg.D_shell_eff)
        self.N_tubes = int(cfg.N_tubes)
        self.Kw = float(cfg.Kw)

        # Radii (m)
        self.R_in = float(self.D_in / 2.0)
        self.R_out = float(self.D_out / 2.0)

        # Single cross-sectional flow areas per tube (strictly self-consistent)
        self.A_tube = (np.pi / 4.0) * (self.D_in ** 2)

        # Bundle pitch and shell cross-sectional areas:
        # In triangular pitch bundle geometry: A_cell = (sqrt(3)/2) * p_pitch^2
        # Equivalent circular cell diameter: D_shell_eff = sqrt(2*sqrt(3)/pi) * p_pitch
        # Giving p_pitch = sqrt(pi / (2 * sqrt(3))) * D_shell_eff ~ 1.315 * D_out (within standard 1.25-1.33 TEMA range)
        self.p_pitch = float(np.sqrt(np.pi / (2.0 * np.sqrt(3.0))) * self.D_shell_eff)
        self.pitch_ratio = float(self.p_pitch / self.D_out)

        # Total shell bundle vessel cross-sectional area enclosed by N_tubes pitch cells:
        self.A_shell_bundle = float(self.N_tubes * (np.pi / 4.0) * (self.D_shell_eff ** 2))
        self.D_shell_vessel = float(np.sqrt(self.N_tubes) * self.D_shell_eff)

        # Net cross-sectional shell flow area across the entire bundle:
        # A_flow,shell = A_shell_bundle - N_tubes * (pi/4 * D_out^2)
        self.A_flow_shell = float(self.A_shell_bundle - self.N_tubes * (np.pi / 4.0) * (self.D_out ** 2))

        # Single consistent net shell flow area per tube
        self.A_shell = float(self.A_flow_shell / self.N_tubes)
        self.A_flow_shell_per_tube = self.A_shell

        # Inter-stream heat transfer perimeter based on tube inner diameter (m²/m)
        self.a_H = float(np.pi * self.D_in)

        # Hydraulic diameters (m):
        # 1. Tube inner hydraulic diameter:
        self.Dh_tube = float(self.D_in)

        # 2. Shell annular equivalent hydraulic diameter per tube:
        self.Dh_shell_annular = float(self.D_shell_eff - self.D_out)

        # 3. Shell bundle hydraulic diameter (4 * A_flow / P_wetted across the bundle):
        # Wetted perimeter includes tube outer walls and outer shell vessel wall:
        self.P_wetted_shell = float(self.N_tubes * np.pi * self.D_out + np.pi * self.D_shell_vessel)
        self.Dh_shell_bundle = float(4.0 * self.A_flow_shell / self.P_wetted_shell)

        # Default shell hydraulic diameter maintained for consistency with annular/bundle correlations:
        self.Dh_shell = float(self.Dh_shell_annular)

        # Total plant catalyst volumes (m³)
        self.V_tube_total = float(self.A_tube * self.L * self.N_tubes)
        self.V_shell_total = float(self.A_shell * self.L * self.N_tubes)

    @property
    def R_tube(self) -> float:
        """Inner tube radius in meters (R_in = D_in / 2)."""
        return self.R_in

    def calc_superficial_velocity(self, F_total: float, T: float, P: float, side: str = "tube") -> float:
        """
        Computes superficial gas velocity u_s (m/s) per tube from total molar flow F_total (mol/s)
        at temperature T (K) and pressure P (Pa):
          u_s = (F_total * R_GAS * T / P) / A_flow
        """
        from .thermo import R_GAS
        P_safe = max(float(P), 1.0)
        T_safe = max(float(T), 100.0)
        A = self.A_tube if side.lower() in ("tube", "exo") else self.A_shell
        return float((F_total * R_GAS * T_safe / P_safe) / A)

    def get_summary(self) -> dict:
        return {
            "L_m": self.L,
            "D_in_mm": self.D_in * 1000.0,
            "D_out_mm": self.D_out * 1000.0,
            "R_in_mm": self.R_in * 1000.0,
            "D_shell_eff_mm": self.D_shell_eff * 1000.0,
            "D_shell_vessel_m": self.D_shell_vessel,
            "p_pitch_mm": self.p_pitch * 1000.0,
            "pitch_ratio": self.pitch_ratio,
            "N_tubes": self.N_tubes,
            "A_tube_m2": self.A_tube,
            "A_shell_m2": self.A_shell,
            "A_flow_shell_total_m2": self.A_flow_shell,
            "a_H_m2_per_m": self.a_H,
            "Dh_tube_mm": self.Dh_tube * 1000.0,
            "Dh_shell_annular_mm": self.Dh_shell_annular * 1000.0,
            "Dh_shell_bundle_mm": self.Dh_shell_bundle * 1000.0,
            "V_tube_total_m3": self.V_tube_total,
            "V_shell_total_m3": self.V_shell_total
        }

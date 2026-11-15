"""
Reactor geometry module: unified tube and shell dimensions with zero undocumented multipliers.
Based on industrial shell-and-tube specifications in Khademi et al. (2009), Table 6.
"""

import numpy as np
from .config import GeometryConfig

class ReactorGeometry:
    def __init__(self, cfg: GeometryConfig):
        self.L = cfg.L
        self.D_in = cfg.D_in
        self.D_out = cfg.D_out
        self.D_shell_eff = cfg.D_shell_eff
        self.N_tubes = cfg.N_tubes
        self.Kw = cfg.Kw

        # Single cross-sectional flow areas per tube (strictly self-consistent)
        self.A_tube = (np.pi / 4.0) * (self.D_in ** 2)
        self.A_shell = (np.pi / 4.0) * (self.D_shell_eff ** 2 - self.D_out ** 2)
        
        # Inter-stream heat transfer perimeter based on tube inner diameter (m²/m)
        self.a_H = np.pi * self.D_in

        # Hydraulic diameters (m)
        self.Dh_tube = self.D_in
        self.Dh_shell = self.D_shell_eff - self.D_out

        # Total plant catalyst volumes (m³)
        self.V_tube_total = self.A_tube * self.L * self.N_tubes
        self.V_shell_total = self.A_shell * self.L * self.N_tubes

    def get_summary(self) -> dict:
        return {
            "L_m": self.L,
            "D_in_mm": self.D_in * 1000.0,
            "D_out_mm": self.D_out * 1000.0,
            "D_shell_eff_mm": self.D_shell_eff * 1000.0,
            "N_tubes": self.N_tubes,
            "A_tube_m2": self.A_tube,
            "A_shell_m2": self.A_shell,
            "a_H_m2_per_m": self.a_H,
            "V_tube_total_m3": self.V_tube_total,
            "V_shell_total_m3": self.V_shell_total
        }

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

        # Bundle pitch and shell cross-sectional areas via Triangular Tube Pitch Unit-Cell Method:
        # For unbaffled axial packed-bed shell reactors (Khademi et al., 2009 Fig. 1):
        # Unit-cell area: A_cell = (sqrt(3)/2) * p_t^2
        # Equivalent circular cell diameter: D_shell_eff = p_t * sqrt(2*sqrt(3)/pi)
        # Inverting: p_t = D_shell_eff * sqrt(pi / (2*sqrt(3))) ~ 0.05523 m (pt/D_out = 1.315, within TEMA 1.25-1.33)
        if getattr(cfg, 'p_t', None) is not None and cfg.p_t is not None:
            self.p_pitch = float(cfg.p_t)
        else:
            self.p_pitch = float(np.sqrt(np.pi / (2.0 * np.sqrt(3.0))) * self.D_shell_eff)
        self.p_t = self.p_pitch
        self.pitch_ratio = float(self.p_pitch / self.D_out)

        # Unit-cell cross-sectional area and net shell flow area per tube:
        # A_cell = (sqrt(3)/2) * p_t^2
        # A_shell = A_cell - (pi/4) * D_out^2 == (pi/4) * (D_shell_eff^2 - D_out^2)
        self.A_cell = float((np.sqrt(3.0) / 2.0) * (self.p_pitch ** 2))
        self.A_shell = float(self.A_cell - (np.pi / 4.0) * (self.D_out ** 2))
        self.A_flow_shell_per_tube = self.A_shell

        # Total shell bundle vessel cross-sectional area enclosed by N_tubes pitch cells:
        self.A_shell_bundle = float(self.N_tubes * self.A_cell)
        self.D_shell_vessel = float(np.sqrt(self.N_tubes) * self.D_shell_eff)
        self.A_flow_shell = float(self.A_shell * self.N_tubes)

        # Inter-stream heat transfer perimeter based on tube inner diameter (m²/m)
        self.a_H = float(np.pi * self.D_in)

        # Hydraulic diameters (m):
        # 1. Tube inner hydraulic diameter:
        self.Dh_tube = float(self.D_in)

        # 2. Shell unit-cell hydraulic diameter: Dh = 4 * A_shell / (pi * D_out)
        self.Dh_shell_unit_cell = float(4.0 * self.A_shell / (np.pi * self.D_out))

        # 3. Shell annular equivalent hydraulic diameter per tube:
        self.Dh_shell_annular = float(self.D_shell_eff - self.D_out)

        # 4. Shell bundle hydraulic diameter (4 * A_flow / P_wetted across the bundle):
        self.P_wetted_shell = float(self.N_tubes * np.pi * self.D_out + np.pi * self.D_shell_vessel)
        self.Dh_shell_bundle = float(4.0 * self.A_flow_shell / self.P_wetted_shell)

        # Default shell hydraulic diameter maintained for consistency:
        self.Dh_shell = float(self.Dh_shell_unit_cell)

        # Total plant catalyst volumes (m³)
        self.V_tube_total = float(self.A_tube * self.L * self.N_tubes)
        self.V_shell_total = float(self.A_shell * self.L * self.N_tubes)

    @property
    def R_tube(self) -> float:
        """Inner tube radius in meters (R_in = D_in / 2)."""
        return self.R_in

    def calc_superficial_velocity(self, F_total: float, T: float, P: float, side: str = "tube",
                                  Z: float = 1.0) -> float:
        """
        Computes superficial gas velocity u_s (m/s) per tube from total molar flow F_total (mol/s)
        at temperature T (K) and pressure P (Pa), accounting for compressibility Z:
          u_s = (Z * F_total * R_GAS * T / P) / A_flow
        """
        from .thermo import R_GAS
        P_safe = max(float(P), 1.0)
        T_safe = max(float(T), 100.0)
        Z_safe = max(float(Z), 0.1)
        A = self.A_tube if side.lower() in ("tube", "exo") else self.A_shell
        return float((Z_safe * F_total * R_GAS * T_safe / P_safe) / A)

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

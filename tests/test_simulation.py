"""Smoke tests for the corrected two-pass DOP853 workflow."""

import unittest
from dataclasses import replace

import numpy as np

from pasha_toomre.config import venus_config
from pasha_toomre.diagnostics import angular_momentum_conservation, diagnostics_for_tau
from pasha_toomre.simulation import build_report_tau_grid, run_simulation


class SimulationSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = replace(
            venus_config(),
            tau_final=0.5,
            probe_rtol=1.0e-9,
            probe_atol=1.0e-11,
            main_rtol=1.0e-10,
            main_atol=1.0e-12,
        )
        cls.result = run_simulation(cls.config)

    def test_detects_upward_crossing(self) -> None:
        crossing = self.result.evaluate_tau([0.0])[:, 0]
        self.assertLess(abs(crossing[0]), 1.0e-7)
        self.assertGreater(crossing[1], 0.0)
        self.assertLess(self.result.tau_start, -10.0)

    def test_report_grid_spans_central_encounter(self) -> None:
        grid = build_report_tau_grid(self.result)
        self.assertLess(grid[0], 0.0)
        self.assertGreaterEqual(grid[-1], 0.49)
        self.assertLess(np.min(np.abs(grid)), 5.0e-4)

    def test_dense_diagnostics_are_valid(self) -> None:
        tau = np.array([-1.0, 0.0, 0.2, 0.5])
        values = diagnostics_for_tau(self.result, tau)
        for name in ("Z", "zeta", "eta", "r", "rxy", "vr3d", "Lz"):
            self.assertTrue(np.all(np.isfinite(values[name])), name)

    def test_lz_is_conserved(self) -> None:
        check = angular_momentum_conservation(self.result)
        self.assertLess(check.absolute_delta_lz, 1.0e-8)


if __name__ == "__main__":
    unittest.main()

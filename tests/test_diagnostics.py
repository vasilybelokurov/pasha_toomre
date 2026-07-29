"""Tests for reconstructed orbital diagnostics."""

import unittest

import numpy as np

from pasha_toomre.config import venus_config
from pasha_toomre.diagnostics import (
    BASE_DIAGNOSTICS,
    EXTENDED_DIAGNOSTICS,
    compute_diagnostics,
)
from pasha_toomre.dynamics import initial_state


class DiagnosticTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = venus_config()
        state = initial_state(self.config)
        self.states = np.column_stack((state, state))
        self.tau = np.array([0.2, 1.0])
        self.values = compute_diagnostics(self.states, self.tau, self.config)

    def test_layout_sizes(self) -> None:
        self.assertEqual(len(BASE_DIAGNOSTICS), 12)
        self.assertEqual(len(EXTENDED_DIAGNOSTICS), 13)
        self.assertIn("eta", EXTENDED_DIAGNOSTICS)
        self.assertNotIn("eta", BASE_DIAGNOSTICS)

    def test_relative_vertical_coordinates(self) -> None:
        np.testing.assert_allclose(self.values["zeta"], 0.0)
        np.testing.assert_allclose(self.values["eta"], 2.0 * self.config.Z0)

    def test_initial_circular_osculating_orbit(self) -> None:
        np.testing.assert_allclose(self.values["e"], 0.0, atol=1.0e-8)
        np.testing.assert_allclose(self.values["ra"], self.config.r0, atol=1.0e-8)
        np.testing.assert_allclose(
            self.values["Lz"], np.sqrt(self.config.r0), atol=1.0e-12
        )

    def test_all_core_diagnostics_are_finite(self) -> None:
        for name, values in self.values.items():
            self.assertTrue(np.all(np.isfinite(values)), name)


if __name__ == "__main__":
    unittest.main()

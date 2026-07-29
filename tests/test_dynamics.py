"""Tests for initial conditions, equations, and event functions."""

import unittest

import numpy as np

from pasha_toomre.config import venus_config
from pasha_toomre.dynamics import (
    find_apocenter,
    find_pericenter,
    initial_state,
    make_rhs,
    sun_crosses_zero,
)


class DynamicsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = venus_config()
        self.state = initial_state(self.config)

    def test_documented_initial_state(self) -> None:
        self.assertEqual(self.state.shape, (8,))
        self.assertAlmostEqual(self.state[0], self.config.Z0)
        self.assertAlmostEqual(self.state[2], self.config.r0)
        self.assertAlmostEqual(self.state[4] - self.state[0], 0.0)
        self.assertAlmostEqual(self.state[7] - self.state[1], 0.0)
        self.assertAlmostEqual(self.state[6], np.sqrt(1.0 / self.config.r0))

    def test_rhs_is_finite_eight_vector(self) -> None:
        derivative = make_rhs(self.config)(0.0, self.state)
        self.assertEqual(derivative.shape, (8,))
        self.assertTrue(np.all(np.isfinite(derivative)))
        self.assertGreater(derivative[1], 0.0)

    def test_crossing_event_is_scalar_and_terminal(self) -> None:
        value = sun_crosses_zero(0.0, self.state)
        self.assertIsInstance(value, float)
        self.assertEqual(value, self.config.Z0)
        self.assertTrue(sun_crosses_zero.terminal)
        self.assertEqual(sun_crosses_zero.direction, 1.0)

    def test_extremum_event_directions(self) -> None:
        self.assertEqual(find_pericenter.direction, 1.0)
        self.assertEqual(find_apocenter.direction, -1.0)


if __name__ == "__main__":
    unittest.main()

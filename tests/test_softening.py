"""Tests for the documented Sun-intruder softening law."""

import unittest

import numpy as np

from pasha_toomre.softening import get_softening, plummer4


class Plummer4Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.R = 0.00464913034

    def test_is_finite_and_zero_at_origin(self) -> None:
        value = plummer4(0.0, self.R)
        self.assertTrue(np.isfinite(value))
        self.assertEqual(value, 0.0)

    def test_is_odd_and_attractive(self) -> None:
        coordinates = np.array([0.01, 0.1, 1.0])
        positive = plummer4(coordinates, self.R)
        negative = plummer4(-coordinates, self.R)
        np.testing.assert_allclose(negative, -positive, rtol=1.0e-14, atol=0.0)
        self.assertTrue(np.all(positive < 0.0))
        self.assertTrue(np.all(negative > 0.0))

    def test_has_expected_far_field_limit(self) -> None:
        coordinate = 1.0e4 * self.R
        expected = -1.0 / (4.0 * coordinate**2)
        self.assertAlmostEqual(plummer4(coordinate, self.R) / expected, 1.0, places=7)

    def test_registry_returns_function(self) -> None:
        self.assertIs(get_softening("plummer4"), plummer4)


if __name__ == "__main__":
    unittest.main()

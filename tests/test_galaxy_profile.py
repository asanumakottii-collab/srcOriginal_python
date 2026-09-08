import math
import unittest

import numpy as np

from galaxy_profile import exponential_scale_arcsec


class ExponentialGalaxyProfileTests(unittest.TestCase):
    def test_circular_half_light_radius_matches_analytic_curve_of_growth(self):
        h = exponential_scale_arcsec(300.0, 1.0, 10.0, 25.0)
        t = 25.0 / h
        self.assertAlmostEqual(0.5, 1.0 - (1.0 + t) * math.exp(-t), places=12)

    def test_elliptical_model_contains_half_the_light_in_circular_aperture(self):
        # 生成側の楕円座標積分とは独立に、空の位置角phiで光量を積分する。
        phi = np.arange(32768) * (2.0 * math.pi / 32768)
        for q in (1.2, 4.0, 28.2):
            with self.subTest(axis_ratio=q):
                h = exponential_scale_arcsec(300.0, q, 10.0, 25.0)
                k2 = np.cos(phi) ** 2 + (q * np.sin(phi)) ** 2
                t = 25.0 * np.sqrt(k2) / h
                fraction = np.mean(q / k2 * (1.0 - (1.0 + t) * np.exp(-t)))
                self.assertAlmostEqual(0.5, fraction, places=8)

    def test_half_light_aperture_takes_priority_over_isophote_and_magnitude(self):
        h = exponential_scale_arcsec(300.0, 4.0, 10.0, 25.0)
        self.assertEqual(h, exponential_scale_arcsec(None, 4.0, None, 25.0))
        self.assertEqual(h, exponential_scale_arcsec(900.0, 4.0, 18.0, 25.0))

    def test_isophote_fit_recovers_known_exponential_disk_on_concentrated_branch(self):
        # I25=1に規格化した既知の円盤を観測し、hを復元する。
        for t in (2.0001, 3.0, 8.0, 40.0):
            with self.subTest(a25_in_scales=t):
                h_expected, q = 12.0, 3.0
                a25 = t * h_expected
                flux = 2.0 * math.pi / q * h_expected ** 2 * math.exp(t)
                bmag = 25.0 - 2.5 * math.log10(flux)
                h = exponential_scale_arcsec(a25, q, bmag)
                self.assertAlmostEqual(h_expected, h, places=7)
                self.assertLessEqual(h, a25 / 2.0)

    def test_isophote_maximum_and_no_solution_fallback(self):
        a25, q = 60.0, 2.0
        bmag_at_max = 25.0 + 2.5 * math.log10(4.0 / math.e ** 2 * q
                                             / (2.0 * math.pi * a25 ** 2))
        self.assertEqual(30.0, exponential_scale_arcsec(a25, q, bmag_at_max))
        h = exponential_scale_arcsec(a25, q, bmag_at_max - 1e-10)
        self.assertAlmostEqual(30.0, h, delta=0.001)
        self.assertEqual(20.0, exponential_scale_arcsec(a25, q, bmag_at_max + 0.01))

    def test_missing_information_uses_documented_fallback_or_skips(self):
        self.assertEqual(20.0, exponential_scale_arcsec(60.0, 2.0, None))
        self.assertIsNone(exponential_scale_arcsec(None, 2.0, 10.0))
        self.assertEqual(exponential_scale_arcsec(60.0, 2.0, 10.0),
                         exponential_scale_arcsec(60.0, 2.0, 10.0, float('nan')))

    def test_invalid_axis_ratio_is_rejected(self):
        for q in (0.0, 0.5, float('inf'), float('nan')):
            with self.subTest(axis_ratio=q):
                with self.assertRaises(ValueError):
                    exponential_scale_arcsec(60.0, q, 10.0)


if __name__ == '__main__':
    unittest.main()

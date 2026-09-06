import importlib.util
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

from galaxy_transformer import GalaxyTransformer, _init_plate_writer, main
from models import SphereConstellation, SpherePosition, SphereStar
from plate_writer import PlateWriterType


def _position(longitude, latitude=0.):
    return SpherePosition.from_galactic(longitude, latitude)


def _constellation(first, last):
    constellation = SphereConstellation()
    constellation.ll = [[_position(*first), _position(*last)]]
    return constellation


class GalacticCoordinatesTests(unittest.TestCase):
    def test_canonical_galactic_pole_and_equator_node(self):
        # Hipparcosの定義角。行列の逆変換を使わずにICRS座標から検証する。
        position = SpherePosition()
        position.radeg, position.dedeg = 192.85948, 27.12825
        self.assertAlmostEqual(90., position.to_galactic()[1], places=10)
        position.radeg, position.dedeg = 282.85948, 0.
        longitude, latitude = position.to_galactic()
        self.assertAlmostEqual(32.93192, longitude, places=10)
        self.assertAlmostEqual(0., latitude, places=10)

        equatorial = _position(32.93192)
        self.assertAlmostEqual(282.85948, equatorial.radeg % 360., places=10)
        self.assertAlmostEqual(0., equatorial.dedeg, places=10)

    def test_galactic_coordinates_round_trip(self):
        for longitude in (-360., -120., 0., 60., 180., 359., 720.):
            for latitude in (-80., -20., 0., 20., 80.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    actual_l, actual_b = _position(longitude, latitude).to_galactic()
                    self.assertAlmostEqual(0., (actual_l - longitude + 180.) % 360. - 180., places=10)
                    self.assertAlmostEqual(latitude, actual_b, places=10)


class GalaxyTransformerTests(unittest.TestCase):
    def setUp(self):
        self.transformer = GalaxyTransformer(0.125, 6500., 355., 200., 50.)

    def test_each_region_contains_both_galactic_latitude_signs(self):
        for longitude, unit in ((30., (0, 0)), (90., (0, 1)),
                                (330., (1, 0)), (210., (1, 1))):
            for latitude in (-19., 0., 19.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    position = _position(longitude, latitude)
                    self.assertEqual(unit, self.transformer.assigned_unit(position))
                    transformed = self.transformer.transform(position)
                    self.assertEqual(unit, (transformed.dir, transformed.index))
                    position.radeg += 360.
                    self.assertEqual(unit, self.transformer.assigned_unit(position))

    def test_longitude_boundaries_have_one_owner_and_gaps_are_excluded(self):
        for longitude, unit in (
            (0., (0, 0)), (60. - 1e-6, (0, 0)), (60., (0, 1)),
            (120. - 1e-6, (0, 1)), (120., None), (150., None),
            (180. - 1e-6, None), (180., (1, 1)),
            (240. - 1e-6, (1, 1)), (240., None), (270., None),
            (300. - 1e-6, None), (300., (1, 0)), (300. + 1e-6, (1, 0)),
            (360. - 1e-6, (1, 0)), (360., (0, 0)), (-60., (1, 0)),
        ):
            with self.subTest(longitude=longitude):
                self.assertEqual(unit, self.transformer.assigned_unit(_position(longitude)))

    def test_latitude_limits_are_inclusive_and_outside_is_excluded(self):
        for longitude in (30., 90., 210., 300.):
            for latitude in (-20., 20.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNotNone(self.transformer.transform(_position(longitude, latitude)))
            for latitude in (-90., -20.000001, 20.000001, 90.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNone(self.transformer.transform(_position(longitude, latitude)))

    def test_each_region_center_projects_to_plate_center_with_projector_offset(self):
        for longitude in (30., 90., 330., 210.):
            with self.subTest(longitude=longitude):
                position = self.transformer.transform(_position(longitude))
                self.assertAlmostEqual(0., position.xmm, places=10)
                self.assertAlmostEqual(0., position.ymm, places=10)

    def test_cylindrical_projection_scale_and_direction(self):
        transformer = GalaxyTransformer(0.125, 6500., 0., 0., 50.)
        for center in (30., 90., 330., 210.):
            for delta_l, latitude in ((-10., -15.), (10., 15.)):
                with self.subTest(center=center, delta_l=delta_l):
                    position = transformer.transform(_position(center + delta_l, latitude))
                    self.assertAlmostEqual(50. * math.tan(math.radians(latitude)), position.xmm, places=10)
                    self.assertAlmostEqual(-50. * math.radians(delta_l), position.ymm, places=10)

    def test_entire_band_fits_default_frames(self):
        for longitude in range(360):
            for latitude in (-20., 0., 20.):
                position = self.transformer.transform(_position(longitude, latitude))
                if position is not None:
                    with self.subTest(longitude=longitude, latitude=latitude):
                        self.assertTrue(math.isfinite(position.xmm))
                        self.assertTrue(math.isfinite(position.ymm))
                        self.assertLess(abs(position.xmm), 69.25)
                        self.assertLess(abs(position.ymm), 69.25)

    def test_process_stars_skips_outside_band_and_longitude_gaps(self):
        stars = []
        for longitude, latitude in ((30., 0.), (30., 21.), (150., 0.), (270., 0.), (300., -10.)):
            star = SphereStar()
            star.p, star.vmag = _position(longitude, latitude), 7.5
            stars.append(star)
        reader = MagicMock()
        reader.read_star.side_effect = [*stars, None]
        writer = MagicMock()
        with patch("sys.stdout", new_callable=StringIO):
            self.transformer.process_stars(reader, writer)
        positions = [call.args[0].p for call in writer.write_star.call_args_list]
        self.assertEqual([(0, 0), (1, 0)], [(p.dir, p.index) for p in positions])

    def _assert_line(self, actual, first, last, unit):
        for position, expected in zip(actual, (first, last)):
            target = self.transformer.transform_unit(_position(*expected), *unit)
            self.assertEqual(unit, (position.dir, position.index))
            self.assertAlmostEqual(target.xmm, position.xmm, places=9)
            self.assertAlmostEqual(target.ymm, position.ymm, places=9)

    def test_constellation_line_is_split_at_longitude_boundary(self):
        result = self.transformer.transform_constellation(_constellation((40., 0.), (80., 0.)))
        self.assertEqual(2, len(result.ll))
        self._assert_line(result.ll[0], (40., 0.), (60., 0.), (0, 0))
        self._assert_line(result.ll[1], (60., 0.), (80., 0.), (0, 1))

    def test_constellation_line_is_split_at_zero_longitude(self):
        result = self.transformer.transform_constellation(_constellation((350., 0.), (10., 0.)))
        self.assertEqual(2, len(result.ll))
        self._assert_line(result.ll[0], (0., 0.), (10., 0.), (0, 0))
        self._assert_line(result.ll[1], (350., 0.), (360., 0.), (1, 0))

    def test_constellation_line_crossing_band_is_clipped_even_with_both_ends_outside(self):
        result = self.transformer.transform_constellation(_constellation((30., 30.), (30., -30.)))
        self.assertEqual(1, len(result.ll))
        self._assert_line(result.ll[0], (30., 20.), (30., -20.), (0, 0))

    def test_constellation_line_crossing_longitude_gaps_is_clipped(self):
        for first, last, expected in (
            (110., 190., ((110., 120., (0, 1)), (180., 190., (1, 1)))),
            (230., 310., ((300., 310., (1, 0)), (230., 240., (1, 1)))),
        ):
            with self.subTest(first=first, last=last):
                result = self.transformer.transform_constellation(_constellation((first, 0.), (last, 0.)))
                self.assertEqual(2, len(result.ll))
                for line, (start, end, unit) in zip(result.ll, expected):
                    self._assert_line(line, (start, 0.), (end, 0.), unit)

    def test_constellation_on_shared_boundary_is_not_duplicated(self):
        result = self.transformer.transform_constellation(_constellation((60., -10.), (60., 10.)))
        self.assertEqual(1, len(result.ll))
        self._assert_line(result.ll[0], (60., -10.), (60., 10.), (0, 1))

    def test_constellation_on_latitude_boundary_is_retained(self):
        for latitude in (-20., 20.):
            result = self.transformer.transform_constellation(_constellation((10., latitude), (50., latitude)))
            self.assertIsNotNone(result)
            self.assertEqual(1, len(result.ll))

    def test_constellation_name_and_lines_outside_are_excluded(self):
        constellation = _constellation((30., 30.), (40., 30.))
        constellation.name, constellation.p = "Outside", _position(30., 40.)
        self.assertIsNone(self.transformer.transform_constellation(constellation))
        constellation.ll = [[_position(20.), _position(40.)]]
        result = self.transformer.transform_constellation(constellation)
        self.assertIsNone(result.p)
        self.assertEqual(1, len(result.ll))


class GalaxyOutputTests(unittest.TestCase):
    def test_svg_creates_all_four_frames_without_stars(self):
        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            writer = _init_plate_writer({"output.directory": directory}, PlateWriterType.SVG)
            writer.close()
            pages = sorted(Path(directory, "galaxy").glob("*.svg"))
            self.assertEqual(["galaxy-0.svg", "galaxy-1.svg"], [p.name for p in pages])
            for page, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                root = ET.parse(page).getroot()
                self.assertEqual(2, len(root.findall("{http://www.w3.org/2000/svg}rect")))
                self.assertEqual(labels, [e.text for e in root.findall("{http://www.w3.org/2000/svg}text")])

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_pdf_creates_all_four_frames_without_stars(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            writer = _init_plate_writer({"output.directory": directory}, PlateWriterType.PDF)
            writer.close()
            pages = sorted(Path(directory, "galaxy").glob("*.pdf"))
            self.assertEqual(2, len(pages))
            for path, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                pdf = PdfReader(path)
                self.assertEqual(1, len(pdf.pages))
                self.assertEqual(labels, pdf.pages[0].extract_text().split())

    def test_command_filters_catalog_and_clips_constellations(self):
        stars = []
        for longitude, latitude in ((30., 0.), (90., 0.), (300., 0.), (210., 0.),
                                    (150., 0.), (270., 0.), (30., 21.)):
            star = SphereStar()
            star.p, star.vmag = _position(longitude, latitude), 7.5
            stars.append(star)
        reader = MagicMock()
        reader.read_star.side_effect = [*stars, None]
        reader.read_constellation.side_effect = [_constellation((30., 30.), (30., -30.)), None]
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory, "galaxy.properties")
            config.write_text(f"output.directory = {directory}\n", encoding="utf-8")
            with (patch("galaxy_transformer._init_sphere_reader", return_value=reader),
                  patch("builtins.input", return_value="y"),
                  patch("sys.stdout", new_callable=StringIO)):
                main(["-f", str(config)])
            pages = [ET.parse(p).getroot() for p in sorted(Path(directory, "galaxy").glob("*.svg"))]
            self.assertEqual(2, len(pages))
            self.assertEqual([2, 2], [len(p.findall("{http://www.w3.org/2000/svg}circle")) for p in pages])
            self.assertEqual([1, 0], [len(p.findall("{http://www.w3.org/2000/svg}line")) for p in pages])


if __name__ == "__main__":
    unittest.main()

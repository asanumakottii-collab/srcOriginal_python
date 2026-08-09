import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from basic_transformer import BasicTransformer
from mathvector import MathVector
from models import PlatePosition, SpherePosition, SphereStar
from plate_polygon import PlateAssignmentPolygonGenerator
from plate_writer import PlateWriterSVG
from sphere_reader import SphereReader
from transformer import Transformer, _init_enlarge_rate
from unit_arrangement import CustomUnitArrangement


class _IdentityTransformer(BasicTransformer):
    def __init__(self):
        super().__init__()
        self.radius = 0.25

    def transform(self, sp):
        p = PlatePosition()
        p.xmm = sp.radeg
        p.ymm = sp.dedeg
        return p

    def transform_unit(self, sp, dir_, index):
        return self.transform(sp)


class _StarReader:
    def __init__(self, stars, constellation_ids):
        self.stars = iter(stars)
        self.constellation_ids = constellation_ids

    def read_star(self):
        return next(self.stars, None)

    def constellation_star_hip_numbers(self):
        return self.constellation_ids


class _StarWriter:
    def __init__(self):
        self.stars = []

    def write_star(self, star):
        self.stars.append(star)


def _star(hip_number, magnitude=4.0):
    star = SphereStar()
    star.hip_number = hip_number
    star.vmag = magnitude
    star.p = SpherePosition()
    return star


class ConstellationStarTests(unittest.TestCase):
    def test_legacy_and_hlc_hip_numbers_are_combined(self):
        hip_numbers = SphereReader.constellation_star_hip_numbers()
        self.assertEqual(893, len(hip_numbers))
        self.assertIn(677, hip_numbers)
        self.assertIn(2072, hip_numbers)  # 補助CSVだけに存在

    def test_enlarge_formula_matches_fourth_magnitude_rule(self):
        base_radius = 0.25
        rmm = base_radius * math.pow(math.sqrt(2.512), -4)
        actual = BasicTransformer.enlarge_constellation_star(rmm, base_radius, 1.0)
        self.assertAlmostEqual(rmm * math.pow(100, 0.1), actual)
        self.assertEqual(rmm, BasicTransformer.enlarge_constellation_star(rmm, base_radius, 0.0))

    def test_only_constellation_stars_are_enlarged(self):
        transformer = _IdentityTransformer()
        writer = _StarWriter()
        transformer.process_stars(
            _StarReader([_star(10), _star(20)], {20}), writer, enlarge_rate=1.0
        )
        self.assertEqual(2, len(writer.stars))
        self.assertAlmostEqual(0.25 * math.pow(math.sqrt(2.512), -4), writer.stars[0].rmm)
        self.assertGreater(writer.stars[1].rmm, writer.stars[0].rmm)

    def test_legacy_enlarge_property_name_is_supported(self):
        self.assertEqual(1.5, _init_enlarge_rate({"star-EnlargeRate": "1.5"}))


class AssignmentPolygonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        north = [
            MathVector.from_mag_lng_lat(1, 0, math.radians(45)),
            MathVector.from_mag_lng_lat(1, math.pi, math.radians(45)),
        ]
        south = [
            MathVector.from_mag_lng_lat(1, 0, math.radians(-45)),
            MathVector.from_mag_lng_lat(1, math.pi, math.radians(-45)),
        ]
        cls.transformer = Transformer(
            0.25, 6000, 300, 200, 50, CustomUnitArrangement(north, south)
        )

    def test_inverse_transform_round_trip(self):
        for dir_ in range(2):
            for index in range(2):
                for xmm, ymm in ((0, 0), (5, 10), (-12, 8)):
                    sphere_position = self.transformer.inverse_transform_unit(
                        dir_, index, xmm, ymm
                    )
                    plate_position = self.transformer.transform_unit(
                        sphere_position, dir_, index
                    )
                    self.assertAlmostEqual(xmm, plate_position.xmm, places=10)
                    self.assertAlmostEqual(ymm, plate_position.ymm, places=10)

    def test_polygon_points_stay_inside_frame_and_belong_to_plate(self):
        generator = PlateAssignmentPolygonGenerator(
            self.transformer, samples=24, iterations=18
        )
        points = generator.generate(0, 0, 37.5)
        self.assertEqual(24, len(points))
        for x, y in points:
            self.assertLessEqual(math.hypot(x, y), 37.5 + 1e-12)
            # 二分探索点の丸め誤差を避け、境界の少し内側を検証する。
            position = self.transformer.inverse_transform_unit(0, 0, x * 0.999, y * 0.999)
            self.assertEqual((0, 0), self.transformer.assigned_unit(position))

    def test_svg_writer_marks_polygon_with_plate_name(self):
        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterSVG(1, 1, 37.5, True, "polygon-", True, directory)
            writer.write_assignment_polygon(0, 0, [(0, 0), (10, 0), (0, 10)])
            writer.close()
            output = Path(directory, "polygon-0.svg")
            root = ET.parse(output).getroot()
            polygons = root.findall("{http://www.w3.org/2000/svg}polygon")
            self.assertEqual("N0", polygons[0].attrib["data-plate"])
            self.assertEqual("white", polygons[0].attrib["fill"])
            first_x, first_y = map(float, polygons[0].attrib["points"].split()[0].split(","))
            self.assertAlmostEqual(95.5 * 96 / 25.4, first_x)
            self.assertAlmostEqual(69.25 * 96 / 25.4, first_y)


if __name__ == "__main__":
    unittest.main()

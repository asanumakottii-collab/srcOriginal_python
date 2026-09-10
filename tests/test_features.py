import importlib.util
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

from basic_transformer import BasicTransformer
from mathvector import MathVector
from models import PlateConstellation, PlatePosition, PlateStar, SpherePosition, SphereStar
from plate_polygon import PlateAssignmentPolygonGenerator
from plate_writer import PlateWriterPDF, PlateWriterSVG, PlateWriterType, categorized_output_dir
from sphere_reader import SphereReader
from transformer import (
    _PlateWriterGroup,
    Transformer,
    _init_enlarge_constellation_names,
    _init_enlarge_rate,
    _init_plate_writer,
    _init_sphere_reader,
    _select_interactive_output_types,
    _parse_constellation_names,
    _write_assignment_polygons,
)
from unit_arrangement import CustomUnitArrangement


_HAS_REPORTLAB = importlib.util.find_spec("reportlab") is not None
_HAS_PYPDF = importlib.util.find_spec("pypdf") is not None


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


class _NamedStarReader(_StarReader):
    def __init__(self, stars, constellation_ids):
        super().__init__(stars, constellation_ids)
        self.requested_names = None

    def constellation_star_hip_numbers(self, constellation_names=None):
        self.requested_names = constellation_names
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

    def test_constellation_names_are_loaded_from_hlc(self):
        names = SphereReader.constellation_names()
        self.assertEqual(89, len(names))
        self.assertEqual("Andromeda", names[0])
        self.assertIn("Serpens (Head)", names)

    def test_only_selected_constellation_hip_numbers_are_returned(self):
        hip_numbers = SphereReader.constellation_star_hip_numbers(
            ["Andromeda"]
        )
        self.assertIn(677, hip_numbers)
        self.assertIn(7607, hip_numbers)
        self.assertNotIn(53502, hip_numbers)  # Antlia
        self.assertNotIn(2072, hip_numbers)  # 星座名のない補助CSVの星

    def test_multiple_selected_constellations_are_combined(self):
        andromeda = SphereReader.constellation_star_hip_numbers(["Andromeda"])
        antlia = SphereReader.constellation_star_hip_numbers(["Antlia"])
        combined = SphereReader.constellation_star_hip_numbers(
            ["Andromeda", "Antlia", "Andromeda"]
        )
        self.assertEqual(andromeda | antlia, combined)

    def test_unknown_constellation_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown"):
            SphereReader.constellation_star_hip_numbers(["Unknown"])

    def test_enlarge_formula_matches_fourth_magnitude_rule(self):
        base_radius = 0.25
        rmm = base_radius * math.pow(math.sqrt(2.512), -4)
        actual = BasicTransformer.enlarge_constellation_star(
            rmm, base_radius, 1.0, magnitude=4.0
        )
        self.assertAlmostEqual(rmm * math.pow(10.0, 0.2), actual)
        self.assertEqual(
            rmm,
            BasicTransformer.enlarge_constellation_star(rmm, base_radius, 0.0),
        )

    def test_bright_stars_use_fixed_width_gaussian_taper(self):
        rmm = 0.25 * math.pow(math.sqrt(2.512), -2)
        taper = math.exp(-math.pow((4.0 - 2.0) / 2.0, 2))
        expected = rmm * math.pow(10.0, 0.2 * 1.5 * taper)
        actual = BasicTransformer.enlarge_constellation_star(
            rmm, 0.25, 1.5, magnitude=2.0
        )
        self.assertAlmostEqual(expected, actual)

    def test_faint_stars_receive_the_full_magnitude_shift(self):
        rmm = 0.25 * math.pow(math.sqrt(2.512), -5)
        actual = BasicTransformer.enlarge_constellation_star(
            rmm, 0.25, 2.0, magnitude=5.0
        )
        self.assertAlmostEqual(rmm * math.pow(10.0, 0.4), actual)

    def test_enlarge_formula_keeps_legacy_call_signature(self):
        base_radius = 0.25
        rmm = base_radius * math.pow(10.0, -0.2 * 2.0)
        explicit = BasicTransformer.enlarge_constellation_star(
            rmm, base_radius, 1.5, magnitude=2.0
        )
        inferred = BasicTransformer.enlarge_constellation_star(rmm, base_radius, 1.5)
        self.assertAlmostEqual(explicit, inferred)

    def test_enlarged_radius_preserves_magnitude_order_at_maximum_rate(self):
        enlarged_radii = []
        for magnitude in [value / 10 for value in range(-10, 81)]:
            rmm = 0.25 * math.pow(math.sqrt(2.512), -magnitude)
            enlarged_radii.append(
                BasicTransformer.enlarge_constellation_star(
                    rmm, 0.25, 2.0, magnitude=magnitude
                )
            )
        self.assertTrue(all(a > b for a, b in zip(enlarged_radii, enlarged_radii[1:])))

    def test_only_constellation_stars_are_enlarged(self):
        transformer = _IdentityTransformer()
        writer = _StarWriter()
        transformer.process_stars(
            _StarReader([_star(10), _star(20)], {20}), writer, enlarge_rate=1.0
        )
        self.assertEqual(2, len(writer.stars))
        self.assertAlmostEqual(0.25 * math.pow(math.sqrt(2.512), -4), writer.stars[0].rmm)
        self.assertAlmostEqual(
            0.25 * math.pow(math.sqrt(2.512), -4) * math.pow(10.0, 0.2),
            writer.stars[1].rmm,
        )

    def test_selected_constellation_names_are_passed_to_reader(self):
        transformer = _IdentityTransformer()
        writer = _StarWriter()
        reader = _NamedStarReader([_star(10), _star(20)], {20})
        transformer.process_stars(
            reader,
            writer,
            enlarge_rate=1.0,
            constellation_names=("Orion",),
        )
        self.assertEqual(("Orion",), reader.requested_names)
        self.assertAlmostEqual(
            0.25 * math.pow(math.sqrt(2.512), -4) * math.pow(10.0, 0.2),
            writer.stars[1].rmm,
        )

    def test_constellation_name_list_accepts_commas_newlines_and_comments(self):
        self.assertEqual(
            ("Orion", "Canis Major", "Ursa Major"),
            _parse_constellation_names(
                "# selected constellations\nOrion, Canis Major\n\nUrsa Major\n"
            ),
        )

    def test_constellation_names_are_loaded_from_properties(self):
        names = _init_enlarge_constellation_names(
            {"star.enlarge-constellations": "Orion, Canis Major"}
        )
        self.assertEqual(("Orion", "Canis Major"), names)

    def test_constellation_names_are_loaded_from_interactive_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "constellations.txt"
            path.write_text("Orion\nCanis Major\n", encoding="utf-8")
            with patch("builtins.input", return_value=str(path)):
                names = _init_enlarge_constellation_names(None)
        self.assertEqual(("Orion", "Canis Major"), names)

    def test_legacy_enlarge_property_name_is_supported(self):
        self.assertEqual(1.5, _init_enlarge_rate({"star-EnlargeRate": "1.5"}))

    def test_enlarge_rate_must_be_between_zero_and_two(self):
        self.assertEqual(0.0, _init_enlarge_rate({"star.enlarge-rate": "0.0"}))
        self.assertEqual(2.0, _init_enlarge_rate({"star.enlarge-rate": "2.0"}))
        for value in ("-0.1", "2.1", "nan"):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    _init_enlarge_rate({"star.enlarge-rate": value})


class SphereReaderTests(unittest.TestCase):
    @staticmethod
    def _rc3_record(btmag, bv, log_d25, bt_code=" ", log_ae=None):
        record = [" "] * 256
        record[0:2] = "12"
        record[2:4] = "00"
        record[4:8] = "00.0"
        record[9] = "+"
        record[10:12] = "20"
        record[12:14] = "00"
        record[14:16] = "00"
        if log_d25 is not None:
            record[151:155] = f"{log_d25:4.2f}"
        record[157:160] = ".00"
        record[161:165] = "0.00"
        if log_ae is not None:
            record[176:180] = f"{log_ae:4.2f}"
        record[185:188] = "  0"
        record[189:194] = f"{btmag:5.2f}"
        record[194] = bt_code
        if bv is not None:
            record[252:256] = f"{bv:4.2f}"
        return "".join(record) + "\n"

    def test_rc3_axis_ratio_uses_log_r25_not_diameter_error(self):
        for diameter_error, log_r25 in ((".00", "0.60"), (".12", "0.60"),
                                         ("   ", "0.60"), (".12", "    ")):
            with self.subTest(diameter_error=diameter_error, log_r25=log_r25):
                record = list(self._rc3_record(9.0, 1.0, 1.0))
                record[157:160] = diameter_error
                record[161:165] = log_r25
                streams = [StringIO(), StringIO(), StringIO("".join(record)), StringIO()]
                with patch("builtins.open", side_effect=streams):
                    reader = SphereReader(False, 1.5, 10.0, False, rc3_enabled=True)
                try:
                    star = reader.read_star()
                    if not log_r25.strip():
                        self.assertIsNone(star)
                        continue
                    self.assertIsNotNone(star)
                    major, minor = reader._galaxy_vectors[:2]
                    self.assertAlmostEqual(10 ** 0.60, major.get_mag() / minor.get_mag())
                finally:
                    reader.close()

    def test_catalog_magnitude_is_preserved_within_limits(self):
        tycho_record = [" "] * 216
        tycho_record[41:46] = f"{4.5:5.2f}"
        tycho_record[51:63] = f"{12.0:12.6f}"
        tycho_record[64:76] = f"{34.0:12.6f}"
        streams = [
            StringIO(),
            StringIO("".join(tycho_record) + "\n"),
            StringIO(),
            StringIO(),
        ]

        with patch("builtins.open", side_effect=streams):
            reader = SphereReader(False, 1.5, 7.5, False)
            star = reader.read_star()

        self.assertEqual(4.5, star.vmag)

    def test_inconsistent_isophote_uses_fallback_and_continues(self):
        invalid_record = self._rc3_record(10.0, 0.5, 3.0)
        valid_record = self._rc3_record(10.0, 0.5, 1.0)
        streams = [
            StringIO(),
            StringIO(),
            StringIO(invalid_record + valid_record),
            StringIO(),
        ]

        with patch("builtins.open", side_effect=streams):
            reader = SphereReader(False, 7.5, 10.0, True, rc3_enabled=True)
        try:
            star = reader.read_star()
            self.assertIsNotNone(star)
            # D25=6000秒角なので、フォールバックのh=a25/3は1000秒角。
            self.assertAlmostEqual(1000.0, math.degrees(
                reader._galaxy_vectors[0].get_mag()) * 3600.0)
            self.assertIsNotNone(reader.read_star())
            self.assertIsNone(reader.read_star())
        finally:
            reader.close()

    def test_rc3_ae_is_a_circular_diameter_and_d25_is_optional(self):
        for log_d25 in (None, 3.0):
            with self.subTest(log_d25=log_d25):
                streams = [StringIO(), StringIO(), StringIO(
                    self._rc3_record(9.0, None, log_d25, "V", log_ae=1.0)), StringIO()]
                with patch("builtins.open", side_effect=streams):
                    reader = SphereReader(False, 1.5, 10.0, False, rc3_enabled=True)
                try:
                    self.assertIsNotNone(reader.read_star())
                    h = math.degrees(reader._galaxy_vectors[0].get_mag()) * 3600.0
                    # Ae=60秒角（直径）。円形銀河の半光量半径は30秒角。
                    self.assertAlmostEqual(30.0 / 1.6783469900166605, h)
                finally:
                    reader.close()

    def test_rc3_b_band_controls_shape_and_v_band_controls_star_count(self):
        results = []
        # 同じB=9で色だけ変更した場合と、同じ天体のV値を格納した場合。
        for magnitude, color, flag in ((9.0, 1.0, " "), (9.0, 0.0, " "),
                                       (8.0, 1.0, "V"), (8.0, None, "V")):
            streams = [StringIO(), StringIO(), StringIO(
                self._rc3_record(magnitude, color, 1.0, flag)), StringIO()]
            with patch("builtins.open", side_effect=streams):
                reader = SphereReader(False, 1.5, 10.0, False, rc3_enabled=True)
            try:
                stars = list(iter(reader.read_star, None))
                results.append((len(stars), reader._galaxy_vectors[0].get_mag()))
            finally:
                reader.close()
        self.assertEqual([6, 2, 6, 6], [count for count, _ in results])
        self.assertAlmostEqual(results[0][1], results[1][1])
        self.assertAlmostEqual(results[0][1], results[2][1])
        self.assertAlmostEqual(math.radians(10.0 / 3600.0), results[3][1])

    def test_generated_positions_follow_exponential_surface_brightness(self):
        record = list(self._rc3_record(0.0, 0.0, 1.0, log_ae=0.0))
        record[161:165] = "0.60"
        streams = [StringIO(), StringIO(), StringIO("".join(record)), StringIO()]
        with patch("builtins.open", side_effect=streams):
            reader = SphereReader(False, 1.5, 10.0, False, rc3_enabled=True)
        try:
            stars = list(iter(reader.read_star, None))
            major, minor, center = reader._galaxy_vectors
            major_unit, minor_unit = major.unit_vector(), minor.unit_vector()
            radii = []
            inside_circle = 0
            for star in stars:
                direction = MathVector.from_mag_lng_lat(
                    1.0, math.radians(star.p.radeg), math.radians(star.p.dedeg))
                # 接平面へ戻して、天球への変換を含めた実際の点の分布を調べる。
                cosine = direction.x * center.x + direction.y * center.y + direction.z * center.z
                offset = direction.mult_scalar(1.0 / cosine).minus(center)
                x = offset.x * major_unit.x + offset.y * major_unit.y + offset.z * major_unit.z
                y = offset.x * minor_unit.x + offset.y * minor_unit.y + offset.z * minor_unit.z
                radii.append(math.hypot(x / major.get_mag(), y / minor.get_mag()))
                inside_circle += offset.get_mag() <= math.radians(3.0 / 3600.0)
            # Gamma(2, 1): 平均半径2h、2h内の光量1-3exp(-2)。
            self.assertAlmostEqual(2.0, sum(radii) / len(radii), delta=0.05)
            self.assertAlmostEqual(1.0 - 3.0 * math.exp(-2.0),
                                   sum(r < 2.0 for r in radii) / len(radii), delta=0.02)
            self.assertAlmostEqual(0.5, inside_circle / len(radii), delta=0.02)
        finally:
            reader.close()

    def test_rc3_magnitude_flag_controls_color_conversion(self):
        # 最微等級10に対し、V=9なら2点、B=9・B-V=1ならV=8で6点。
        cases = [
            ("V", 1.0, 2),
            ("V", None, 2),
            (" ", 1.0, 6),
            ("M", 1.0, 6),
            ("S", 1.0, 6),
            ("v", 1.0, 6),
            (" ", None, 0),
        ]
        for bt_code, bv, expected_count in cases:
            with self.subTest(bt_code=bt_code, bv=bv):
                streams = [
                    StringIO(), StringIO(),
                    StringIO(self._rc3_record(9.0, bv, 1.0, bt_code)),
                    StringIO(),
                ]
                with patch("builtins.open", side_effect=streams):
                    reader = SphereReader(False, 1.5, 10.0, False, rc3_enabled=True)
                try:
                    stars = list(iter(reader.read_star, None))
                    self.assertEqual(expected_count, len(stars))
                    self.assertTrue(all(star.vmag == 10.0 for star in stars))
                finally:
                    reader.close()

    def test_faint_stars_and_rc3_can_be_enabled_independently(self):
        tycho_record = [" "] * 216
        tycho_record[41:46] = f"{8.5:5.2f}"
        tycho_record[51:63] = f"{12.0:12.6f}"
        tycho_record[64:76] = f"{34.0:12.6f}"
        hip_record = tycho_record.copy()
        hip_record[8:14] = "     1"
        for under_minimum in (False, True):
            for rc3_enabled in (False, True):
                with self.subTest(under_minimum=under_minimum, rc3_enabled=rc3_enabled):
                    streams = [
                        StringIO("".join(hip_record) + "\n"),
                        StringIO("".join(tycho_record) + "\n"),
                        StringIO(self._rc3_record(7.5, 0.5, 1.0)),
                        StringIO(),
                    ]
                    with patch("builtins.open", side_effect=streams):
                        reader = SphereReader(
                            False, 1.5, 7.5, under_minimum, rc3_enabled=rc3_enabled
                        )
                    try:
                        # 8.5等の恒星が確実に採用される乱数で分岐を確認する。
                        with patch.object(reader._random, "random", return_value=0.1):
                            stars = list(iter(reader.read_star, None))
                        self.assertEqual(2 * int(under_minimum) + int(rc3_enabled), len(stars))
                        self.assertTrue(all(star.vmag == 7.5 for star in stars))
                        self.assertEqual(
                            2 * int(under_minimum),
                            sum(star.p.radeg == 12.0 for star in stars),
                        )
                    finally:
                        reader.close()


class SphereReaderSettingsTests(unittest.TestCase):
    def test_config_switches_and_logs_are_independent(self):
        for under_minimum in (False, True):
            for rc3_enabled in (False, True):
                with self.subTest(under_minimum=under_minimum, rc3_enabled=rc3_enabled):
                    props = {
                        "star.above-maximum": "no",
                        "star.under-minimum": "yes" if under_minimum else "no",
                        "rc3.enabled": "yes" if rc3_enabled else "no",
                    }
                    output = StringIO()
                    with patch("transformer.SphereReader") as reader, patch("sys.stdout", output):
                        _init_sphere_reader(props)
                    reader.assert_called_once_with(
                        False, 1.5, 7.5, under_minimum, rc3_enabled=rc3_enabled
                    )
                    self.assertIn("最輝星より明るい星を、最輝星で疑似的に表現しません。", output.getvalue())
                    self.assertIn(
                        "最微星より暗い星を、最微星で疑似的に表現"
                        + ("します。" if under_minimum else "しません。"), output.getvalue()
                    )
                    self.assertIn(
                        "RC3カタログの系外銀河を、最微等級の点の集合で表現"
                        + ("します。" if rc3_enabled else "しません。"), output.getvalue()
                    )

    def test_omitted_rc3_setting_is_disabled_even_with_faint_stars_enabled(self):
        with patch("transformer.SphereReader") as reader, patch("sys.stdout", StringIO()):
            _init_sphere_reader({"star.under-minimum": "yes"})
        reader.assert_called_once_with(True, 1.5, 7.5, True, rc3_enabled=False)

    def test_interactive_switches_are_independent(self):
        for under_minimum in (False, True):
            for rc3_enabled in (False, True):
                with self.subTest(under_minimum=under_minimum, rc3_enabled=rc3_enabled):
                    answers = ["n", "", "", "y" if under_minimum else "n", "y" if rc3_enabled else ""]
                    with patch("builtins.input", side_effect=answers), \
                            patch("transformer.SphereReader") as reader, \
                            patch("sys.stdout", StringIO()):
                        _init_sphere_reader(None)
                    reader.assert_called_once_with(
                        False, 1.5, 7.5, under_minimum, rc3_enabled=rc3_enabled
                    )


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
            self.assertEqual("210.0mm", root.attrib["width"])
            self.assertEqual("297.0mm", root.attrib["height"])
            first_x, first_y = map(float, polygons[0].attrib["points"].split()[0].split(","))
            self.assertAlmostEqual(105 * 96 / 25.4, first_x)
            self.assertAlmostEqual(74.25 * 96 / 25.4, first_y)

    def test_polygon_generation_uses_svg_output_category(self):
        with tempfile.TemporaryDirectory() as directory, patch(
            "sys.stdout", new_callable=StringIO
        ):
            star_writer = PlateWriterSVG(
                1, 1, 37.5, True, "star-", False, Path(directory, "star_SVG")
            )
            try:
                _write_assignment_polygons(
                    self.transformer,
                    star_writer,
                    {
                        "output.directory": directory,
                        "polygon.enabled": "yes",
                        "polygon.samples": "12",
                    },
                )
            finally:
                star_writer.close()

            self.assertTrue(Path(directory, "polygon_SVG", "polygon-0.svg").is_file())
            self.assertFalse(Path(directory, "star_SVG", "polygon-0.svg").exists())

    @unittest.skipUnless(_HAS_REPORTLAB, "ReportLab is required for PDF writer tests")
    def test_polygon_generation_uses_pdf_output_category_with_pdf_writer(self):
        with tempfile.TemporaryDirectory() as directory, patch(
            "sys.stdout", new_callable=StringIO
        ):
            star_writer = PlateWriterPDF(
                1, 1, 37.5, True, "star-", False, Path(directory, "star_pdf")
            )
            try:
                _write_assignment_polygons(
                    self.transformer,
                    star_writer,
                    {
                        "output.directory": directory,
                        "polygon.enabled": "yes",
                        "polygon.samples": "12",
                    },
                )
            finally:
                star_writer.close()

            output = Path(directory, "polygon_pdf", "polygon-0.pdf")
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))
            self.assertFalse(Path(directory, "polygon_SVG", "polygon-0.svg").exists())


def _plate_star(dir_=0, index=0, xmm=0.0, ymm=0.0, radius=0.5):
    star = PlateStar()
    star.p = PlatePosition()
    star.p.dir = dir_
    star.p.index = index
    star.p.xmm = xmm
    star.p.ymm = ymm
    star.rmm = radius
    return star


class PlateFrameBoundaryTests(unittest.TestCase):
    def test_circular_writer_rejects_square_corner_outside_radius(self):
        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterSVG(1, 1, 40, True, "circle-", False, directory)
            writer.write_star(_plate_star(xmm=35, ymm=35))
            writer.write_star(_plate_star(xmm=24, ymm=32))
            writer.close()

            root = ET.parse(Path(directory, "circle-0.svg")).getroot()
            star_holes = [
                element
                for element in root.findall("{http://www.w3.org/2000/svg}circle")
                if "stroke" not in element.attrib
            ]
            self.assertEqual(1, len(star_holes))
            self.assertEqual("129.0mm", star_holes[0].attrib["cx"])
            self.assertEqual("106.25mm", star_holes[0].attrib["cy"])

    def test_rectangular_writer_keeps_square_corner(self):
        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterSVG(1, 1, 40, False, "rectangle-", False, directory)
            writer.write_star(_plate_star(xmm=35, ymm=35))
            writer.close()

            root = ET.parse(Path(directory, "rectangle-0.svg")).getroot()
            star_holes = root.findall("{http://www.w3.org/2000/svg}circle")
            self.assertEqual(1, len(star_holes))


@unittest.skipUnless(_HAS_REPORTLAB, "ReportLab is required for PDF writer tests")
class PDFWriterTests(unittest.TestCase):
    def test_pdf_replaces_existing_file_without_opening_it_for_writing(self):
        import builtins
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, "print-0.pdf")
            output.write_bytes(b"existing PDF")
            original_open = builtins.open

            def guarded_open(file, mode="r", *args, **kwargs):
                if str(file) == str(output) and "w" in mode:
                    raise TimeoutError(60, "Operation timed out")
                return original_open(file, mode, *args, **kwargs)

            writer = PlateWriterPDF(1, 1, 37.5, True, "print-", False, directory)
            writer.write_star(_plate_star())
            with patch("builtins.open", side_effect=guarded_open):
                writer.close()
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))
            self.assertEqual([output], list(Path(directory).iterdir()))

    def test_failed_pdf_replace_preserves_existing_file_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, "print-0.pdf")
            output.write_bytes(b"existing PDF")
            writer = PlateWriterPDF(1, 1, 37.5, True, "print-", False, directory)
            writer.write_star(_plate_star())
            with patch("plate_writer.os.replace", side_effect=TimeoutError(60, "Operation timed out")):
                with self.assertRaises(OSError) as error:
                    writer.close()
            self.assertEqual(str(output), error.exception.filename)
            self.assertEqual(b"existing PDF", output.read_bytes())
            self.assertEqual([output], list(Path(directory).iterdir()))
            self.assertIsNone(writer.outs)

    def test_pdf_writer_replaces_postscript_writer(self):
        self.assertEqual({"SVG", "PDF"}, {writer_type.name for writer_type in PlateWriterType})

        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterPDF(1, 1, 37.5, True, "print-", False, directory)
            writer.write_star(_plate_star())
            writer.close()

            output = Path(directory, "print-0.pdf")
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))
            self.assertFalse(Path(directory, "print-0.ps").exists())

    def test_pdf_factory_is_used_by_star_transformer(self):
        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            props = {"output.directory": directory}
            star_writer = _init_plate_writer(props, PlateWriterType.PDF)
            try:
                self.assertIsInstance(star_writer, PlateWriterPDF)
                self.assertEqual(Path(directory, "star_pdf"), Path(star_writer.output_dir))
            finally:
                star_writer.close()

    @unittest.skipUnless(_HAS_PYPDF, "pypdf is required for PDF structure tests")
    def test_pdf_writer_draws_assignment_polygon(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterPDF(1, 1, 37.5, True, "polygon-", True, directory)
            writer.write_assignment_polygon(0, 0, [(0, 0), (10, 0), (0, 10)])
            writer.close()

            output = Path(directory, "polygon-0.pdf")
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))
            page = PdfReader(output).pages[0]
            operators = [operator for _, operator in page.get_contents().operations]
            self.assertIn(b"m", operators)
            self.assertIn(b"l", operators)
            self.assertIn(b"h", operators)
            self.assertIn(b"B*", operators)
            self.assertEqual([], list(page.images))

    @unittest.skipUnless(_HAS_PYPDF, "pypdf is required for PDF structure tests")
    def test_pdf_circular_writer_uses_radial_boundary(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterPDF(1, 1, 40, True, "print-", False, directory)
            writer.write_star(_plate_star(xmm=35, ymm=35))
            writer.write_star(_plate_star(xmm=24, ymm=32))
            writer.close()

            operations = PdfReader(Path(directory, "print-0.pdf")).pages[0].get_contents().operations
            # 北天・南天の枠2円と、半径40 mm上の星穴1円だけが残る。
            self.assertEqual(12, sum(operator == b"c" for _, operator in operations))

    @unittest.skipUnless(_HAS_PYPDF, "pypdf is required for PDF structure tests")
    def test_pdf_is_print_sized_vector_and_supports_constellations(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterPDF(1, 1, 37.5, True, "print-", True, directory)
            writer.write_star(_plate_star(dir_=0, xmm=-5.0, ymm=3.0))
            writer.write_star(_plate_star(dir_=1, xmm=4.0, ymm=-2.0))

            constellation = PlateConstellation()
            constellation.name = "Orion"
            constellation.p = PlatePosition()
            constellation.p.xmm = 0.0
            constellation.p.ymm = 10.0
            line_start = PlatePosition()
            line_start.xmm = -12.0
            line_end = PlatePosition()
            line_end.xmm = 12.0
            constellation.ll = [[line_start, line_end]]
            writer.write_constellation(constellation)
            writer.close()

            reader = PdfReader(Path(directory, "print-0.pdf"))
            self.assertEqual(1, len(reader.pages))
            page = reader.pages[0]
            point_to_mm = 25.4 / 72
            self.assertAlmostEqual(210, float(page.mediabox.width) * point_to_mm, places=3)
            self.assertAlmostEqual(297, float(page.mediabox.height) * point_to_mm, places=3)
            self.assertEqual(["N0", "S0", "Orion"], page.extract_text().split())
            self.assertEqual([], list(page.images))
            self.assertIn(b"l", [operator for _, operator in page.get_contents().operations])


class OutputFormatCLITests(unittest.TestCase):
    def test_interactive_output_types_can_select_all_four_outputs(self):
        with (
            patch("builtins.input", side_effect=["", "yes", "true", "y"]),
            patch("sys.stdout", new_callable=StringIO),
        ):
            star_types, polygon_types = _select_interactive_output_types()

        self.assertEqual(
            [PlateWriterType.SVG, PlateWriterType.PDF], star_types
        )
        self.assertEqual(
            [PlateWriterType.SVG, PlateWriterType.PDF], polygon_types
        )

    def test_interactive_output_types_default_to_star_svg_only(self):
        with (
            patch("builtins.input", side_effect=["", "", "", ""]),
            patch("sys.stdout", new_callable=StringIO),
        ):
            star_types, polygon_types = _select_interactive_output_types()

        self.assertEqual([PlateWriterType.SVG], star_types)
        self.assertEqual([], polygon_types)

    def test_interactive_main_routes_stars_and_polygons_to_both_formats(self):
        from transformer import main as star_main

        transformer = MagicMock()
        reader = MagicMock()
        svg_writer = MagicMock()
        pdf_writer = MagicMock()
        writers = {
            PlateWriterType.SVG: svg_writer,
            PlateWriterType.PDF: pdf_writer,
        }
        both_types = [PlateWriterType.SVG, PlateWriterType.PDF]
        with (
            patch("transformer._init_transformer", return_value=transformer),
            patch("transformer._init_sphere_reader", return_value=reader),
            patch(
                "transformer._select_interactive_output_types",
                return_value=(both_types, both_types),
            ),
            patch("transformer._init_plate_writers", return_value=writers),
            patch("transformer._init_enlarge_rate", return_value=0.0),
            patch("transformer._write_assignment_polygons") as write_polygons,
            patch("builtins.input", return_value="n"),
            patch("sys.stdout", new_callable=StringIO),
        ):
            star_main([])

        output_writer = transformer.process_stars.call_args.args[1]
        self.assertIsInstance(output_writer, _PlateWriterGroup)
        self.assertEqual([svg_writer, pdf_writer], output_writer.output_writers)
        self.assertEqual(
            [
                ((transformer, svg_writer, None), {"force": True}),
                ((transformer, pdf_writer, None), {"force": True}),
            ],
            [(call.args, call.kwargs) for call in write_polygons.call_args_list],
        )
        svg_writer.close.assert_called_once_with()
        pdf_writer.close.assert_called_once_with()

    def test_interactive_svg_writer_can_invert_plate_and_star_colors(self):
        with (
            tempfile.TemporaryDirectory() as directory,
            patch("transformer.DEFAULT_OUTPUT_DIR", directory),
            patch("builtins.input", side_effect=["1", "1", "0", "yes"]),
            patch("sys.stdout", new_callable=StringIO) as stdout,
        ):
            writer = _init_plate_writer(None, PlateWriterType.SVG)
            try:
                self.assertTrue(writer.invert_color)
                writer.write_star(_plate_star())
            finally:
                writer.close()

            root = ET.parse(Path(directory, "star_SVG", "star-0.svg")).getroot()
            circles = root.findall("{http://www.w3.org/2000/svg}circle")
            self.assertEqual("black", circles[0].attrib["fill"])
            self.assertEqual("white", circles[-1].attrib["fill"])
            self.assertIn("原板を黒色、星を白色", stdout.getvalue())

    def test_interactive_svg_writer_keeps_non_inverted_default(self):
        with (
            tempfile.TemporaryDirectory() as directory,
            patch("transformer.DEFAULT_OUTPUT_DIR", directory),
            patch("builtins.input", side_effect=["", "", "", ""]),
            patch("sys.stdout", new_callable=StringIO),
        ):
            writer = _init_plate_writer(None, PlateWriterType.SVG)
            try:
                self.assertFalse(writer.invert_color)
            finally:
                writer.close()

    def test_output_categories_are_created_below_the_configured_root(self):
        self.assertEqual(
            Path("custom-output", "star_SVG"),
            Path(categorized_output_dir("custom-output", "star_SVG")),
        )
        with tempfile.TemporaryDirectory() as directory, patch(
            "sys.stdout", new_callable=StringIO
        ):
            writer = _init_plate_writer(
                {"output.directory": directory}, PlateWriterType.SVG
            )
            try:
                self.assertEqual(Path(directory, "star_SVG"), Path(writer.output_dir))
                writer.write_star(_plate_star())
            finally:
                writer.close()
            self.assertTrue(Path(directory, "star_SVG", "star-0.svg").is_file())

    def test_deprecated_postscript_option_stops_star_command(self):
        from transformer import main as star_main

        with (patch("sys.stdout", new_callable=StringIO),
              patch("sys.stderr", new_callable=StringIO) as stderr):
            self.assertEqual(2, star_main(["-PS"]))
            self.assertIn("-PDF", stderr.getvalue())

    def test_pdf_mode_processes_constellations_in_star_command(self):
        from transformer import main as star_main

        transformer = MagicMock()
        reader = MagicMock()
        writer = MagicMock()
        with (patch("transformer._init_transformer", return_value=transformer),
              patch("transformer._init_sphere_reader", return_value=reader),
              patch("transformer._init_plate_writer", return_value=writer),
              patch("transformer._init_enlarge_rate", return_value=0.0),
              patch("transformer._write_assignment_polygons"),
              patch("builtins.input", return_value="y"),
              patch("sys.stdout", new_callable=StringIO) as stdout):
            star_main(["-PDF"])

        self.assertIn("星座を処理しますか。(y/N)", stdout.getvalue())
        transformer.process_constellations.assert_called_once_with(reader, writer)
if __name__ == "__main__":
    unittest.main()

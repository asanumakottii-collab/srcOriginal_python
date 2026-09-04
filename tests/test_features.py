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
from galaxy_transformer import _init_plate_writer as _init_galaxy_plate_writer
from models import PlateConstellation, PlatePosition, PlateStar, SpherePosition, SphereStar
from plate_polygon import PlateAssignmentPolygonGenerator
from plate_writer import PlateWriterPDF, PlateWriterSVG, PlateWriterType, categorized_output_dir
from sphere_reader import SphereReader
from transformer import (
    _PlateWriterGroup,
    Transformer,
    _init_enlarge_rate,
    _init_plate_writer,
    _select_interactive_output_types,
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


class SphereReaderTests(unittest.TestCase):
    @staticmethod
    def _rc3_record(btmag, bv, log_d25):
        record = [" "] * 256
        record[0:2] = "12"
        record[2:4] = "00"
        record[4:8] = "00.0"
        record[9] = "+"
        record[10:12] = "20"
        record[12:14] = "00"
        record[14:16] = "00"
        record[151:155] = f"{log_d25:4.2f}"
        record[157:160] = ".00"
        record[185:188] = "  0"
        record[189:194] = f"{btmag:5.2f}"
        record[252:256] = f"{bv:4.2f}"
        return "".join(record) + "\n"

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

    def test_invalid_galaxy_scale_is_skipped(self):
        invalid_record = self._rc3_record(10.0, 0.5, 3.0)
        valid_record = self._rc3_record(10.0, 0.5, 1.0)
        streams = [
            StringIO(),
            StringIO(),
            StringIO(invalid_record + valid_record),
            StringIO(),
        ]

        with patch("builtins.open", side_effect=streams):
            reader = SphereReader(False, 7.5, 10.0, True)
            star = reader.read_star()

        self.assertIsNotNone(star)
        self.assertTrue(math.isfinite(star.p.radeg))
        self.assertTrue(math.isfinite(star.p.dedeg))


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
    def test_pdf_writer_replaces_postscript_writer(self):
        self.assertEqual({"SVG", "PDF"}, {writer_type.name for writer_type in PlateWriterType})

        with tempfile.TemporaryDirectory() as directory:
            writer = PlateWriterPDF(1, 1, 37.5, True, "print-", False, directory)
            writer.write_star(_plate_star())
            writer.close()

            output = Path(directory, "print-0.pdf")
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))
            self.assertFalse(Path(directory, "print-0.ps").exists())

    def test_pdf_factories_are_used_by_both_transformers(self):
        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            props = {"output.directory": directory}
            star_writer = _init_plate_writer(props, PlateWriterType.PDF)
            galaxy_writer = _init_galaxy_plate_writer(props, PlateWriterType.PDF)
            try:
                self.assertIsInstance(star_writer, PlateWriterPDF)
                self.assertIsInstance(galaxy_writer, PlateWriterPDF)
                self.assertEqual(Path(directory, "star_pdf"), Path(star_writer.output_dir))
                self.assertEqual(Path(directory, "galaxy"), Path(galaxy_writer.output_dir))
            finally:
                star_writer.close()
                galaxy_writer.close()

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

    def test_deprecated_postscript_option_stops_both_commands(self):
        from galaxy_transformer import main as galaxy_main
        from transformer import main as star_main

        for main in (star_main, galaxy_main):
            with (patch("sys.stdout", new_callable=StringIO),
                  patch("sys.stderr", new_callable=StringIO) as stderr):
                self.assertEqual(2, main(["-PS"]))
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

    def test_pdf_mode_processes_constellations_in_galaxy_command(self):
        from galaxy_transformer import main as galaxy_main

        transformer = MagicMock()
        reader = MagicMock()
        writer = MagicMock()
        with (patch("galaxy_transformer._init_galaxy_transformer", return_value=transformer),
              patch("galaxy_transformer._init_sphere_reader", return_value=reader),
              patch("galaxy_transformer._init_plate_writer", return_value=writer),
              patch("builtins.input", return_value="y"),
              patch("sys.stdout", new_callable=StringIO) as stdout):
            galaxy_main(["-PDF"])

        self.assertIn("星座を処理しますか。(y/N)", stdout.getvalue())
        transformer.process_constellations.assert_called_once_with(reader, writer)


if __name__ == "__main__":
    unittest.main()

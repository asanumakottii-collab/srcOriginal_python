# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012 東京大学地文研究会天文部
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
"""天球の星から、複数の投影機ユニットに配置された原板の星への変換。"""

import math
import sys

from basic_transformer import BasicTransformer
from config import load_properties
from mathvector import MathVector
from models import PlatePosition, SpherePosition
from plate_polygon import PlateAssignmentPolygonGenerator
from plate_writer import DEFAULT_OUTPUT_DIR, PlateWriterPS, PlateWriterSVG, PlateWriterType, resolve_output_path
from sphere_reader import SphereReader
from unit_arrangement import CustomUnitArrangement, InvalidUnitArrangementException, StandardUnitArrangement

# 浮動小数点の計算で生じる誤差よりも十分大きな数 ここでは 2^-16
_EPS = 0.0000152587890625


def _get_output_dir(props):
    if props is None:
        return DEFAULT_OUTPUT_DIR
    return props.get("output.directory", DEFAULT_OUTPUT_DIR)


class _Unit:
    """ユニットの座標系の単位ベクトルとレンズ位置。

    vectors[0]:原板のX軸方向, vectors[1]:原板のY軸方向, vectors[2]:原板の法線ベクトル(Z軸方向)
    """

    __slots__ = ("vectors", "lens")

    def __init__(self):
        self.vectors = [None, None, None]
        self.lens = 0.0


class Transformer(BasicTransformer):
    def __init__(self, radius, sphere, projector, plate, focus, ar):
        """
        座標変換のパラメータを設定します。
        :param radius: 0等星の半径
        :param sphere: ドームの半径
        :param projector: ドームの中心と投影機中心との間の距離
        :param plate: 投影機中心と原版中心の間の距離
        :param focus: 焦点の距離
        :param ar: ユニットの配置
        """
        super().__init__()
        self.radius = radius
        self.sphere = sphere
        self.projector = projector
        self.plate = plate
        self.base_magnitude = 0.

        y_axis = MathVector(0, 1, 0)
        z_axis = MathVector(0, 0, 1)
        self._units = [None, None]
        for i in range(2):
            number_of_units = ar.number_of_units(i)
            self._units[i] = [None] * number_of_units
            for j in range(number_of_units):
                unit = _Unit()
                self._units[i][j] = unit
                unit.vectors[2] = ar.get_unit_vector(i, j)
                if unit.vectors[2].z == (1. if i == 0 else -1):  # 北極と南極の場合
                    unit.vectors[1] = y_axis  # 原板のY軸方向: (0,1,0)
                    unit.vectors[0] = unit.vectors[1].cross(unit.vectors[2])  # 原板のX軸方向
                else:  # 北極と南極以外の場合
                    unit.vectors[0] = z_axis.cross(unit.vectors[2]).unit_vector()  # 原板のX軸方向
                    unit.vectors[1] = unit.vectors[2].cross(unit.vectors[0])  # 原板のY軸方向

                # ドームの中心をO、投影機中心をP、像の位置をIとしたとき、三角形OPIの∠OPIについて
                # 余弦定理を適用すると得られる方程式 OI^2=OP^2+PI^2-2OP*PD*cos∠OPI を解く。
                a = 1.
                b = -projector * math.cos((math.pi / 2 if i == 0 else -math.pi / 2) + unit.vectors[2].get_lat())
                c = projector * projector - sphere * sphere
                d = b * b - a * c
                image = (-b + math.sqrt(d)) / a - plate  # 原板の中心と像の中心との間の距離

                # unit.lens(未知数)をxとしたとき、レンズの公式より得られる方程式
                # 1/focus = 1/x + 1/(image-x) を解く。
                a = 1.
                b = -image / 2
                c = image * focus
                d = b * b - a * c
                unit.lens = (-b - math.sqrt(d)) / a
        self._sort_units()

    def _sort_units(self):
        """ユニットを、緯度、経度の偶奇、経度の順でソートします。"""
        for i in range(len(self._units)):
            units_i = self._units[i]
            for j in range(len(units_i)):
                for k in range(j + 1, len(units_i)):
                    unit_j, unit_k = units_i[j], units_i[k]
                    sign = 1. if i == 0 else -1.
                    if (sign * unit_j.vectors[2].get_lat()) - _EPS > (sign * unit_k.vectors[2].get_lat()):
                        continue
                    if (sign * unit_j.vectors[2].get_lat()) + _EPS > (sign * unit_k.vectors[2].get_lat()):
                        odd_j = int(math.floor((unit_j.vectors[2].get_lng() + _EPS + math.pi * 2) * 5 / math.pi)) % 2
                        odd_k = int(math.floor((unit_k.vectors[2].get_lng() + _EPS + math.pi * 2) * 5 / math.pi)) % 2
                        if odd_j < odd_k:
                            continue
                        if odd_j == odd_k:
                            if ((unit_j.vectors[2].get_lng() + _EPS + math.pi * 2) % (math.pi * 2)
                                    < (unit_k.vectors[2].get_lng() + _EPS + math.pi * 2) % (math.pi * 2)):
                                continue
                    units_i[j], units_i[k] = unit_k, unit_j

    def transform_unit(self, sp, dir_, index):
        """レンズの中心を通る光は直進する性質を利用して、天球上の位置を指定されたユニットの原板上の位置に変換します。"""
        unit = self._units[dir_][index]
        pp = PlatePosition()
        vector = sp.to_vector(self.sphere)
        vector.z += -self.projector if dir_ == 0 else self.projector
        vector = vector.minus(unit.vectors[2].mult_scalar(self.plate + unit.lens))
        vector = vector.transform(unit.vectors)
        vector = vector.mult_scalar(-unit.lens / abs(vector.z))
        pp.xmm = vector.x
        pp.ymm = vector.y
        pp.dir = dir_
        pp.index = index
        # レンズによる反転と座標系の変更が相殺されるため、そのままでよい。
        if vector.z < 0:
            return pp
        return None

    def transform(self, sp):
        """レンズの中心を通る光は直進する性質を利用して、天球上の位置をもっとも適切と判定されたユニットの原板上の位置に変換します。"""
        pp = None
        e = float("-inf")
        for i in range(len(self._units)):
            for j in range(len(self._units[i])):
                pp1 = self.transform_unit(sp, i, j)
                e1 = self._units[pp1.dir][pp1.index].lens / pp1.distance_from_center() if pp1 is not None else float("-inf")
                if pp1 is not None and e1 >= e:
                    pp = pp1
                    e = e1
        return pp

    def inverse_transform_unit(self, dir_, index, xmm, ymm):
        """原板上の座標から、その光線が到達する天球上の位置を逆算します。

        担当星域ポリゴンの計算に使うほか、transform_unit() の逆変換として
        利用できます。指定点からレンズ中心を通る光線とドーム球面の交点を
        返します。
        """
        unit = self._units[dir_][index]
        projector_center = MathVector(0, 0, self.projector if dir_ == 0 else -self.projector)
        lens_center = projector_center.plus(
            unit.vectors[2].mult_scalar(self.plate + unit.lens)
        )
        ray = (
            unit.vectors[0].mult_scalar(-xmm)
            .plus(unit.vectors[1].mult_scalar(-ymm))
            .plus(unit.vectors[2].mult_scalar(unit.lens))
            .unit_vector()
        )
        lens_dot_ray = lens_center.dot(ray)
        discriminant = (
            lens_dot_ray * lens_dot_ray
            + self.sphere * self.sphere
            - lens_center.dot(lens_center)
        )
        if discriminant < 0:
            return None
        distance = -lens_dot_ray + math.sqrt(discriminant)
        return SpherePosition.from_vector(lens_center.plus(ray.mult_scalar(distance)))

    def number_of_units(self, dir_):
        """指定半球に配置されたユニット数を返します。"""
        return len(self._units[dir_])

    def assigned_unit(self, sp):
        """天球上の位置を担当する ``(半球, ユニット番号)`` を返します。

        transform() と同じ選択基準を、原板座標オブジェクトを生成せずに評価
        します。大量の境界点を調べる担当星域ポリゴン生成向けの高速版です。
        """
        sphere_vector = sp.to_vector(self.sphere)
        best = None
        best_score = float("-inf")
        for dir_, units in enumerate(self._units):
            projector_center = MathVector(
                0, 0, self.projector if dir_ == 0 else -self.projector
            )
            for index, unit in enumerate(units):
                relative = sphere_vector.minus(projector_center).minus(
                    unit.vectors[2].mult_scalar(self.plate + unit.lens)
                )
                x = relative.dot(unit.vectors[0])
                y = relative.dot(unit.vectors[1])
                z = relative.dot(unit.vectors[2])
                if z <= 0:
                    continue
                distance = math.hypot(x, y)
                score = float("inf") if distance == 0 else z / distance
                if score >= best_score:
                    best = (dir_, index)
                    best_score = score
        return best

    def write_unit_position(self, size):
        """ユニットの配置を画像化します。"""
        from PIL import Image, ImageDraw

        radius = size * 0.8 / 2
        image = Image.new("RGB", (size, 2 * size), "white")
        draw = ImageDraw.Draw(image)
        for i in range(2):
            cx, cy = size / 2, size / 2 + i * size
            draw.ellipse([cx - radius, cy - radius, cx + radius, cy + radius], outline="black")
            prefix = "N" if i == 0 else "S"
            s = 1 if i == 0 else -1
            draw.line([cx - radius, cy, cx - radius - 10, cy], fill="black")
            draw.line([cx, cy - s * radius, cx, cy - s * (radius + 10)], fill="black")
            draw.line([cx + radius, cy, cx + radius + 10, cy], fill="black")
            draw.line([cx, cy + s * radius, cx, cy + s * (radius + 10)], fill="black")
            draw.text((cx - radius - 20, cy), "0h", fill="black", anchor="rm")
            draw.text((cx, cy + s * (radius + 20)), "6h", fill="black", anchor="ma" if s > 0 else "mb")
            draw.text((cx + radius + 20, cy), "12h", fill="black", anchor="lm")
            draw.text((cx, cy - s * (radius + 20)), "18h", fill="black", anchor="mb" if s > 0 else "ma")
            for j, unit in enumerate(self._units[i]):
                x = cx - radius * unit.vectors[2].x
                y = cy + s * radius * unit.vectors[2].y
                draw.line([x - 10, y - 10, x + 10, y + 10], fill="black")
                draw.line([x - 10, y + 10, x + 10, y - 10], fill="black")
                draw.text((x, y), f"{prefix}{j}", fill="black")
        return image


def _usage():
    """このプログラムの使い方を表示します。"""
    print("Usage:")
    print("\tpython transformer.py [options]")
    print("Options:")
    print("\t-PS\t原板データをPostScript形式で出力します。")
    print("\t-f [configFile]\tプラネタリウムの設定を configFile から読み込みます。")
    print("\t--polygons\t原盤ごとの担当星域ポリゴンをSVG出力します。")
    print("\t-h,-help\tこのメッセージを表示します。")


def _init_transformer(props):
    """設定を読み込み、Transformer のインスタンスを作成します。"""
    ar = None
    if props is None:  # interactive mode
        print("後で何分の一の倍率にしますか。(default=1) ")
        scale = BasicTransformer.parse_double_with_default(input(), 1.)
        print(f"0等星の半径は何mmですか。{scale}分の一の倍率をかける前の値を指定してください。(default=0.25) ")
        radius = BasicTransformer.parse_double_with_default(input(), 0.25)
        print("ドームの半径は何mmですか。(default=5000)")
        sphere = BasicTransformer.parse_double_with_default(input(), 5000.)
        print("ドームの中心と投影機の中心との間の距離は何mmですか。(default=300)")
        projector = BasicTransformer.parse_double_with_default(input(), 300.)
        print("投影機の中心と原板の中心の間の距離は何mmですか。(default=200)。")
        plate = BasicTransformer.parse_double_with_default(input(), 200.)
        print("焦点距離は何mmですか。(default=50)")
        focus = BasicTransformer.parse_double_with_default(input(), 50.)
        print("ユニットを正十二面体の面方向12箇所に配置しますか。(Y/n)")
        has_unit_12 = input().lower() != "n"
        print("ユニットを正十二面体の頂点方向20箇所に配置しますか。(y/N)")
        has_unit_20 = input().lower() == "y"
        print("ユニットを正十二面体の辺方向30箇所に配置しますか。(Y/n)")
        has_unit_30 = input().lower() != "n"
        print("ユニットを正十二面体の対称性をもつ最適な方向60箇所に配置しますか。(y/N)")
        has_unit_60 = input().lower() == "y"
        print("加えて、ユニットを赤道付近に配置しますか。(Y/n)")
        has_unit_extra = input().lower() != "n"
        ar = StandardUnitArrangement(has_unit_12, has_unit_20, has_unit_30, has_unit_60, has_unit_extra)
    else:  # non-interactive mode
        scale = BasicTransformer.parse_double_with_default(props.get("scale", ""), 1.)
        radius = BasicTransformer.parse_double_with_default(props.get("star-radius"), 0.25)
        sphere = BasicTransformer.parse_double_with_default(props.get("dome-radius"), 5000.)
        projector = BasicTransformer.parse_double_with_default(props.get("projector-distance"), 300.)
        plate = BasicTransformer.parse_double_with_default(props.get("plate-distance"), 200.)
        focus = BasicTransformer.parse_double_with_default(props.get("focal-length"), 50.)
        customunits = BasicTransformer.parse_boolean_with_default(props.get("unit.custom", ""), False)
        print("設定を読み込みました:")
        print(f"\t後で {scale} 分の一の倍率にします。")
        print(f"\t0等星の半径は {radius} mm です。")
        print(f"\tドームの半径は {sphere} mm です。")
        print(f"\tドームの中心と投影機の中心の間の距離は {projector} mm です。")
        print(f"\t投影機の中心と原板の中心の間の距離は {plate} mm です。")
        print(f"\t焦点距離は {focus} mm です。")
        if not customunits:
            has_unit_12 = BasicTransformer.parse_boolean_with_default(props.get("unit-12"), True)
            has_unit_20 = BasicTransformer.parse_boolean_with_default(props.get("unit-20"), False)
            has_unit_30 = BasicTransformer.parse_boolean_with_default(props.get("unit-30"), True)
            has_unit_60 = BasicTransformer.parse_boolean_with_default(props.get("unit-60"), True)
            has_unit_extra = BasicTransformer.parse_boolean_with_default(props.get("unit-extra"), True)
            ar = StandardUnitArrangement(has_unit_12, has_unit_20, has_unit_30, has_unit_60, has_unit_extra)
            print("\tユニットを正十二面体の面方向12箇所に配置" + ("します。" if has_unit_12 else "しません。"))
            print("\tユニットを正十二面体の頂点方向20箇所に配置" + ("します。" if has_unit_20 else "しません。"))
            print("\tユニットを正十二面体の辺方向30箇所に配置" + ("します。" if has_unit_30 else "しません。"))
            print("\tユニットを正十二面体の対称性をもつ最適な方向60箇所に配置" + ("します。" if has_unit_60 else "しません。"))
            print("\t加えて、ユニットを赤道付近に配置" + ("します。" if has_unit_extra else "しません。"))
        else:
            try:
                ar = CustomUnitArrangement.from_properties(props, "units")
                print("\t設定ファイルに記述されたユニットの位置を採用します。")
                print(f"\t北半球:{ar.number_of_units(0)}個")
                print(f"\t南半球:{ar.number_of_units(1)}個")
            except InvalidUnitArrangementException as e:
                print(f"設定ファイルに記述されたユニットの配置がおかしいよ（；；）\n{e}", file=sys.stderr)
                return None
    return Transformer(radius, sphere * scale, projector * scale, plate * scale, focus * scale, ar)


def _init_sphere_reader(props):
    """設定を読み込み、SphereReader のインスタンスを作成します。"""
    excluding_stars = None
    if props is None:  # interactive mode
        print("最輝星より明るい星を、最輝星で疑似的に表現しますか。(Y/n) ")
        above_maximum = input().lower() != "n"
        print("最輝星は何等星ですか。(default=1.5) ")
        maximum = BasicTransformer.parse_double_with_default(input(), 1.5)
        print("青セロハンを使わない最微星は何等星ですか。(default=3.5) ")
        extremum = BasicTransformer.parse_double_with_default(input(), 3.5)
        print("青セロハンを使った最微星は何等星ですか。(default=7.5) ")
        minimum = BasicTransformer.parse_double_with_default(input(), 7.5)
        print("最微星より暗い星を、最微星で疑似的に表現しますか。(y/N) ")
        under_minimum = input().lower() == "y"
    else:  # non-interactive mode
        above_maximum = BasicTransformer.parse_boolean_with_default(props.get("star.above-maximum", ""), True)
        maximum = BasicTransformer.parse_double_with_default(props.get("star.maximum"), 1.5)
        extremum = BasicTransformer.parse_double_with_default(props.get("star.extremum"), 3.5)
        minimum = BasicTransformer.parse_double_with_default(props.get("star.minimum"), 7.5)
        under_minimum = BasicTransformer.parse_boolean_with_default(props.get("star.under-minimum"), False)
        excluding_stars = props.get("star.excluding", "").strip()
        print("\t最輝星より明るい星を、最輝星で疑似的に表現" + ("します。" if under_minimum else "しません。"))
        print(f"\t最輝星は {maximum} 等星です。")
        print(f"\t青セロハンを使わない最微星は {extremum} 等星です。")
        print(f"\t青セロハンを使った最微星は {minimum} 等星です。")
        print("\t最微星より暗い星を、最微星で疑似的に表現" + ("します。" if under_minimum else "しません。"))
        if excluding_stars == "":
            print("\t除外する星はありません。")
        else:
            print(f"\t除外する星は {excluding_stars} です。")
    reader = SphereReader(above_maximum, maximum, extremum, minimum, under_minimum)
    if excluding_stars is not None:
        for hip_num in excluding_stars.split(" "):
            if hip_num != "":
                print(f"{hip_num} を除外")
                reader.exclude_star(hip_num)
    return reader


def _init_plate_writer(props, writer_type):
    """設定を読み込み、PlateWriter のインスタンスを作成します。"""
    filename_prefix = "star-"
    invert_color = False
    output_dir = _get_output_dir(props)
    if props is None:  # interactive mode
        print("横に各天のユニットをいくつ配置しますか。(default=1) ")
        column = BasicTransformer.parse_int_with_default(input(), 1)
        print("縦に各天のユニットをいくつ配置しますか。(default=1) ")
        row = BasicTransformer.parse_int_with_default(input(), 1)
        print("枠の半径は何mmですか。(default=0.) ")
        frame = BasicTransformer.parse_double_with_default(input(), 0.)
    else:  # non-interactive mode
        column = BasicTransformer.parse_int_with_default(props.get("plate.column"), 1)
        row = BasicTransformer.parse_int_with_default(props.get("plate.row"), 1)
        frame = BasicTransformer.parse_double_with_default(props.get("plate.frame-radius"), 0.)
        filename_prefix = props.get("file.prefix", filename_prefix)
        invert_color = BasicTransformer.parse_boolean_with_default(props.get("color.invert", ""), False)
        print(f"\t横に各天のユニットを {column} 個配置します。")
        print(f"\t縦に各天のユニットを {row} 個配置します。")
        print(f"\t枠の半径は {frame} mm です。")
        print(f"\t出力フォルダは {output_dir} です。")
        print("\t原板を" + ("黒色" if invert_color else "白色") + "、星を" + ("白色" if invert_color else "黒色") + "で書き出します。")
    if writer_type == PlateWriterType.SVG:
        return PlateWriterSVG(column, row, frame, True, filename_prefix, invert_color, output_dir)
    if writer_type == PlateWriterType.POSTSCRIPT:
        return PlateWriterPS(column, row, frame, True, filename_prefix, invert_color, output_dir)
    return None


def _init_enlarge_rate(props):
    """星座線を構成する星の拡大率を読み込みます。"""
    if props is None:
        print("星座線を構成する星の拡大率は幾つですか。"
              "(1大きくなるごとに4等星相当の星が1等級明るくなります。default=0.0)")
        enlarge_rate = BasicTransformer.parse_double_with_default(input(), 0.0)
    else:
        # Java版のキーを維持しつつ、Python版向けの表記も受け付ける。
        value = props.get("star.enlarge-rate", props.get("star-EnlargeRate"))
        enlarge_rate = BasicTransformer.parse_double_with_default(value, 0.0)
        print(f"\t星座線を構成する星の拡大率は {enlarge_rate} です。")
    if enlarge_rate < 0:
        raise ValueError("星座線構成星の拡大率には0以上の値を指定してください。")
    return enlarge_rate


def _write_assignment_polygons(transformer, plate_writer, props, force=False):
    """設定に応じて全原盤の担当星域ポリゴンをSVG出力します。"""
    if props is None:
        if force:
            enabled = True
        else:
            print("原盤ごとの担当星域ポリゴンをSVG出力しますか。(y/N) ")
            enabled = input().lower() == "y"
        samples = 180
        filename_prefix = "polygon-"
    else:
        enabled = force or BasicTransformer.parse_boolean_with_default(
            props.get("polygon.enabled", ""), False
        )
        samples = BasicTransformer.parse_int_with_default(
            props.get("polygon.samples"), 180
        )
        filename_prefix = props.get("polygon.file-prefix", "polygon-")
    if not enabled:
        return

    print(f"担当星域ポリゴンを計算しています（1原盤あたり{samples}頂点）。")
    generator = PlateAssignmentPolygonGenerator(transformer, samples=samples)
    polygon_writer = PlateWriterSVG(
        plate_writer.column,
        plate_writer.row,
        plate_writer.r,
        plate_writer.shape,
        filename_prefix,
        True,  # ポリゴン原盤は常に黒地に白い投影領域で出力する。
        plate_writer.output_dir,
    )
    try:
        count = 0
        for dir_, index, points in generator.generate_all(plate_writer.r):
            polygon_writer.write_assignment_polygon(dir_, index, points)
            count += 1
    finally:
        polygon_writer.close()
    print(f"{count}枚の原盤の担当星域ポリゴンをSVG出力しました。")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    BasicTransformer.print_version_and_license()
    writer_type = PlateWriterType.SVG
    config_file_name = None
    force_polygons = False
    # コマンドライン引数を解釈
    i = 0
    while i < len(argv):
        if argv[i].lower() == "-ps":
            writer_type = PlateWriterType.POSTSCRIPT
        elif argv[i] == "--polygons":
            force_polygons = True
        elif argv[i] == "-f":
            if i + 1 >= len(argv):
                print("入力ファイルが指定されていません。", file=sys.stderr)
                return
            if config_file_name is not None:
                print("入力ファイルが2回以上指定されました。2番目以降は無視されます。", file=sys.stderr)
            else:
                config_file_name = argv[i + 1]
            i += 1
        elif argv[i] in ("-h", "-help", "--help"):
            _usage()
            return
        i += 1
    props = load_properties(config_file_name) if config_file_name is not None else None
    t = _init_transformer(props)
    if t is None:
        print("エラーが発生したようなのでOTLを終了します。", file=sys.stderr)
        return
    r = _init_sphere_reader(props)
    w = _init_plate_writer(props, writer_type)
    enlarge_rate = _init_enlarge_rate(props)
    if props is not None and props.get("unit-position-file") is not None:
        print("ユニットの配置を画像に出力しています。")
        unit_position_file = resolve_output_path(_get_output_dir(props), props["unit-position-file"])
        t.write_unit_position(500).save(unit_position_file)
    _write_assignment_polygons(t, w, props, force_polygons)
    print("ユニットを配置しています。。。")
    t.process_stars(r, w, enlarge_rate)
    print("星座を処理しますか。(y/N) ")
    if input().lower() == "y":
        t.process_constellations(r, w)
    print("データを発行しています。")
    w.close()
    print("完了しました。(`･ω･´) ｼｬｷｰﾝ")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2014 東京大学地文研究会天文部
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
"""天球の星から、天の川専用原板への変換。"""

import math
import sys

from basic_transformer import BasicTransformer
from config import load_properties
from models import PlatePosition
from plate_writer import (
    DEFAULT_OUTPUT_DIR,
    PlateWriterPDF,
    PlateWriterSVG,
    PlateWriterType,
    categorized_output_dir,
)
from sphere_reader import SphereReader


def _get_output_dir(props):
    if props is None:
        return DEFAULT_OUTPUT_DIR
    return props.get("output.directory", DEFAULT_OUTPUT_DIR)


class GalaxyTransformer(BasicTransformer):
    def __init__(self, radius, sphere, projector_horizontal, projector_vertical, plate):
        """
        座標変換のパラメータを設定します。
        :param radius: 7.5等星の半径
        :param sphere: ドームの半径
        :param projector_horizontal: ドームの中心から投影機中心へのベクトルの赤道面に平行な成分
        :param projector_vertical: ドームの中心から投影機中心へのベクトルの赤道面に垂直な成分
        :param plate: 投影機中心と原版中心の間の距離
        """
        super().__init__()
        self.radius = radius
        self.sphere = sphere
        self.projector_vertical = projector_vertical
        self.projector_horizontal = projector_horizontal
        self.plate = plate
        self.base_magnitude = 7.5

    def transform(self, sp):
        """レンズの中心を通る光は直進する性質を利用して、天球上の位置をもっとも適切と判定されたユニットの原板上の位置に変換します。"""
        index = 0 if 0. <= sp.radeg <= 180. else 1
        dir_ = 0 if sp.dedeg >= 0 else 1
        return self.transform_unit(sp, dir_, index)

    def transform_unit(self, sp, dir_, index):
        pp = PlatePosition()
        vector = sp.to_vector(self.sphere)
        pp.index = index
        if index == 0:
            vector.y -= self.projector_horizontal
        else:
            vector.y += self.projector_horizontal
        pp.dir = dir_
        if dir_ == 0:
            vector.z -= self.projector_vertical
        else:
            vector.z += self.projector_vertical
        vector.x, vector.y, vector.z = vector.z, vector.x, vector.y
        pp.xmm = self.plate * math.tan(vector.get_lat())
        pp.ymm = -self.plate * vector.get_lng()
        if dir_ == 0:
            if index == 0:
                pp.ymm += self.plate * math.pi / 4
            else:
                pp.ymm += self.plate * math.pi * 3 / 4
        else:
            if index == 0:
                pp.ymm -= self.plate * math.pi / 4
            else:
                pp.ymm -= self.plate * math.pi * 3 / 4
        if pp.ymm < -self.plate * math.pi:
            pp.ymm += self.plate * math.pi * 2
        if pp.ymm > self.plate * math.pi:
            pp.ymm -= self.plate * math.pi * 2
        return pp


def _usage():
    """このプログラムの使い方を表示します。"""
    print("Usage:")
    print("\tpython galaxy_transformer.py [options]")
    print("Options:")
    print("\t-PDF\t原板データを印刷用PDF形式で出力します。")
    print("\t-f [configFile]\tプラネタリウムの設定を configFile から読み込みます。")
    print("\t-h,-help\tこのメッセージを表示します。")


def _init_galaxy_transformer(props):
    if props is None:  # interactive mode
        print("後で何分の一の倍率にしますか。(default=1) ")
        scale = BasicTransformer.parse_double_with_default(input(), 1.)
        print(f"7.5等星の半径は何mmですか。{scale}分の一の倍率をかける前の値を指定してください。(default=0.125) ")
        radius = BasicTransformer.parse_double_with_default(input(), 0.125)
        print("ドームの半径は何mmですか。(default=5000)")
        sphere = BasicTransformer.parse_double_with_default(input(), 5000.)
        print("ドームの中心から投影機中心へのベクトルの赤道面に平行な成分は何mmですか。(default=300)")
        projector_horizontal = BasicTransformer.parse_double_with_default(input(), 300.)
        print("ドームの中心から投影機中心へのベクトルの赤道面に垂直な成分は何mmですか。(default=200)")
        projector_vertical = BasicTransformer.parse_double_with_default(input(), 200.)
        print("投影機の中心と原板の中心の間の距離は何mmですか。(default=50)。")
        plate = BasicTransformer.parse_double_with_default(input(), 50.)
    else:  # non-interactive mode
        scale = BasicTransformer.parse_double_with_default(props.get("scale"), 1.)
        radius = BasicTransformer.parse_double_with_default(props.get("star-radius-7.5"), 0.125)
        sphere = BasicTransformer.parse_double_with_default(props.get("dome-radius"), 5000.)
        projector_horizontal = BasicTransformer.parse_double_with_default(props.get("projector-horizontal"), 300.)
        projector_vertical = BasicTransformer.parse_double_with_default(props.get("projector-vertical"), 200.)
        plate = BasicTransformer.parse_double_with_default(props.get("plate-distance"), 50.)
        print("設定を読み込みました:")
        print(f"\t後で {scale} 分の一の倍率にします。")
        print(f"\t7.5等星の半径は {radius} mm です。")
        print(f"\tドームの半径は {sphere} mm です。")
        print(f"\tドームの中心から投影機中心へのベクトルの赤道面に平行な成分は {projector_horizontal} mm です。")
        print(f"\tドームの中心から投影機中心へのベクトルの赤道面に垂直な成分は {projector_vertical} mm です。")
        print(f"\t投影機の中心と原板の中心の間の距離は {plate} mm です。")
    return GalaxyTransformer(radius, sphere * scale, projector_horizontal * scale, projector_vertical * scale, plate * scale)


def _init_sphere_reader(props):
    if props is None:  # interactive mode
        print("最輝星は何等星ですか。(default=7.5) ")
        maximum = BasicTransformer.parse_double_with_default(input(), 7.5)
        print("最微星は何等星ですか。(default=10.) ")
        minimum = BasicTransformer.parse_double_with_default(input(), 10.)
    else:  # non-interactive mode
        maximum = BasicTransformer.parse_double_with_default(props.get("galaxy.star.maximum"), 7.5)
        minimum = BasicTransformer.parse_double_with_default(props.get("galaxy.star.minimum"), 10.)
        print(f"\t最輝星は {maximum} 等星です。")
        print(f"\t最微星は {minimum} 等星です。")
    return SphereReader(False, maximum, minimum, True)


def _init_plate_writer(props, writer_type):
    filename_prefix = "galaxy-"
    invert_color = True
    output_dir = categorized_output_dir(_get_output_dir(props), "galaxy")
    if props is None:  # interactive mode
        print("横に各天のユニットをいくつ配置しますか。(default=1) ")
        column = BasicTransformer.parse_int_with_default(input(), 1)
        print("縦に各天のユニットをいくつ配置しますか。(default=1) ")
        row = BasicTransformer.parse_int_with_default(input(), 1)
    else:  # non-interactive mode
        column = BasicTransformer.parse_int_with_default(props.get("plate.column"), 1)
        row = BasicTransformer.parse_int_with_default(props.get("plate.row"), 1)
        filename_prefix = props.get("galaxy.file.prefix", filename_prefix)
        invert_color = BasicTransformer.parse_boolean_with_default(props.get("color.invert", ""), False)
        print(f"\t横に各天のユニットを {column} 個配置します。")
        print(f"\t縦に各天のユニットを {row} 個配置します。")
        print(f"\t出力フォルダは {output_dir} です。")
        print("\t原板を" + ("黒色" if invert_color else "白色") + "、星を" + ("白色" if invert_color else "黒色") + "で書き出します。")
    if writer_type == PlateWriterType.SVG:
        return PlateWriterSVG(column, row, 0., False, filename_prefix, invert_color, output_dir)
    if writer_type == PlateWriterType.PDF:
        return PlateWriterPDF(column, row, 0., False, filename_prefix, invert_color, output_dir)
    return None


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    BasicTransformer.print_version_and_license()
    writer_type = PlateWriterType.SVG
    config_file_name = None
    i = 0
    while i < len(argv):
        if argv[i].lower() in ("-pdf", "--pdf"):
            writer_type = PlateWriterType.PDF
        elif argv[i].lower() == "-ps":
            print("PostScript出力は廃止されました。-PDFを使用してください。", file=sys.stderr)
            return 2
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
    t = _init_galaxy_transformer(props)
    r = _init_sphere_reader(props)
    w = _init_plate_writer(props, writer_type)
    print("ユニットを配置しています。。。")
    t.process_stars(r, w)
    print("星座を処理しますか。(y/N) ")
    if input().lower() == "y":
        t.process_constellations(r, w)
    print("データを発行しています。")
    w.close()
    print("完了しました。(`･ω･´) ｼｬｷｰﾝ")


if __name__ == "__main__":
    raise SystemExit(main())

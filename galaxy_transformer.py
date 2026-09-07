# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2014 東京大学地文研究会天文部
# Copyright (C) 2026 東京大学地文研究会天文部
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
import json
import os
from pathlib import Path
import sys

from basic_transformer import BasicTransformer
from config import load_properties
from mathvector import MathVector
from models import PlateConstellation, PlatePosition, SpherePosition
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
    # (N0, N1, S0, S1) の順。N/Sは出力用紙上の位置を表し、銀緯の正負ではない。
    LONGITUDE_REGIONS = ((0., 60.), (60., 120.), (300., 360.), (180., 240.))
    LATITUDE_LIMIT = 20.
    _ANGLE_EPS = 1e-10  # 度。座標回転の丸め誤差を吸収する。

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

        # 投影機の位置は従来どおり赤道座標の(0, ±水平距離, ±垂直距離)。
        # 各担当領域の銀経中央・銀緯0度への光線を原盤中心に合わせる。
        galactic_north = SpherePosition.from_galactic(0., 90.).to_vector(1.)
        self._plate_geometry = {}
        for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
            dir_, index = divmod(region, 2)
            projector = MathVector(
                0., projector_horizontal * (1 if index == 0 else -1),
                projector_vertical * (1 if dir_ == 0 else -1),
            )
            center = SpherePosition.from_galactic((start + end) / 2., 0.)
            forward = center.to_vector(sphere).minus(projector).unit_vector()
            east = galactic_north.cross(forward).unit_vector()
            north = forward.cross(east).unit_vector()
            self._plate_geometry[dir_, index] = (projector, forward, east, north)

    @classmethod
    def _normalize_longitude(cls, longitude):
        longitude %= 360.
        for boundary in (0., 60., 120., 180., 240., 300., 360.):
            if abs(longitude - boundary) <= cls._ANGLE_EPS:
                return boundary % 360.
        return longitude

    def assigned_unit(self, sp):
        """銀緯±20度の担当 ``(出力位置, 番号)`` を返す。対象外はNone。"""
        longitude, latitude = sp.to_galactic()
        if not abs(latitude) <= self.LATITUDE_LIMIT + self._ANGLE_EPS:
            return None
        longitude = self._normalize_longitude(longitude)
        for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
            if start <= longitude < end:
                return divmod(region, 2)
        return None

    def transform(self, sp):
        """銀河座標で担当原盤を選び、領域外の位置は除外します。"""
        unit = self.assigned_unit(sp)
        if unit is None:
            return None
        return self.transform_unit(sp, *unit)

    def transform_unit(self, sp, dir_, index):
        """担当領域の中央を基準に、従来の円筒投影で原盤座標へ変換します。"""
        projector, forward, east, north = self._plate_geometry[dir_, index]
        ray = sp.to_vector(self.sphere).minus(projector)
        x, y, z = ray.dot(forward), ray.dot(east), ray.dot(north)
        pp = PlatePosition()
        pp.index = index
        pp.dir = dir_
        pp.xmm = self.plate * z / math.hypot(x, y)
        pp.ymm = -self.plate * math.atan2(y, x)
        return pp

    def inverse_transform_unit(self, xmm, ymm, dir_, index):
        """原盤上の点から、同じ円筒投影の光線とドームの交点を求める。"""
        projector, forward, east, north = self._plate_geometry[dir_, index]
        angle = -ymm / self.plate
        ray = (forward.mult_scalar(math.cos(angle))
               .plus(east.mult_scalar(math.sin(angle)))
               .plus(north.mult_scalar(xmm / self.plate))).unit_vector()
        dot = projector.dot(ray)
        discriminant = dot * dot + self.sphere ** 2 - projector.dot(projector)
        if discriminant < 0:
            raise ValueError("原盤の光線がドームと交差しません。")
        distance = -dot + math.sqrt(discriminant)
        return SpherePosition.from_vector(projector.plus(ray.mult_scalar(distance)))

    @staticmethod
    def _clip_line(longitude, latitude, delta_l, delta_b, start, end, limit):
        """銀経・銀緯で線分を長方形領域へクリップし、線分上の区間を返す。"""
        lower, upper = 0., 1.
        for origin, delta, minimum, maximum in (
            (longitude, delta_l, start, end),
            (latitude, delta_b, -limit, limit),
        ):
            if abs(delta) <= GalaxyTransformer._ANGLE_EPS:
                if (origin < minimum - GalaxyTransformer._ANGLE_EPS
                        or origin > maximum + GalaxyTransformer._ANGLE_EPS):
                    return None
                continue
            entry, leave = sorted(((minimum - origin) / delta, (maximum - origin) / delta))
            lower, upper = max(lower, entry), min(upper, leave)
            if lower >= upper:
                return None
        return lower, upper

    def transform_constellation(self, sc):
        """星座名を領域で選別し、星座線は銀経・銀緯の境界で分割します。"""
        pc = PlateConstellation()
        pc.hip_numbers = list(sc.hip_numbers)
        if sc.name is not None and sc.p is not None:
            pc.name, pc.p = sc.name, self.transform(sc.p)
        pc.ll = []
        for first, last in sc.ll or ():
            longitude, latitude = first.to_galactic()
            end_l, end_b = last.to_galactic()
            longitude = self._normalize_longitude(longitude)
            delta_l = (self._normalize_longitude(end_l) - longitude + 180.) % 360. - 180.
            delta_b = end_b - latitude
            for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
                dir_, index = divmod(region, 2)
                # 0/360度をまたぐ線分も、それぞれの原盤の境界で分割する。
                for shift in (-360., 0., 360.):
                    left, right = start + shift, end + shift
                    if abs(delta_l) <= self._ANGLE_EPS and not left <= longitude < right:
                        continue
                    interval = self._clip_line(
                        longitude, latitude, delta_l, delta_b,
                        left, right, self.LATITUDE_LIMIT,
                    )
                    if interval is None:
                        continue
                    pc.ll.append([
                        self.transform_unit(SpherePosition.from_galactic(
                            longitude + t * delta_l, latitude + t * delta_b
                        ), dir_, index)
                        for t in interval
                    ])
        return pc if pc.p is not None or pc.ll else None


def _usage():
    """このプログラムの使い方を表示します。"""
    print("Usage:")
    print("\tpython galaxy_transformer.py [options]")
    print("Options:")
    print("\t-PDF\t原板データを印刷用PDF形式で出力します。")
    print("\t-f [configFile]\tプラネタリウムの設定を configFile から読み込みます。")
    print("\t--download-gaia\tGaia DR3の光量マップを取得・キャッシュして原盤を生成します。")
    print("\t--gaia-query\t取得用ADQLクエリを表示して終了します（通信なし）。")
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
    print(f"\t銀緯±{GalaxyTransformer.LATITUDE_LIMIT:g}度を4枚の原盤に割り当てます。")
    for region, (start, end) in enumerate(GalaxyTransformer.LONGITUDE_REGIONS):
        dir_, index = divmod(region, 2)
        label = ("N" if dir_ == 0 else "S") + str(index)
        print(f"\t{label}: 銀経{start:g}度以上{end:g}度未満")
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


def _init_plate_writer(props, writer_type, write_frames=True):
    filename_prefix = "galaxy-"
    invert_color = True
    output_dir = categorized_output_dir(_get_output_dir(props), "galaxy")
    if props is None:  # interactive mode
        print("横に各段の原盤をいくつ配置しますか。(default=1) ")
        column = BasicTransformer.parse_int_with_default(input(), 1)
        print("縦に各段の原盤をいくつ配置しますか。(default=1) ")
        row = BasicTransformer.parse_int_with_default(input(), 1)
    else:  # non-interactive mode
        column = BasicTransformer.parse_int_with_default(props.get("plate.column"), 1)
        row = BasicTransformer.parse_int_with_default(props.get("plate.row"), 1)
        filename_prefix = props.get("galaxy.file.prefix", filename_prefix)
        invert_color = BasicTransformer.parse_boolean_with_default(props.get("color.invert", ""), False)
        print(f"\t横に各段の原盤を {column} 個配置します。")
        print(f"\t縦に各段の原盤を {row} 個配置します。")
        print(f"\t出力フォルダは {output_dir} です。")
        print("\t原板を" + ("黒色" if invert_color else "白色") + "、星を" + ("白色" if invert_color else "黒色") + "で書き出します。")
    if writer_type == PlateWriterType.SVG:
        writer = PlateWriterSVG(column, row, 0., False, filename_prefix, invert_color, output_dir)
    elif writer_type == PlateWriterType.PDF:
        writer = PlateWriterPDF(column, row, 0., False, filename_prefix, invert_color, output_dir)
    else:
        return None
    if write_frames:
        writer.write_frames(2)
    return writer


def _process_gaia(props, config_directory, writer_type, fetch=False, query_only=False):
    from gaia_catalog import (GaiaSettings, build_queries, download_gaia,
                              iter_flux_samples, positive_number)
    from galaxy_etching import EtchingBuilder, EtchingSettings

    gaia = GaiaSettings.from_properties(props, config_directory)
    etching = EtchingSettings.from_properties(props)
    regions, limit = GalaxyTransformer.LONGITUDE_REGIONS, GalaxyTransformer.LATITUDE_LIMIT
    if query_only:
        for index, query in enumerate(build_queries(gaia, regions, limit)):
            print(f"-- region {index}\n{query};\n")
        return 0
    positive_number(props, "star-radius-7.5", .125)
    sphere = positive_number(props, "dome-radius", 5000.)
    positive_number(props, "plate-distance", 50.)
    horizontal = float(props.get("projector-horizontal", 300.))
    vertical = float(props.get("projector-vertical", 200.))
    if (not math.isfinite(horizontal) or not math.isfinite(vertical)
            or math.hypot(horizontal, vertical) >= sphere):
        raise ValueError("Gaiaモードの半径・距離は有限の正値、投影機位置はドームの内側が必要です。")
    transformer = _init_galaxy_transformer(props)
    if fetch:
        download_gaia(gaia, regions, limit)
    output_props = dict(props)
    output_props["galaxy.file.prefix"] = props.get("galaxy.gaia.file.prefix", "gaia-etch-")
    # 読み込み・配置に失敗したとき、既存原盤を途中結果で上書きしない。
    writer = _init_plate_writer(output_props, writer_type, write_frames=False)
    try:
        builder = EtchingBuilder(transformer, etching, writer.r)
        input_stats = {}
        print(f"Gaia G={gaia.bright:g}より暗く{gaia.faint:g}以下の光量を集計します。", flush=True)
        print(f"最終穴径 {etching.min_diameter_mm:g} mm、金属幅 {etching.min_web_mm:g} mm（仮の加工条件）", flush=True)
        for position, flux, count in iter_flux_samples(gaia, regions, limit, input_stats):
            builder.add_sample(position, flux, count)
            if input_stats["rows"] % 100000 == 0:
                print(f"  {input_stats['rows']:,}レコード / {input_stats['stars']:,}天体", flush=True)
        if not sum(builder.counts.values()):
            raise ValueError("Gaia入力に担当領域・等級範囲の天体がありません。")
        print("最小穴径・金属幅を満たす穴を配置しています。", flush=True)
        holes, report = builder.build()
        if not holes:
            raise ValueError("出力できる穴がありません。光量倍率・最小穴径・用紙配置を確認してください。")
        report["input"] = input_stats
        report["geometry"] = {k: getattr(transformer, k) for k in (
            "radius", "sphere", "projector_horizontal", "projector_vertical", "plate")}
        writer.write_frames(2)
        for hole in holes:
            writer.write_star(hole)
    finally:
        writer.close()
    report_path = Path(writer.output_dir) / f"{writer.filename_prefix}report.json"
    partial = report_path.with_suffix(".json.part")
    partial.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, report_path)
    for item in report["plates"]:
        requested = item["requested_area_mm2"]
        loss = item["capacity_loss_area_mm2"] / requested if requested else 0.
        print(f"{item['plate']}: {item['source_stars']:,}天体 → {item['holes']:,}穴、容量不足の光量 {loss:.2%}")
        if loss > .01:
            print("  注意: 濃淡が飽和しています。galaxy.etch.flux-gainを下げるか、加工条件・集計セルを見直してください。")
    print(f"生成完了。加工条件と光量の記録: {report_path}")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    BasicTransformer.print_version_and_license()
    writer_type = PlateWriterType.SVG
    config_file_name = None
    fetch_gaia = query_gaia = False
    i = 0
    while i < len(argv):
        if argv[i].lower() in ("-pdf", "--pdf"):
            writer_type = PlateWriterType.PDF
        elif argv[i].lower() == "-ps":
            print("PostScript出力は廃止されました。-PDFを使用してください。", file=sys.stderr)
            return 2
        elif argv[i] == "--download-gaia":
            fetch_gaia = True
        elif argv[i] == "--gaia-query":
            query_gaia = True
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
    mode = (props or {}).get("galaxy.mode", "legacy").strip().lower()
    if mode not in {"legacy", "gaia-etch"}:
        print("galaxy.mode は legacy または gaia-etch を指定してください。", file=sys.stderr)
        return 2
    if mode == "gaia-etch" or fetch_gaia or query_gaia:
        if props is None:
            print("Gaiaモードでは -f で設定ファイルを指定してください。", file=sys.stderr)
            return 2
        try:
            return _process_gaia(props, Path(config_file_name).resolve().parent,
                                 writer_type, fetch_gaia, query_gaia)
        except (OSError, ValueError, RuntimeError) as exc:
            print(f"Gaia原盤の生成に失敗しました: {exc}", file=sys.stderr)
            return 2
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

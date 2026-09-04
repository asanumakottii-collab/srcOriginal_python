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
"""hip_main.dat / tyc_main.dat / rc3.dat / IAU88.hlc を読み込むモジュールです。

NOTE: 乱数生成には Python 標準の random モジュールを使用しています。java.util.Random とはアルゴリズムが異なるため、同じシード値(0)でも生成される乱数列そのものは Java 版と一致しません(統計的性質は同等です)。
"""

import math
import os
import random

from mathvector import MathVector
from models import SphereConstellation, SpherePosition, SphereStar

# hip_main.dat 等のカタログはこのモジュールと同じフォルダに置く。実行時のカレントディレクトリに関わらず python フォルダ単体で完結させるため、カレントディレクトリではなくこのファイルの場所を基準に探す。
_DATA_DIR = os.path.dirname(os.path.abspath(__file__))

# OTL_prog の hip_constellation_line_star.csv にだけ含まれ、IAU88.hlc のHIPブロックには無い星。補助CSVを同梱しない配布形態でも旧版と同じ対象を拡大できるよう、HIP番号だけを保持する。
_LEGACY_CONSTELLATION_STAR_HIP_NUMBERS = {
    2072, 3760, 3881, 5348, 7083, 8645, 8833, 8837, 10324, 10559, 10826, 11345,
    11783, 12093, 12390, 12484, 13209, 13254, 17448, 17847, 18246, 18505, 18614,
    20455, 21060, 21594, 21881, 21949, 22549, 22845, 23123, 23453, 24244, 24845,
    24873, 25110, 25918, 26207, 27288, 27890, 28103, 28691, 28910, 29151, 30883,
    31416, 31592, 32362, 32759, 32768, 33165, 33347, 33449, 33977, 34769, 35037,
    35264, 36046, 36145, 36962, 37740, 40843, 42402, 42536, 42568, 43234, 44382,
    44700, 46509, 46776, 48402, 49641, 50335, 51232, 51437, 52468, 53253, 54463,
    56211, 56480, 56561, 57363, 57380, 59196, 59199, 60030, 60823, 61932, 63125,
    64166, 65109, 65936, 66249, 68520, 69427, 69701, 70576, 71683, 71795, 71957,
    72220, 72622, 73334, 75141, 76333, 76552, 77760, 77853, 81126, 82514, 82671,
    84606, 85112, 85267, 85755, 85829, 87072, 88866, 90595, 90887, 91875, 91971,
    92202, 92814, 92946, 92953, 92989, 93085, 93244, 93683, 94160, 94779, 94820,
    95168, 95294, 95347, 96406, 97649, 97804, 98032, 98110, 98412, 98543, 98688,
    98920, 100345, 102618, 105858, 107310, 107608, 108661, 109352, 111123, 112440,
    112447, 112716, 113136, 113638, 114131, 114855, 115438,
}


def _readline(f):
    """java.util.BufferedReader.readLine() 相当。EOFではNoneを返す。"""
    line = f.readline()
    if line == "":
        return None
    return line.rstrip("\r\n")


class SphereReader:
    def __init__(self, above_maximum, maximum, minimum, under_minimum):
        self._in_hip = open(os.path.join(_DATA_DIR, "hip_main.dat"), encoding="ascii")
        self._in_tyc = open(os.path.join(_DATA_DIR, "tyc_main.dat"), encoding="ascii")
        self._in_rc3 = open(os.path.join(_DATA_DIR, "rc3.dat"), encoding="ascii")
        self._in_hlc = open(os.path.join(_DATA_DIR, "IAU88.hlc"), encoding="utf-8")
        self._above_maximum = above_maximum
        self._maximum = maximum
        self._minimum = minimum
        self._under_minimum = under_minimum
        self._random = random.Random(0)
        self._galaxy_count = 0
        self._galaxy_scale = 0.0
        self._galaxy_vectors = [None, None, None]
        self._excluding_stars = set()

    def exclude_star(self, hip_number):
        self._excluding_stars.add(str(hip_number))

    @staticmethod
    def _constellation_hip_numbers_by_name():
        """IAU88.hlc の星座名とHIP番号の対応を読み込みます。"""
        hip_numbers_by_name = {}
        current_name = None
        in_hip_block = False
        with open(os.path.join(_DATA_DIR, "IAU88.hlc"), encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.strip()
                if line.startswith("Name:"):
                    current_name = line.split(":", 1)[1].strip()
                    hip_numbers_by_name.setdefault(current_name, set())
                    in_hip_block = False
                elif line == "HIP: {":
                    in_hip_block = True
                elif in_hip_block and line == "}":
                    in_hip_block = False
                elif in_hip_block and line and line != "-1":
                    if current_name is None:
                        raise ValueError("IAU88.hlcのHIPブロックに星座名がありません。")
                    hip_numbers_by_name[current_name].add(int(line))
        return hip_numbers_by_name

    @staticmethod
    def constellation_names():
        """IAU88.hlc に記載された ``Name`` をファイル順で返します。"""
        return tuple(SphereReader._constellation_hip_numbers_by_name())

    @staticmethod
    def constellation_star_hip_numbers(constellation_names=None):
        """星座線を構成する恒星のHIP番号を返します。

        ``constellation_names`` に IAU88.hlc の ``Name`` を指定すると、
        指定星座のHIP番号だけを返します。``None`` の場合は従来どおり
        全星座を対象とし、Java版の補助CSVにだけ含まれる星も併合します。
        補助CSVには星座名がないため、星座を指定した場合は併合しません。
        """
        hip_numbers_by_name = SphereReader._constellation_hip_numbers_by_name()

        if constellation_names is not None:
            if isinstance(constellation_names, str):
                constellation_names = [constellation_names]
            selected_names = tuple(dict.fromkeys(constellation_names))
            unknown_names = [
                name for name in selected_names if name not in hip_numbers_by_name
            ]
            if unknown_names:
                raise ValueError(
                    "IAU88.hlcに存在しない星座名です: "
                    + ", ".join(unknown_names)
                )
            hip_numbers = set()
            for name in selected_names:
                hip_numbers.update(hip_numbers_by_name[name])
            return hip_numbers

        hip_numbers = set(_LEGACY_CONSTELLATION_STAR_HIP_NUMBERS)
        for constellation_hip_numbers in hip_numbers_by_name.values():
            hip_numbers.update(constellation_hip_numbers)

        csv_path = os.path.join(_DATA_DIR, "hip_constellation_line_star.csv")
        if os.path.exists(csv_path):
            with open(csv_path, encoding="utf-8") as f:
                for raw_line in f:
                    fields = raw_line.split(",", 1)
                    if fields and fields[0].strip():
                        hip_numbers.add(int(fields[0]))
        return hip_numbers

    def read_star(self):
        s = SphereStar()
        s.p = SpherePosition()
        while True:
            line = _readline(self._in_tyc)
            if line is not None:
                # Hipparcos HIP number: HIP番号が振られている星は無視する
                value = line[210:216]
                if value != "      ":
                    continue

                # Magnitude in Johnson V
                value = line[41:46]
                if value == "     ":
                    continue
                s.vmag = float(value)
                if s.vmag <= self._maximum:
                    if not self._above_maximum:
                        continue
                    s.vmag = self._maximum
                elif s.vmag > self._minimum:
                    if (not self._under_minimum
                            or s.vmag > self._minimum - math.log(self._random.random()) / math.log(2.512)):
                        continue
                    s.vmag = self._minimum
                # alpha, degrees (ICRS, Epoch=J1991.25)
                value = line[51:63]
                if value == "            ":
                    continue
                s.p.radeg = float(value)

                # delta, degrees (ICRS, Epoch=J1991.25)
                value = line[64:76]
                if value == "            ":
                    continue
                s.p.dedeg = float(value)
                return s

            line = _readline(self._in_hip)
            if line is not None:
                # HIP number
                value = line[8:14]
                if value.strip() in self._excluding_stars:
                    continue
                s.hip_number = int(value)

                # Magnitude in Johnson V
                value = line[41:46]
                if value == "     ":
                    continue
                s.vmag = float(value)
                if s.vmag <= self._maximum:
                    if not self._above_maximum:
                        continue
                    s.vmag = self._maximum
                elif s.vmag > self._minimum:
                    if (not self._under_minimum
                            or s.vmag > self._minimum - math.log(self._random.random()) / math.log(2.512)):
                        continue
                    s.vmag = self._minimum
                # alpha, degrees (ICRS, Epoch=J1991.25)
                value = line[51:63]
                if value == "            ":
                    continue
                s.p.radeg = float(value)

                # delta, degrees (ICRS, Epoch=J1991.25)
                value = line[64:76]
                if value == "            ":
                    continue
                s.p.dedeg = float(value)
                return s

            if self._galaxy_count != 0:
                self._galaxy_count -= 1
                r = math.exp(self._random.gauss(0.0, 1.0)) / self._galaxy_scale
                theta = self._random.random() * math.pi * 2
                vector = (self._galaxy_vectors[2]
                          .plus(self._galaxy_vectors[0].mult_scalar(r * math.cos(theta)))
                          .plus(self._galaxy_vectors[1].mult_scalar(r * math.sin(theta))))
                s.vmag = self._minimum
                s.p = SpherePosition.from_vector(vector)
                return s

            line = _readline(self._in_rc3) if self._under_minimum else None
            if line is not None:
                # BT (total B magnitude)
                value = line[189:194]
                if value == "     ":
                    continue
                btmag = float(value)

                # (B-V)T (total (B-V))
                value = line[252:256]
                if value == "    ":
                    continue
                vtmag = btmag - float(value)

                # Right Ascension B2000 (hours)
                value = line[0:2]
                if value == "  ":
                    continue
                radeg = float(value) * 15

                # Right Ascension B2000 (minutes)
                value = line[2:4]
                if value == "  ":
                    continue
                radeg += float(value) / 4

                # Right Ascension B2000 (sec. or min.)
                value = line[4:8]
                if value == "    ":
                    continue
                radeg += float(value) / 240

                # Declination B2000 (degrees)
                value = line[10:12]
                if value == "   ":
                    continue
                dedeg = float(value)

                # Declination B2000 (minutes)
                value = line[12:14]
                if value == "  ":
                    continue
                dedeg += float(value) / 60

                # Declination B2000 (seconds)
                value = line[14:16]
                if value == "  ":
                    continue
                dedeg += float(value) / 3600

                # Sign of declination
                if line[9] == "-":
                    dedeg *= -1

                # Log D25
                value = line[151:155]
                if value == "    ":
                    continue
                d25 = 10 ** float(value)

                # Mean error on log D25
                value = line[157:160]
                if value == "   ":
                    continue
                r25 = 10 ** float(value)

                # Position angle of the major axis
                value = line[185:188]
                if value == "   ":
                    continue
                pa = float(value)

                galaxy_count = int(math.pow(2.512, self._minimum - vtmag))
                if galaxy_count <= 0:
                    continue

                scale_argument = (0.5 * math.pi * math.sqrt(2. * math.pi) * d25 * d25 / r25
                                  * math.pow(2.512, -25 + vtmag) * 36. * 0.25)
                # この分布モデルでは scale_argument が (0, 1] の範囲にある必要がある。
                # RC3 にはこの条件を満たさないレコードがあり、Python では負数の平方根が
                # ValueError になるため、生成不能な銀河を読み飛ばす。
                if not 0. < scale_argument <= 1.:
                    continue
                self._galaxy_scale = math.exp(math.sqrt(-2. * math.log(scale_argument)))
                gv = self._galaxy_vectors
                gv[2] = MathVector.from_mag_lng_lat(1., math.radians(radeg), math.radians(dedeg))
                gv[0] = MathVector(0, 0, 1.)
                gv[1] = gv[2].cross(gv[0]).unit_vector()
                gv[0] = gv[1].cross(gv[2])
                gv[0] = gv[0].mult_scalar(math.cos(math.radians(pa))).plus(gv[1].mult_scalar(math.sin(-math.radians(pa))))
                gv[1] = gv[2].cross(gv[0])
                gv[0] = gv[0].mult_scalar(math.radians(d25 / 600 / 2))
                gv[1] = gv[1].mult_scalar(math.radians(d25 / 600 / 2 / r25))
                self._galaxy_count = galaxy_count
                continue

            return None

    def read_constellation(self):
        c = SphereConstellation()
        c.p = SpherePosition()
        line = _readline(self._in_hlc)
        while line is not None and line == "":
            line = _readline(self._in_hlc)
        if line is None:
            return None
        while True:
            separator = line.find(":")
            if separator != -1:
                key = line[:separator].strip()
                value = line[separator + 1:].strip()
                if key == "Name":
                    c.name = value
                elif key == "Pos":
                    separator = value.find(" ")
                    if separator != -1:
                        c.p.radeg = float(value[:separator])
                        c.p.dedeg = float(value[separator + 1:])
                elif key == "Lines":
                    p0 = None
                    p1 = None
                    ll = []
                    while True:
                        line = _readline(self._in_hlc)
                        if line is None:
                            break
                        line = line.strip()
                        if line == "}":
                            break
                        p0 = p1
                        p1 = None
                        separator = line.find(" ")
                        if separator == -1:
                            continue
                        y = float(line[:separator].strip())
                        line = line[separator + 1:].strip()
                        separator = line.find(" ")
                        if separator == -1:
                            continue
                        z = float(line[:separator].strip())
                        line = line[separator + 1:].strip()
                        separator = line.find(" ")
                        if separator == -1:
                            continue
                        x = float(line[:separator].strip())
                        p1 = SpherePosition.from_vector(MathVector(x, y, z))
                        if p0 is not None and p0 != p1:
                            ll.append([p0, p1])
                    c.ll = ll
                elif key == "HIP":
                    while True:
                        line = _readline(self._in_hlc)
                        if line is None or line.strip() == "}":
                            break
                        value = line.strip()
                        if value and value != "-1":
                            c.hip_numbers.append(int(value))
                elif value == "{":
                    while True:
                        line = _readline(self._in_hlc)
                        if line is None or line.strip() == "}":
                            break
            # do-whileの条件部に相当(continue時もこの行に到達する)
            line = _readline(self._in_hlc)
            if line is None or line == "":
                break
        return c

    def close(self):
        self._in_hip.close()
        self._in_tyc.close()
        self._in_rc3.close()
        self._in_hlc.close()

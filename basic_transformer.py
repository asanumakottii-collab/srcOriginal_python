# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012-2014 東京大学地文研究会天文部
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
"""Transformer と GalaxyTransformer に共通する処理をまとめた抽象基底クラス。"""

import math
from abc import ABC, abstractmethod

from models import PlateConstellation, PlateStar


class BasicTransformer(ABC):
    CONSTELLATION_TAPER_BOUNDARY_MAGNITUDE = 4.0
    CONSTELLATION_TAPER_WIDTH = 2.0
    MAX_CONSTELLATION_ENLARGE_RATE = 2.0

    def __init__(self):
        self.radius = 0.0  # baseMagnitude等星の原板上の半径
        self.base_magnitude = 0.0  # 原板上の穴の大きさがradiusとなるような明るさ

    @abstractmethod
    def transform(self, sp):
        """天球上の位置をもっとも適切と判定されたユニットの原板上の位置に変換します。"""
        raise NotImplementedError

    @abstractmethod
    def transform_unit(self, sp, dir_, index):
        """天球上の位置を指定されたユニットの原板上の位置に変換します。"""
        raise NotImplementedError

    def transform_star(self, ss):
        """天球上の星を原板上の星に変換します。"""
        ps = PlateStar()
        ps.rmm = self.radius * math.pow(math.sqrt(2.512), -ss.vmag + self.base_magnitude)
        ps.p = self.transform(ss.p)
        ps.hip_number = ss.hip_number
        return ps

    def transform_constellation(self, sc):
        """天球上の星座を原板上の星座に変換します。"""
        pc = PlateConstellation()
        if sc.name is not None and sc.p is not None:
            pc.name = sc.name
            pc.p = self.transform(sc.p)
            pc.hip_numbers = list(sc.hip_numbers)
        if sc.ll is not None:
            ll = []
            for line in sc.ll:
                l1 = [None, None]
                l1[0] = self.transform(line[0])
                l1[1] = self.transform_unit(line[1], l1[0].dir, l1[0].index)
                ll.append(l1)
                l2 = [None, None]
                l2[0] = self.transform(line[1])
                if not l2[0].is_same_plate(l1[1]):
                    l2[1] = self.transform_unit(line[0], l2[0].dir, l2[0].index)
                    ll.append(l2)
            pc.ll = ll
        return pc

    @staticmethod
    def enlarge_constellation_star(rmm, base_radius, enlarge_rate, magnitude=None):
        """等級に対するガウス型テーパーで星座線構成星の穴半径を拡大します。

        magnitudeを省略した旧形式の呼び出しでは、0等星のbase_radiusから等級を逆算します。
        """
        if not (
            0.0
            <= enlarge_rate
            <= BasicTransformer.MAX_CONSTELLATION_ENLARGE_RATE
        ):
            raise ValueError("星座線構成星の拡大率には0以上2以下の値を指定してください。")
        if enlarge_rate == 0:
            return rmm

        if magnitude is None:
            if rmm <= 0 or base_radius <= 0:
                raise ValueError("星の穴半径と0等星の基準半径には正の値を指定してください。")
            magnitude = -5.0 * math.log10(rmm / base_radius)

        boundary = BasicTransformer.CONSTELLATION_TAPER_BOUNDARY_MAGNITUDE
        if magnitude >= boundary:
            taper = 1.0
        else:
            distance = (boundary - magnitude) / BasicTransformer.CONSTELLATION_TAPER_WIDTH
            taper = math.exp(-(distance * distance))

        effective_magnitude_shift = enlarge_rate * taper
        return rmm * math.pow(10.0, 0.2 * effective_magnitude_shift)

    def process_stars(
        self,
        r,
        w,
        enlarge_rate=0.0,
        constellation_hip_numbers=None,
        constellation_names=None,
    ):
        counter = 0
        enlarged_counter = 0
        if constellation_hip_numbers is None:
            if enlarge_rate == 0:
                constellation_hip_numbers = set()
            elif constellation_names is None:
                constellation_hip_numbers = r.constellation_star_hip_numbers()
            else:
                constellation_hip_numbers = r.constellation_star_hip_numbers(
                    constellation_names
                )
        while True:
            ss = r.read_star()
            if ss is None:
                break
            ps = self.transform_star(ss)
            if ps is None or ps.p is None:
                continue
            if ps.hip_number in constellation_hip_numbers:
                ps.rmm = self.enlarge_constellation_star(
                    ps.rmm, self.radius, enlarge_rate, magnitude=ss.vmag
                )
                enlarged_counter += 1
            if counter % 1000 == 0:
                print(f"{counter}個目の星を変換しています。")
            w.write_star(ps)
            counter += 1
        print(f"{counter}個の星を変換しました。")
        if enlarge_rate != 0:
            print(f"{enlarged_counter}個の星座線構成星を拡大しました。")

    def process_constellations(self, r, w):
        counter = 0
        while True:
            sc = r.read_constellation()
            if sc is None:
                break
            pc = self.transform_constellation(sc)
            if pc is None:
                continue
            if counter % 10 == 0:
                print(f"{counter}個目の星座を変換しています。")
            w.write_constellation(pc)
            counter += 1
        print(f"{counter}個の星座の変換しました。")

    @staticmethod
    def print_version_and_license():
        print("Orb Transform Library(OTL) Version 1.4")
        print("Copyright (C) 2007,2012-2014 東京大学地文研究会天文部")
        print("")
        print("This program is free software; you can redistribute it and/or modify it under "
              "the terms of the GNU General Public License as published by the Free Software "
              "Foundation; either version 2 of the License, or (at your option) any later version.")
        print("")
        print("This program is distributed in the hope that it will be useful, but WITHOUT ANY "
              "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A "
              "PARTICULAR PURPOSE.  See the GNU General Public License for more details.")
        print("")
        print("You should have received a copy of the GNU General Public License along with this "
              "program; if not, write to the Free Software Foundation, Inc., 59 Temple Place - "
              "Suite 330, Boston, MA  02111-1307, USA.")
        print("")

    @staticmethod
    def parse_double_with_default(s, d):
        try:
            return float(s)
        except (TypeError, ValueError):
            return d

    @staticmethod
    def parse_int_with_default(s, d):
        try:
            return int(s)
        except (TypeError, ValueError):
            return d

    @staticmethod
    def parse_boolean_with_default(s, d):
        if s is None:
            return d
        sl = s.lower()
        if sl in ("y", "yes", "true"):
            return True
        if sl in ("n", "no", "false"):
            return False
        return d

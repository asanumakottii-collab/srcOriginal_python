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
"""天球・原板上の位置、恒星、星座を表すデータクラス群。"""

import math

from mathvector import MathVector


class PlatePosition:
    """原板上での位置を扱うクラスです。"""

    def __init__(self):
        self.xmm = 0.0
        self.ymm = 0.0
        self.dir = 0  # 北天(0)，南天(1)の別
        self.index = 0  # ユニット番号

    def __eq__(self, other):
        if not isinstance(other, PlatePosition):
            return NotImplemented
        # NOTE: 原Javaコードの isSamePlate 判定にあった dir と p.index の
        # 比較をそのまま踏襲している(挙動を変えないための忠実な移植)。
        return (self.xmm == other.xmm and self.ymm == other.ymm
                and self.dir == other.index and self.index == other.index)

    def is_same_plate(self, other):
        """同じ原板上にあるかどうか"""
        return self.dir == other.dir and self.index == other.index

    def distance_from_center(self):
        """原板中心からの距離"""
        return math.sqrt(self.xmm * self.xmm + self.ymm * self.ymm)


class SpherePosition:
    """赤道座標を扱うクラスです。"""

    def __init__(self):
        self.dedeg = 0.0  # 赤緯
        self.radeg = 0.0  # 赤経

    def __eq__(self, other):
        if not isinstance(other, SpherePosition):
            return NotImplemented
        return self.dedeg == other.dedeg and self.radeg == other.radeg

    def to_vector(self, radius):
        """直交座標に変換します。"""
        return MathVector.from_mag_lng_lat(radius, math.radians(self.radeg), math.radians(self.dedeg))

    @staticmethod
    def from_vector(vector):
        """直交座標から新しいSpherePositionインスタンスを作成します。"""
        p = SpherePosition()
        p.radeg = math.degrees(vector.get_lng())
        p.dedeg = math.degrees(vector.get_lat())
        return p


class PlateStar:
    """原板上の恒星を扱うクラスです。"""

    def __init__(self):
        self.rmm = 0.0  # 原板上での半径
        self.p = None  # 原板上での位置 (PlatePosition)
        self.hip_number = None  # Hipparcos HIP番号 (Tycho星などはNone)


class SphereStar:
    """天球上の恒星を扱うクラスです。"""

    def __init__(self):
        self.vmag = 0.0  # 恒星のV等級
        self.p = None  # 恒星の赤道座標 (SpherePosition)
        self.hip_number = None  # Hipparcos HIP番号 (Tycho星などはNone)


class PlateConstellation:
    """原板上の星座を扱うクラスです。"""

    def __init__(self):
        self.name = None  # 星座名
        self.p = None  # 原板上での星座の中心位置 (PlatePosition)
        self.ll = None  # 原板上での星座線の端点を格納する配列
        self.hip_numbers = []  # 星座線を構成する恒星のHIP番号


class SphereConstellation:
    """天球上の星座を扱うクラスです。"""

    def __init__(self):
        self.name = None  # 星座名
        self.p = None  # 天球上での星座の中心位置 (SpherePosition)
        self.ll = None  # 天球上での星座線の端点を格納する配列
        self.hip_numbers = []  # 星座線を構成する恒星のHIP番号

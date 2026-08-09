# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2013 東京大学地文研究会天文部
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
"""ユニット(投影機)の配置を扱うモジュールです。"""

import math
from abc import ABC, abstractmethod

from mathvector import MathVector
from geometry import GeneralizedCone, Polyhedron

# 浮動小数点の計算で生じる誤差よりも十分大きな数 ここでは 2^-16
_EPS = 0.0000152587890625


class InvalidUnitArrangementException(Exception):
    """ユニット配置の指定が無効だった時に投げられる例外のクラス"""


class UnitArrangement(ABC):
    """ユニットの配置を表現する抽象クラス"""

    @abstractmethod
    def get_unit_vector(self, i, j):
        """指定された半球上にあるユニットの向きを表すベクトルを返します。"""
        raise NotImplementedError

    @abstractmethod
    def number_of_units(self, i):
        """指定された半球上にあるユニットの個数を返します。"""
        raise NotImplementedError


class StandardUnitArrangement(UnitArrangement):
    """正十二面体に基づいたユニットの配置を表現するクラス"""

    def __init__(self, has_unit_12, has_unit_20, has_unit_30, has_unit_60, has_unit_extra):
        """
        ユニットを配置する方向を計算します。
        :param has_unit_12: ユニットを正十二面体の面方向12箇所に配置する場合True
        :param has_unit_20: ユニットを正十二面体の頂点方向20箇所に配置する場合True
        :param has_unit_30: ユニットを正十二面体の辺方向30箇所に配置する場合True
        :param has_unit_60: ユニットを正十二面体の対称性を持っている最適な方向60箇所に配置する場合True
        :param has_unit_extra: 加えてユニットの赤道付近に配置する場合True
        """
        self._has_unit_12 = has_unit_12
        self._has_unit_20 = has_unit_20
        self._has_unit_30 = has_unit_30
        self._has_unit_60 = has_unit_60
        self._has_unit_extra = has_unit_extra

        # 正二十面体の対称性を持っている凧型六十面体の面を対称軸で分割した三角形を底面、
        # 原点を頂点とする三角錐を考える。
        vector_60 = None  # その立体でunit_60を配置する方向
        vector_12 = MathVector(1. / 20. * math.sqrt(250. + 110. * math.sqrt(5.)),
                                1. / 10. * math.sqrt(25. + 10. * math.sqrt(5.)),
                                0)
        vector_20 = MathVector(1. / 12. * (3. * math.sqrt(3.) + math.sqrt(15.)),
                                0,
                                1. / 6. * math.sqrt(3.))
        vector_30 = MathVector(1, 0, 0)
        vector_12_20 = vector_12.cross(vector_20)
        vector_20_30 = vector_20.cross(vector_30)
        vector_30_12 = vector_30.cross(vector_12)

        gcone_orig = GeneralizedCone()
        gcone_orig.intersect_face(GeneralizedCone.Face(vector_12_20, math.pi))
        gcone_orig.intersect_face(GeneralizedCone.Face(vector_20_30, math.pi))
        gcone_orig.intersect_face(GeneralizedCone.Face(vector_30_12, math.pi))

        angle = math.pi / 2  # 画角
        dangle = angle  # 画角の誤差

        class _Found(Exception):
            pass

        for order in range(16):
            dangle /= 2
            gcone = GeneralizedCone(gcone_orig)
            if has_unit_12:
                gcone.intersect_face(GeneralizedCone.Face(vector_12.reflect(), math.pi * 2 - angle))
            if has_unit_20:
                gcone.intersect_face(GeneralizedCone.Face(vector_20.reflect(), math.pi * 2 - angle))
            if has_unit_30:
                gcone.intersect_face(GeneralizedCone.Face(vector_30.reflect(), math.pi * 2 - angle))
            try:
                if gcone.count_edge() == 0:
                    raise _Found
                if has_unit_60:
                    for i in range(gcone_orig.count_face()):
                        vector_orig = gcone_orig.get_face(i).vector
                        for j in range(gcone.count_edge()):
                            vector_a = gcone.get_edge(j).vector
                            for k in range(j + 1, gcone.count_edge()):
                                vector_b = gcone.get_edge(k).vector
                                vector_b = vector_orig.cross(vector_a.minus(vector_b))
                                if vector_b.x < 0:
                                    vector_b = vector_b.reflect()
                                for l in range(gcone.count_edge()):
                                    if vector_b.ang(gcone.get_edge(l).vector) > angle / 2:
                                        break
                                else:
                                    vector_60 = MathVector(vector_b)
                                    raise _Found
                            else:
                                vector_a_proj = vector_a.project(vector_orig)
                                for k in range(gcone.count_edge()):
                                    if vector_a_proj.ang(gcone.get_edge(k).vector) > angle / 2:
                                        break
                                else:
                                    vector_60 = MathVector(vector_a_proj)
                                    raise _Found
            except _Found:
                angle -= dangle
                continue
            angle += dangle

        self._angle = angle

        # 正十二面体でユニットを配置する方向を計算します。
        n_units = ((12 if has_unit_12 else 0) + (20 if has_unit_20 else 0)
                   + (30 if has_unit_30 else 0) + (60 if has_unit_60 else 0))
        vectors = [None] * n_units
        if has_unit_12:
            for i in range(12):
                vectors[i] = Polyhedron.DODECAHEDRON.get_face(i).vector
        if has_unit_20:
            offset = 12 if has_unit_12 else 0
            for i in range(20):
                vectors[offset + i] = Polyhedron.DODECAHEDRON.get_vertex(i).vector
        if has_unit_30:
            offset = (12 if has_unit_12 else 0) + (20 if has_unit_20 else 0)
            for i in range(30):
                vectors[offset + i] = Polyhedron.DODECAHEDRON.get_edge(i).vector
        if has_unit_60:
            k = vector_60.dot(vector_20_30) / vector_12.unit_vector().dot(vector_20_30)
            l = vector_60.dot(vector_30_12) / vector_20.unit_vector().dot(vector_30_12)
            m = vector_60.dot(vector_12_20) / vector_30.unit_vector().dot(vector_12_20)
            vectors_offset = ((12 if has_unit_12 else 0) + (20 if has_unit_20 else 0)
                              + (30 if has_unit_30 else 0))
            if k < l and k < m:
                for i in range(20):
                    vertex = Polyhedron.DODECAHEDRON.get_vertex(i)
                    for j in range(3):
                        vectors[vectors_offset + i * 3 + j] = (
                            vertex.vector.mult_scalar(l)
                            .plus(vertex.get_edge(j).vector.mult_scalar(m))
                            .unit_vector())
            elif l < m:
                for i in range(30):
                    edge = Polyhedron.DODECAHEDRON.get_edge(i)
                    for j in range(2):
                        vectors[vectors_offset + i * 2 + j] = (
                            edge.vector.mult_scalar(m)
                            .plus(edge.get_face(j).vector.mult_scalar(k))
                            .unit_vector())
            else:
                for i in range(20):
                    vertex = Polyhedron.DODECAHEDRON.get_vertex(i)
                    for j in range(3):
                        vectors[vectors_offset + i * 3 + j] = (
                            vertex.vector.mult_scalar(l)
                            .plus(vertex.get_face(j).vector.mult_scalar(k))
                            .unit_vector())

        # 正十二面体でユニットを配置する方向を経度、緯度の偶奇の順でソートします。
        n = len(vectors)
        for i in range(n):
            for j in range(i + 1, n):
                if vectors[i].get_lat() - _EPS > vectors[j].get_lat():
                    continue
                if vectors[i].get_lat() + _EPS > vectors[j].get_lat():
                    odd_i = int(math.floor((vectors[i].get_lng() + _EPS + math.pi * 2) * 5 / math.pi)) % 2
                    odd_j = int(math.floor((vectors[j].get_lng() + _EPS + math.pi * 2) * 5 / math.pi)) % 2
                    if odd_i < odd_j:
                        continue
                vectors[i], vectors[j] = vectors[j], vectors[i]

        self._vectors = vectors

    def get_unit_vector(self, i, j):
        return self._vectors[j if i == 0 else len(self._vectors) - 1 - j]

    def number_of_units(self, i):
        return ((6 if self._has_unit_12 else 0)
                + (10 if self._has_unit_20 else 0)
                + (15 if self._has_unit_30 else 0)
                + (30 if self._has_unit_60 else 0)
                + (int(math.ceil(math.pi / self._angle / 5)) * 5 if self._has_unit_extra else 0))


class CustomUnitArrangement(UnitArrangement):
    """手入力したユニットの配置を表現するクラス"""

    def __init__(self, north, south):
        self._vectors = [list(north), list(south)]

    def get_unit_vector(self, i, j):
        """指定された半球上にあるユニットの向きを表すベクトルを返します。"""
        return self._vectors[i][j]

    def number_of_units(self, i):
        """指定された半球上にあるユニットの個数を返します。"""
        return len(self._vectors[i])

    @staticmethod
    def _parse_vector(props, key):
        """プロパティファイルからベクトルを読み取ります。"""
        if props.get(key + ".xyz") is not None:
            coords = props[key + ".xyz"].split()
            if len(coords) != 3:
                raise InvalidUnitArrangementException(f"invalid vector format at key {key}")
            x, y, z = (float(c) for c in coords)
            return MathVector(x, y, z)
        if props.get(key + ".x") is not None:
            x = float(props[key + ".x"])
            y = float(props[key + ".y"])
            z = float(props[key + ".z"])
            return MathVector(x, y, z)
        if props.get(key + ".alpha_delta") is not None:
            coords = props[key + ".alpha_delta"].split()
            if len(coords) != 2:
                raise InvalidUnitArrangementException(f"invalid vector format at key {key}")
            alpha = float(coords[0])  # in hours
            delta = float(coords[1])  # in degrees
            return MathVector.from_mag_lng_lat(1., math.radians(alpha / 24. * 360.), math.radians(delta))
        if props.get(key + ".alpha") is not None:
            alpha = float(props[key + ".alpha"])  # in hours
            delta = float(props[key + ".delta"])  # in degrees
            return MathVector.from_mag_lng_lat(1., math.radians(alpha / 24. * 360.), math.radians(delta))
        return None

    @classmethod
    def from_properties(cls, props, key_prefix):
        north = []
        south = []
        i = 0
        while True:
            key = f"{key_prefix}[{i}]"
            v = cls._parse_vector(props, key)
            if v is None:
                break
            hemisphere = props.get(key + ".hemisphere")
            if hemisphere is not None:
                h = hemisphere.lower()
                if h in ("n", "north"):
                    north.append(v)
                    i += 1
                    continue
                if h in ("s", "south"):
                    south.append(v)
                    i += 1
                    continue
            if v.z >= 0:
                north.append(v)
            elif v.z < 0:
                south.append(v)
            else:  # NaN
                raise InvalidUnitArrangementException(f"invalid unit arrangement at {i}")
            i += 1
        return cls(north, south)

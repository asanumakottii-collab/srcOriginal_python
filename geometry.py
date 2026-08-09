# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2013 東京大学地文研究会天文部
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
"""錐体(GeneralizedCone)と多面体(Polyhedron)を扱うモジュールです。"""

import math

from mathvector import MathVector

# 浮動小数点の計算で生じる誤差よりも十分大きな数 ここでは 2^-16
_EPS = 0.0000152587890625


class GeneralizedCone:
    """錐体を扱うクラスです。(Perl5での成果を移植)"""

    class Face:
        """錐体の面を扱うクラスです。"""

        def __init__(self, vector, angle):
            self._vector = vector.unit_vector()
            self._angle = angle

        @property
        def vector(self):
            """中心の方向を取得します。"""
            return MathVector(self._vector)

        @property
        def angle(self):
            """直径の角度を取得します。"""
            return self._angle

        def contains(self, vector):
            """指定された方向を含んでいるかどうかを判定します。"""
            return self._angle / 2 - self._vector.ang(vector) >= 0

    class Edge:
        """錐体の辺を扱うクラスです。"""

        def __init__(self, vector):
            self._vector = vector.unit_vector()

        @property
        def vector(self):
            """辺の方向を取得します。"""
            return MathVector(self._vector)

    def __init__(self, gcone=None):
        self._faces = []
        self._edges = []
        if gcone is not None:
            self._faces.extend(gcone._faces)
            self._edges.extend(gcone._edges)

    def count_face(self):
        """面の個数を取得します。"""
        return len(self._faces)

    def get_face(self, index):
        """面を取得します。"""
        return self._faces[index]

    def count_edge(self):
        """辺の個数を取得します。"""
        return len(self._edges)

    def get_edge(self, index):
        """辺を取得します。"""
        return self._edges[index]

    def intersect_face(self, face):
        """指定された面を重ね合わせます。"""
        self._faces.append(face)
        self._edges = [e for e in self._edges if face.contains(e.vector)]

        class _Found(Exception):
            pass

        for i in range(self.count_face() - 1):
            simultaneous_face = self.get_face(i)
            # 線形方程式系
            system = [[0.0] * 4, [0.0] * 4]
            # 対称群
            sym = [0, 1, 2]
            v = simultaneous_face.vector
            system[0][0], system[0][1], system[0][2] = v.x, v.y, v.z
            system[0][3] = math.cos(simultaneous_face.angle / 2)
            v = face.vector
            system[1][0], system[1][1], system[1][2] = v.x, v.y, v.z
            system[1][3] = math.cos(face.angle / 2)

            try:
                for j in range(2):
                    # 列を入れ替える
                    k = j
                    while abs(system[j][j]) < _EPS:
                        k += 1
                        if k == 3:
                            raise _Found  # 解が求まらないので次のiへ(continue equation)
                        for l in range(2):
                            system[l][j], system[l][k] = system[l][k], system[l][j]
                        sym[j], sym[k] = sym[k], sym[j]
                    # 行を実数倍してほかの行に加える
                    for k in range(3, j - 1, -1):
                        system[j][k] /= system[j][j]
                        for l in range(2):
                            if j != l:
                                system[l][k] -= system[j][k] * system[l][j]

                # 直線の式を球の方程式に代入して交点の極座標を求めます。
                a, b, c = 1.0, 0.0, -1.0  # 2次の係数, 1次の係数の1/2, 定数
                for j in range(2):
                    a += system[j][2] * system[j][2]
                    b -= system[j][2] * system[j][3]
                    c += system[j][3] * system[j][3]
                d = b * b - a * c  # 判別式
                if d < 0:
                    raise _Found  # 解が求まらないので次のiへ(continue equation)

                for j in range(2):
                    components = [0.0, 0.0, 0.0]
                    components[sym[2]] = (-b + math.sqrt(d) * (-1 if j != 0 else 1)) / a
                    components[sym[1]] = -components[sym[2]] * system[1][2] + system[1][3]
                    components[sym[0]] = -components[sym[2]] * system[0][2] + system[0][3]
                    vector = MathVector(components[0], components[1], components[2])
                    for k in range(self.count_face() - 1):
                        if i != k and not self.get_face(k).contains(vector):
                            break
                    else:
                        self._edges.append(GeneralizedCone.Edge(vector))
            except _Found:
                continue


class Polyhedron:
    """多面体を扱うクラスです。"""

    class Face:
        """多面体の面を扱うクラスです。"""

        def __init__(self, vector):
            self._vector = vector.unit_vector()

        @property
        def vector(self):
            """原点から面への法線の方向を取得します。"""
            return MathVector(self._vector)

    class Edge:
        """多面体の辺を扱うクラスです。"""

        def __init__(self, vector, faces):
            self._vector = vector.unit_vector()
            self._faces = list(faces)

        @property
        def vector(self):
            """原点から辺の中点への方向を取得します。"""
            return MathVector(self._vector)

        def count_face(self):
            return len(self._faces)

        def get_face(self, index):
            return self._faces[index]

    class Vertex:
        """多面体の頂点を扱うクラスです。"""

        def __init__(self, vector, faces, edges):
            self._vector = vector.unit_vector()
            self._faces = list(faces)
            self._edges = list(edges)

        @property
        def vector(self):
            """原点から頂点への方向を取得します。"""
            return MathVector(self._vector)

        def count_face(self):
            return len(self._faces)

        def get_face(self, index):
            return self._faces[index]

        def count_edge(self):
            return len(self._edges)

        def get_edge(self, index):
            return self._edges[index]

    def __init__(self, faces, edges, vertexes):
        self._faces = faces
        self._edges = edges
        self._vertexes = vertexes

    def count_face(self):
        return len(self._faces)

    def get_face(self, index):
        return self._faces[index]

    def count_edge(self):
        return len(self._edges)

    def get_edge(self, index):
        return self._edges[index]

    def count_vertex(self):
        return len(self._vertexes)

    def get_vertex(self, index):
        return self._vertexes[index]

    @classmethod
    def _build_dodecahedron(cls):
        faces = [None] * 12
        for i in range(2):
            vector = MathVector.from_mag_lng_lat(1., 0., (1. if i % 2 == 0 else -1.) * math.pi / 2)
            faces[i] = cls.Face(vector)
        for i in range(10):
            vector = MathVector.from_mag_lng_lat(
                1., math.pi * i / 5.,
                (1. if i % 2 == 0 else -1.) * (math.pi / 2 - math.acos(1. / math.sqrt(5.))))
            faces[i + 2] = cls.Face(vector)

        edges = [None] * 30
        for i in range(10):
            f = [faces[i % 2], faces[i + 2]]
            edges[i] = cls.Edge(f[0].vector.plus(f[1].vector), f)
        for i in range(10):
            f = [faces[i + 2], faces[(i + 2) % 10 + 2]]
            edges[i + 10] = cls.Edge(f[0].vector.plus(f[1].vector), f)
        for i in range(10):
            f = [faces[i + 2], faces[(i + 1) % 10 + 2]]
            edges[i + 20] = cls.Edge(f[0].vector.plus(f[1].vector), f)

        vertexes = [None] * 20
        for i in range(10):
            f = [faces[i % 2], faces[i + 2], faces[(i + 2) % 10 + 2]]
            e = [edges[i], edges[(i + 2) % 10], edges[i + 10]]
            vertexes[i] = cls.Vertex(f[0].vector.plus(f[1].vector).plus(f[2].vector), f, e)
        for i in range(10):
            f = [faces[i + 2], faces[(i + 2) % 10 + 2], faces[(i + 1) % 10 + 2]]
            e = [edges[i + 10], edges[i + 20], edges[(i + 1) % 10 + 20]]
            vertexes[i + 10] = cls.Vertex(f[0].vector.plus(f[1].vector).plus(f[2].vector), f, e)

        return cls(faces, edges, vertexes)


# 正十二面体
Polyhedron.DODECAHEDRON = Polyhedron._build_dodecahedron()

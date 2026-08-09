# -*- coding: utf-8 -*-
"""各原盤が担当する天球領域を原盤座標上のポリゴンとして生成します。"""

import math


class PlateAssignmentPolygonGenerator:
    """Transformer のユニット選択境界を極座標サンプリングで求めます。

    各原盤の中心から一定角度ごとに外向きへ進み、その位置の光線を天球へ
    逆投影します。担当ユニットが切り替わる半径を二分探索することで、旧
    DrawPoly.ipynb の星分布に依存した手作業を自動化します。
    """

    def __init__(self, transformer, samples=180, iterations=22):
        if samples < 12:
            raise ValueError("polygon.samples には12以上を指定してください。")
        if iterations < 1:
            raise ValueError("iterations には1以上を指定してください。")
        self.transformer = transformer
        self.samples = samples
        self.iterations = iterations

    def _is_assigned(self, dir_, index, xmm, ymm):
        position = self.transformer.inverse_transform_unit(dir_, index, xmm, ymm)
        return (
            position is not None
            and self.transformer.assigned_unit(position) == (dir_, index)
        )

    def generate(self, dir_, index, frame_radius):
        """1枚の原盤について、中心相対mm座標の頂点列を返します。"""
        if frame_radius <= 0:
            raise ValueError("frame_radius には正の値を指定してください。")
        if not self._is_assigned(dir_, index, 0., 0.):
            raise ValueError(f"原盤中心が担当領域に含まれません: {dir_}, {index}")

        points = []
        for sample in range(self.samples):
            theta = math.tau * sample / self.samples
            cos_theta = math.cos(theta)
            sin_theta = math.sin(theta)
            if self._is_assigned(
                dir_, index, frame_radius * cos_theta, frame_radius * sin_theta
            ):
                boundary = frame_radius
            else:
                inside = 0.
                outside = frame_radius
                for _ in range(self.iterations):
                    middle = (inside + outside) / 2
                    if self._is_assigned(
                        dir_, index, middle * cos_theta, middle * sin_theta
                    ):
                        inside = middle
                    else:
                        outside = middle
                boundary = inside
            points.append((boundary * cos_theta, boundary * sin_theta))
        return points

    def generate_all(self, frame_radius):
        """全原盤の ``(半球, 番号, 頂点列)`` を順に返します。"""
        for dir_ in range(2):
            for index in range(self.transformer.number_of_units(dir_)):
                yield dir_, index, self.generate(dir_, index, frame_radius)

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
"""原板データをSVGまたはPostScriptとして書き出すモジュールです。"""

import math
import os
from abc import ABC, abstractmethod
from enum import Enum

DEFAULT_OUTPUT_DIR = "output"

# SVGのpoints属性にはmmなどの単位を付けられないため、絶対長で記述した円・文字と座標を揃える際にCSS標準の96 dpiへ変換する。
_CSS_PIXELS_PER_MM = 96. / 25.4


def resolve_output_path(output_dir, filename):
    """出力フォルダを作成し、ファイルの出力パスを返します。"""
    output_dir = (output_dir or DEFAULT_OUTPUT_DIR).strip() or DEFAULT_OUTPUT_DIR
    if output_dir == "." or os.path.isabs(filename):
        path = filename
    else:
        path = os.path.join(output_dir, filename)
    dirname = os.path.dirname(path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    return path


class PlateWriterType(Enum):
    SVG = "SVG"
    POSTSCRIPT = "PostScript"


class PlateWriter(ABC):
    @abstractmethod
    def write_star(self, s):
        raise NotImplementedError

    @abstractmethod
    def write_constellation(self, c):
        raise NotImplementedError

    @abstractmethod
    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class _PlateWriterBase(PlateWriter):
    """PlateWriterPS と PlateWriterSVG に共通するレイアウト計算。"""

    WIDTH = 191.  # 出力用紙の幅
    HEIGHT = 277. / 2.  # 出力用紙の高さ

    def __init__(self, column, row, r, shape, filename_prefix, invert_color, output_dir=DEFAULT_OUTPUT_DIR):
        if column <= 0:
            raise ValueError("横の値が不正です。")
        if row <= 0:
            raise ValueError("縦の値が不正です。")
        self.outs = []
        self.number_of_plate = 0  # 出力用紙1ページ中の原板の数
        self.number_of_page = 0  # 出力用紙のページ数
        self.column = column
        self.row = row
        self.r = (min(self.WIDTH / column, self.HEIGHT / row) - 10.) / 2 if r == 0 else r
        self.shape = shape  # 原板が円形ならTrue, 原板が長方形ならFalse
        self.filename_prefix = filename_prefix
        self.invert_color = invert_color  # Falseなら原板を白で星を黒, Trueなら原板を黒で星を白
        self.output_dir = output_dir

    def _get_cx(self, index):
        """枠の中心のx座標"""
        return self.WIDTH / (self.column * 2) * ((index % (self.column * self.row)) // self.row * 2 + 1)

    def _get_cy(self, dir_, index):
        """枠の中心のy座標"""
        return self.HEIGHT / (self.row * 2) * ((index % self.row) * 2 + 1) + self.HEIGHT * dir_

    def _output_path(self, extension):
        return resolve_output_path(self.output_dir, f"{self.filename_prefix}{self.number_of_page}.{extension}")

    def __del__(self):
        if getattr(self, "outs", None) is not None:
            try:
                self.close()
            except Exception:
                pass


class PlateWriterPS(_PlateWriterBase):
    def _write_frame(self, index):
        page = index // (self.column * self.row)
        while page >= self.number_of_page:
            out = open(self._output_path("ps"), "w", encoding="ascii")
            out.write("%!PS-Adobe-3.0\n")
            out.write("%%Pages: 1\n")
            out.write("%%Page: 1page 1\n")
            out.write("/circ { 0 360 arc } def\n")
            out.write("/strmiddle { dup stringwidth pop 2 div neg 0 rmoveto } def\n")
            out.write("/showmiddle { strmiddle show } def\n")
            out.write("2.834646 -2.834646 scale\n")  # mm単位、y軸は下方向
            out.write(f"0 {-self.HEIGHT * 2} translate\n")
            out.write("/Verdana findfont 3 scalefont setfont\n")
            self.outs.append(out)
            self.number_of_page += 1
        while index >= self.number_of_plate:
            out = self.outs[self.number_of_plate // (self.column * self.row)]
            cx = self._get_cx(self.number_of_plate)
            for i in range(2):
                cy = self._get_cy(i, self.number_of_plate)
                label = ("N" if i == 0 else "S") + str(self.number_of_plate)
                out.write(f"gsave {cx} {cy - self.r} translate 1 -1 scale 0 0 moveto ({label}) showmiddle grestore\n")
                out.write("0.1 setlinewidth\n")
                if self.shape:
                    out.write(f"{cx} {cy} {self.r} circ "
                              f"{'0 setgray fill 1 setgray' if self.invert_color else 'stroke'}\n")
                    out.write(f"newpath {cx - self.r * 3. / 4.} {cy - self.r * 3. / 4.} moveto "
                              f"{cx - self.r} {cy - self.r} lineto stroke\n")
                    out.write(f"newpath {cx - self.r * 3. / 4.} {cy - self.r} moveto "
                              f"{cx - self.r} {cy - self.r * 3. / 4.} lineto stroke\n")
                else:
                    x, y = cx - self.r, cy - self.r
                    w, h = self.r * 2, self.r * 2
                    out.write(f"newpath {x} {y} moveto {w} 0 rlineto 0 {h} rlineto {-w} 0 rlineto closepath "
                              f"{'0 setgray fill' if self.invert_color else 'stroke'}\n")
            self.number_of_plate += 1

    def write_star(self, s):
        if self.outs is None:
            raise IOError("writer is closed")
        self._write_frame(s.p.index)
        if (abs(s.p.xmm) > self.r or abs(s.p.ymm) > self.r
                or math.isnan(s.p.xmm) or math.isnan(s.p.ymm)):
            return
        out = self.outs[s.p.index // (self.column * self.row)]
        cx = self._get_cx(s.p.index)
        cy = self._get_cy(s.p.dir, s.p.index)
        if self.invert_color:
            out.write("1 setgray ")  # white
        out.write(f"{cx + s.p.xmm} {cy + s.p.ymm} {s.rmm} circ fill\n")

    def write_constellation(self, c):
        if self.outs is None:
            raise IOError("writer is closed")
        if c.name is not None and c.p is not None:
            self._write_frame(c.p.index)
            out = self.outs[c.p.index // (self.column * self.row)]
            cx = self._get_cx(c.p.index)
            cy = self._get_cy(c.p.dir, c.p.index)
            if abs(c.p.xmm) <= self.r and abs(c.p.ymm) <= self.r:
                out.write(f"gsave {cx + c.p.xmm} {cy + c.p.ymm} translate 1 -1 scale 0 0 moveto "
                          f"({c.name}) showmiddle grestore\n")
        if c.ll is not None:
            for line in c.ll:
                self._write_frame(line[0].index)
                out = self.outs[line[0].index // (self.column * self.row)]
                cx = self._get_cx(line[0].index)
                cy = self._get_cy(line[0].dir, line[0].index)
                if ((abs(line[0].xmm) <= self.r and abs(line[0].ymm) <= self.r)
                        or (abs(line[1].xmm) <= self.r and abs(line[1].ymm) <= self.r)):
                    if not self.shape:
                        x, y = cx - self.r, cy - self.r
                        w, h = self.r * 2, self.r * 2
                        out.write(f"initclip newpath {x} {y} moveto {w} 0 rlineto 0 {h} rlineto "
                                  f"{-w} 0 rlineto closepath clip\n")
                    if self.invert_color:
                        out.write("1 setgray\n")  # white
                    out.write(f"0.1 setlinewidth newpath {cx + line[0].xmm} {cy + line[0].ymm} moveto "
                              f"{cx + line[1].xmm} {cy + line[1].ymm} lineto stroke\n")
                    if not self.shape:
                        out.write("initclip\n")

    def close(self):
        if self.outs is None:
            raise IOError("writer is closed")
        for out in self.outs:
            out.write("showpage\n")
            out.write("%%EOF\n")
            out.close()
        self.outs = None


class PlateWriterSVG(_PlateWriterBase):
    def _write_frame(self, index):
        page = index // (self.column * self.row)
        while page >= self.number_of_page:
            out = open(self._output_path("svg"), "w", encoding="utf-8")
            out.write("<?xml version=\"1.0\"?>")
            out.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
                      "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">")
            out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{self.WIDTH}mm\" "
                      f"height=\"{self.HEIGHT * 2}mm\" version=\"1.1\">")
            self.outs.append(out)
            self.number_of_page += 1
        while index >= self.number_of_plate:
            out = self.outs[self.number_of_plate // (self.column * self.row)]
            cx = self._get_cx(self.number_of_plate)
            for i in range(2):
                cy = self._get_cy(i, self.number_of_plate)
                label = ("N" if i == 0 else "S") + str(self.number_of_plate)
                out.write(f"<text x=\"{cx}mm\" y=\"{cy - self.r}mm\" text-anchor=\"middle\" "
                          f"font-size=\"3mm\" font-familiy=\"Verdana\">{label}</text>")
                if self.shape:
                    out.write(f"<circle cx=\"{cx}mm\" cy=\"{cy}mm\" r=\"{self.r}mm\" "
                              f"fill=\"{'black' if self.invert_color else 'none'}\" "
                              f"stroke=\"{'none' if self.invert_color else 'black'}\" "
                              f"stroke-width=\"0.1mm\" />")
                    out.write(f"<line x1=\"{cx - self.r * 3. / 4.}mm\" y1=\"{cy - self.r * 3. / 4.}mm\" "
                              f"x2=\"{cx - self.r}mm\" y2=\"{cy - self.r}mm\" "
                              f"stroke=\"black\" stroke-width=\"0.1mm\" />")
                    out.write(f"<line x1=\"{cx - self.r * 3. / 4.}mm\" y1=\"{cy - self.r}mm\" "
                              f"x2=\"{cx - self.r}mm\" y2=\"{cy - self.r * 3. / 4.}mm\" "
                              f"stroke=\"black\" stroke-width=\"0.1mm\" />")
                else:
                    out.write(f"<rect x=\"{cx - self.r}mm\" y=\"{cy - self.r}mm\" "
                              f"width=\"{self.r * 2}mm\" height=\"{self.r * 2}mm\" "
                              f"fill=\"{'black' if self.invert_color else 'none'}\" "
                              f"stroke=\"{'white' if self.invert_color else 'black'}\" "
                              f"stroke-width=\"0.1mm\" />")
            self.number_of_plate += 1

    def write_star(self, s):
        if self.outs is None:
            raise IOError("writer is closed")
        self._write_frame(s.p.index)
        if abs(s.p.xmm) > self.r or abs(s.p.ymm) > self.r:
            return
        out = self.outs[s.p.index // (self.column * self.row)]
        cx = self._get_cx(s.p.index)
        cy = self._get_cy(s.p.dir, s.p.index)
        out.write(f"<circle cx=\"{cx + s.p.xmm}mm\" cy=\"{cy + s.p.ymm}mm\" r=\"{s.rmm}mm\" "
                  f"fill=\"{'white' if self.invert_color else 'black'}\"/>")

    def write_constellation(self, c):
        if self.outs is None:
            raise IOError("writer is closed")
        if c.name is not None and c.p is not None:
            self._write_frame(c.p.index)
            out = self.outs[c.p.index // (self.column * self.row)]
            cx = self._get_cx(c.p.index)
            cy = self._get_cy(c.p.dir, c.p.index)
            fill_attr = " fill=\"white\"" if self.invert_color else ""
            out.write(f"<text x=\"{cx + c.p.xmm}mm\" y=\"{cy + c.p.ymm}mm\" text-anchor=\"middle\" "
                      f"stroke-width=\"3mm\" font-size=\"3mm\" font-familiy=\"Verdana\"{fill_attr}>"
                      f"{c.name}</text>")
        if c.ll is not None:
            for line in c.ll:
                self._write_frame(line[0].index)
                out = self.outs[line[0].index // (self.column * self.row)]
                cx = self._get_cx(line[0].index)
                cy = self._get_cy(line[0].dir, line[0].index)
                out.write(f"<line x1=\"{cx + line[0].xmm}mm\" y1=\"{cy + line[0].ymm}mm\" "
                          f"x2=\"{cx + line[1].xmm}mm\" y2=\"{cy + line[1].ymm}mm\" "
                          f"stroke=\"{'white' if self.invert_color else 'black'}\" stroke-width=\"0.1mm\" />")

    def write_assignment_polygon(self, dir_, index, points):
        """指定原盤の担当星域を塗りつぶすポリゴンを書き出します。"""
        if self.outs is None:
            raise IOError("writer is closed")
        if not points:
            return
        self._write_frame(index)
        out = self.outs[index // (self.column * self.row)]
        cx = self._get_cx(index)
        cy = self._get_cy(dir_, index)
        point_text = " ".join(
            f"{(cx + x) * _CSS_PIXELS_PER_MM},{(cy + y) * _CSS_PIXELS_PER_MM}"
            for x, y in points
        )
        label = ("N" if dir_ == 0 else "S") + str(index)
        out.write(
            f"<polygon data-plate=\"{label}\" points=\"{point_text}\" "
            f"fill=\"white\" stroke=\"white\" stroke-width=\"0.1mm\" />"
        )

    def close(self):
        if self.outs is None:
            raise IOError("writer is closed")
        for out in self.outs:
            out.write("</svg>")
            out.close()
        self.outs = None

"""RC3銀河を I(R) = I0 exp(-R/h) の楕円状点群にするための尺度。

長径/短径比を q とし、R = sqrt(x**2 + (q*y)**2) を使う。
距離は秒角単位。点の半径は Gamma(shape=2, scale=h) から生成する。
"""

import math

import numpy as np


# 四分円でのGauss-Legendre積分。平均値なので重みの合計を1にする。
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)
_ANGLES = (_NODES + 1.0) * (math.pi / 4.0)
_WEIGHTS = _WEIGHTS / 2.0
_HALF_LIGHT_RADIUS_IN_SCALES = 1.6783469900166605
_FALLBACK_ISOPHOTE_IN_SCALES = 3.0


def _scale_from_circular_half_light(radius, axis_ratio):
    """円形開口内の光量を1/2にするh。楕円の半光量半長径とは区別する。"""
    if axis_ratio == 1.0:
        return radius / _HALF_LIGHT_RADIUS_IN_SCALES
    # 楕円座標角thetaで、円形開口までの楕円半径はradius/g(theta)。
    g = np.sqrt(np.cos(_ANGLES) ** 2 + (np.sin(_ANGLES) / axis_ratio) ** 2)
    lower = _HALF_LIGHT_RADIUS_IN_SCALES / axis_ratio
    upper = _HALF_LIGHT_RADIUS_IN_SCALES
    for _ in range(60):
        z = (lower + upper) / 2.0  # radius/h
        t = z / g
        # Gamma(2, 1)のCDF。小さいtでも桁落ちを抑える。
        enclosed = float(np.dot(_WEIGHTS, -np.expm1(-t) - t * np.exp(-t)))
        if enclosed < 0.5:
            lower = z
        else:
            upper = z
    return radius / ((lower + upper) / 2.0)


def exponential_scale_arcsec(a25, axis_ratio, bmag, half_light_radius=None):
    """尺度hを返す。大きさの情報がなければNone。

    a25: D25の半長径、half_light_radius: Aeの円形開口半径（いずれも秒角）。
    bmag: 未補正の総B等級。Vしか分からない場合はNone。
    Aeを優先し、なければ総B光量とB=25 mag/arcsec^2の輪郭を合わせる。
    二解のうちa25/h >= 2の集中した解を選ぶ。解なし/B不明ではh=a25/3。
    """
    if not math.isfinite(axis_ratio) or axis_ratio < 1.0:
        raise ValueError("銀河の長径/短径比は1以上の有限値である必要があります。")
    if (half_light_radius is not None and math.isfinite(half_light_radius)
            and half_light_radius > 0.0):
        return _scale_from_circular_half_light(half_light_radius, axis_ratio)
    if a25 is None or not math.isfinite(a25) or a25 <= 0.0:
        return None
    fallback = a25 / _FALLBACK_ISOPHOTE_IN_SCALES
    if bmag is None or not math.isfinite(bmag):
        return fallback

    # t=a25/h: t^2 exp(-t) = (2*pi*a25^2/q) * I25,B/F_B。
    # 対数で解き、暗い/大きな銀河でもオーバーフローを避ける。
    log_c = (math.log(2.0 * math.pi / axis_ratio) + 2.0 * math.log(a25)
             + 0.4 * math.log(10.0) * (bmag - 25.0))
    maximum = 2.0 * math.log(2.0) - 2.0
    if math.isclose(log_c, maximum, rel_tol=0.0, abs_tol=1e-12):
        return a25 / 2.0
    if log_c > maximum:
        return fallback
    lower = 2.0
    upper = max(4.0, -log_c + 2.0)
    while 2.0 * math.log(upper) - upper > log_c:
        upper *= 2.0
    for _ in range(60):
        t = (lower + upper) / 2.0
        if 2.0 * math.log(t) - t > log_c:
            lower = t
        else:
            upper = t
    return a25 / ((lower + upper) / 2.0)

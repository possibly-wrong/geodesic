# Reference: T. Vincenty (1975) Direct and inverse solutions of geodesics on
# the ellipsoid with application of nested equations, Survey Review, 23:176,
# 88-93, with Fortran implementation from U.S. National Geodetic Survey
# (https://geodesy.noaa.gov/PC_PROD/Inv_Fwd/)

# Results are undefined for nearly antipodal points, i.e., WGS-84 geodesics
# longer than about 19000 km.

import math

def direct(lat1, lon1, faz, s, a=6378137, f=1/298.257223563):
    eps = 0.5e-13
    r = 1 - f
    tu = r * math.sin(lat1) / math.cos(lat1)
    sf = math.sin(faz)
    cf = math.cos(faz)
    baz = (0 if cf == 0 else math.atan2(tu, cf) * 2)
    cu = 1 / math.sqrt(tu * tu + 1)
    su = tu * cu
    sa = cu * sf
    c2a = -sa * sa + 1
    x = math.sqrt((1 / r / r - 1) * c2a + 1) + 1
    x = (x - 2) / x
    c = 1 - x
    c = (x * x / 4 + 1) / c
    d = (0.375 * x * x - 1) * x
    tu = s / r / a / c
    y = tu
    while True:
        sy = math.sin(y)
        cy = math.cos(y)
        cz = math.cos(baz + y)
        e = cz * cz * 2 - 1
        c = y
        x = e * cy
        y = e + e - 1
        y = ((((sy * sy * 4 - 3) * y * cz * d / 6 + x) * d / 4 - cz) *
            sy * d + tu)
        if abs(y - c) <= eps:
            break
    baz = cu * cy * cf - su * sy
    c = r * math.sqrt(sa * sa + baz * baz)
    d = su * cy + cu * sy * cf
    lat2 = math.atan2(d, c)
    c = cu * cy - su * sy * cf
    x = math.atan2(sy * sf, c)
    c = ((-3 * c2a + 4) * f + 4) * c2a * f / 16
    d = ((e * cy * c + cz) * sy * c + y) * sa
    lon2 = lon1 + x - (1 - c) * d * f
    return (lat2, lon2)

def inverse(lat1, lon1, lat2, lon2, a=6378137, f=1/298.257223563):
    eps = 1e-12
    r = 1 - f
    U1 = math.atan(r * math.tan(lat1))
    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    U2 = math.atan(r * math.tan(lat2))
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    L = (lon2 - lon1) % (2 * math.pi)
    if L > math.pi:
        L -= 2 * math.pi
    if L < -math.pi:
        L += 2 * math.pi
    lam = L
    for i in range(50):
        sinlambda = math.sin(lam)
        coslambda = math.cos(lam)
        y = cosU2 * sinlambda
        x = cosU1 * sinU2 - sinU1 * cosU2 * coslambda
        sinsigma = math.sqrt(y * y + x * x)
        cossigma = sinU1 * sinU2 + cosU1 * cosU2 * coslambda
        sigma = math.atan2(sinsigma, cossigma)
        sinalpha = cosU1 * cosU2 * sinlambda / max(sinsigma, 1e-15)
        cosalpha2 = 1 - sinalpha * sinalpha
        cos2sigmam = cossigma - 2 * sinU1 * sinU2 / max(cosalpha2, 1e-15)
        C = f / 16 * cosalpha2 * (4 + f * (4 - 3 * cosalpha2))
        prev = lam
        lam = L + (1 - C) * f * sinalpha * (sigma + C * sinsigma * (
            cos2sigmam + C * cossigma * (
                -1 + 2 * cos2sigmam * cos2sigmam)))
        if abs(lam - prev) <= eps:
            break
    sqrt1u2 = math.sqrt(1 + cosalpha2 * (1 / (r * r) - 1))
    k1 = (sqrt1u2 - 1) / (sqrt1u2 + 1)
    A = (1 + k1 * k1 / 4) / (1 - k1)
    B = k1 * (1 - 0.375 * k1 * k1)
    dsigma = B * sinsigma * (cos2sigmam + B / 4 * (cossigma * (
        -1 + 2 * cos2sigmam * cos2sigmam) - B / 6 * cos2sigmam * (
        -3 + 4 * sinsigma * sinsigma) * (
        -3 + 4 * cos2sigmam * cos2sigmam)))
    s = r * a * A * (sigma - dsigma)
    faz = math.atan2(y, x)
    return (faz, s)

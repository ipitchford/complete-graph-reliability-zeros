#!/usr/bin/env python3
"""Certify bounds on Waldvogel's rotation-restricted constant W_0.

All certificate inequalities use outward-rounded integer intervals.
An optional independent backend uses mpmath.iv. Each leaf excludes F(theta,phi)=0 on one phase rectangle at a
fixed radius. At theta=0 the function exp(z) has no zeros. Homotopy in theta
then excludes roots throughout the disk. Complex conjugation supplies the
other half of the theta circle. A separate Rouche disk proves the upper bound.
"""
from functools import lru_cache
from math import factorial
import sys


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)

if '--mpmath' in sys.argv:
    from mpmath import iv
    iv.dps = 42
else:
    from certified_intervals import iv
DEGREE = 25
LOWER = '0.73775075747822433247'
UPPER = '0.73775075747822433248'
RADIUS = iv.mpf(LOWER)
RADIAL_BAND = iv.mpf([LOWER, UPPER])
POWERS = [[r**k / factorial(k) for k in range(DEGREE+1)]
          for r in (RADIUS, RADIAL_BAND)]
THETA = iv.mpf(['2.39224', '2.39425'])
PHI = iv.mpf(['3.02767', '3.02968'])
BINS = [k*(k-1)//2 for k in range(DEGREE+1)]
# All omitted terms on |z|<1 sum to < 2/(d+1)! regardless of phase.
TAIL = iv.mpf(2)/factorial(DEGREE+1)


@lru_cache(maxsize=32768)
def phase(numerator, denominator):
    """cis(pi*n/d); all callers reduce n modulo 2*d first."""
    angle = iv.pi * numerator / denominator
    return iv.cos(angle), iv.sin(angle)


@lru_cache(maxsize=64)
def box_error(d, band):
    # theta halfwidth pi/d, phi halfwidth 2*pi/d.
    h = iv.pi/d
    remainder = TAIL + sum(POWERS[band][k] * ((2*k+BINS[k])*h)**2
                           for k in range(1, DEGREE+1))/2
    return h, remainder


def cross(x, y):
    return x[0]*y[1] - x[1]*y[0]


def norm(z):
    return iv.sqrt(z[0]**2 + z[1]**2)


def excludes(level, i, j, band):
    """First-order Taylor test, optionally projected off the phi tangent."""
    d = 2**(level+1)
    u, v = 2*i+1, 2*(2*j+1)-d
    f = [iv.mpf(1), iv.mpf(0)]
    f_theta = [iv.mpf(0), iv.mpf(0)]
    f_phi = [iv.mpf(0), iv.mpf(0)]
    for k in range(1, DEGREE+1):
        b = BINS[k]
        cs, sn = phase((k*v-b*u)%(2*d), d)
        x, y = POWERS[band][k]*cs, POWERS[band][k]*sn
        f[0] += x
        f[1] += y
        f_theta[0] += b*y
        f_theta[1] -= b*x
        f_phi[0] -= k*y
        f_phi[1] += k*x
    h, remainder = box_error(d, band)
    for v in range(2):
        if abs(f[v]) > h*(abs(f_theta[v])+2*abs(f_phi[v])) + remainder:
            return True
    # The determinant with f_phi eliminates the first-order phi term.
    return abs(cross(f, f_phi)) > h*abs(cross(f_theta, f_phi)) + remainder*norm(f_phi)


def cover(band=False):
    stack = [(0, 0, 0)]
    resolved = 0
    deepest = 0
    while stack:
        level, i, j = stack.pop()
        d = 2**(level+1)
        if band:
            theta = iv.mpf([2*i, 2*i+2])*iv.pi/d
            phi = iv.mpf([4*j-d, 4*j+4-d])*iv.pi/d
            if theta in THETA and phi in PHI:
                resolved += 1
                deepest = max(deepest, level)
                continue
        if excludes(level, i, j, int(band)):
            resolved += 1
            deepest = max(deepest, level)
        else:
            if level == 48:
                raise ArithmeticError('insufficient depth for claimed lower endpoint')
            stack.extend((level+1, 2*i+a, 2*j+b) for a in (0,1) for b in (0,1))
    if band:
        print('All zeros in the enclosed radial band have phases in the stated patch, up to conjugation.')
    else:
        print('Certified no zeros on |z| =', LOWER, 'for all unit phases.')
    print('Terminal phase rectangles:', resolved, 'maximum depth:', deepest)


def derivatives(z, theta):
    # Each derivative of order at most two multiplies term k by at most k^4,
    # apart from reducing powers of z by at most 2. Here |z| <= 1.
    require(abs(z) < 1)
    values = [iv.mpc(0) for _ in range(6)]
    for k in range(DEGREE + 1):
        b = BINS[k]
        argument = -b*theta
        phase = iv.mpc(iv.cos(argument), iv.sin(argument))/factorial(k)
        term = phase*z**k
        first_z = k*phase*z**(k-1) if k else 0
        second_z = k*(k-1)*phase*z**(k-2) if k >= 2 else 0
        values[0] += term
        values[1] += first_z
        values[2] -= 1j*b*term
        values[3] += second_z
        values[4] -= 1j*b*first_z
        values[5] -= b**2*term
    # For k >= d+1 the bound k^4/k! decreases by a factor < 1/2.
    tail = 2*(DEGREE+1)**4 * iv.mpf([-1, 1])/factorial(DEGREE+1)
    return [value + iv.mpc(tail,tail) for value in values]


def minimiser_uniqueness():
    central = iv.mpc('-0.733052523296829803468871187856', '0.083127481904870509963603546822')
    eps = iv.mpf('0.002')
    # Every root admitted by the band covering lies in this same small disk.
    patch = RADIAL_BAND * iv.mpc(iv.cos(PHI), iv.sin(PHI))
    require(abs(patch-central) < eps)
    f, fz, *_ = derivatives(central, THETA)
    require(abs(f) + 3*eps**2/2 < abs(fz)*eps)
    # Rouche gives a unique simple smooth root throughout the theta interval.
    radius = iv.mpf(['-0.002', '0.002'])
    square = central + iv.mpc(radius, radius)
    _, fz, ft, fzz, ftz, ftt = derivatives(square, THETA)
    require(abs(fz)>1)
    zp = -ft/fz
    zpp = -(ftt + 2*ftz*zp + fzz*zp**2)/fz
    curvature = 2*(abs(zp)**2 + zpp.real*square.real + zpp.imag*square.imag)
    require(curvature > 0)
    print('All minimisers with 0 <= theta <= pi have 2.39224 <= theta <= 2.39425.')
    print('Rouche gives one simple smooth branch there. Lower bound on the second')
    print('derivative of its squared modulus:', curvature.a)
    # Check opposite slope signs on either side of the certified minimiser.
    endpoints = [
        ('2.39324227413489', '-0.733052523296829652446548891581', '0.083127481904871841740915893735', -1),
        ('2.39324227413490', '-0.733052523296829882222455924863', '0.083127481904869815481914853325', 1),
    ]
    for angle, re, im, direction in endpoints:
        theta, seed, eps = iv.mpf(angle), iv.mpc(re, im), iv.mpf('1e-19')
        require(abs(seed-central)+eps < iv.mpf('0.002'))
        value, d_z, *_ = derivatives(seed, theta)
        require(abs(value)+3*eps**2/2 < abs(d_z)*eps)
        side = iv.mpf(['-1e-19', '1e-19'])
        square = seed + iv.mpc(side, side)
        _, d_z, d_t, *_ = derivatives(square, theta)
        zp = -d_t/d_z
        slope = 2*(zp.real*square.real+zp.imag*square.imag)
        require(direction*slope > 0)
    # Bound the radial derivative lambda=-z F''/(2F') at the true minimiser.
    theta = iv.mpf([endpoints[0][0], endpoints[1][0]])
    tiny = iv.mpf('1e-14')
    val, first, *_ = derivatives(central, theta)
    require(abs(val)+3*tiny**2/2 < abs(first)*tiny)
    side = iv.mpf(['-1e-14','1e-14'])
    square = central + iv.mpc(side, side)
    _, first, _, second, *_ = derivatives(square, theta)
    correction = (-square*second/(2*first)).real
    require(correction in iv.mpf(['-0.2764138907713','-0.2764138907709']))
    print('Radial coefficient Re(lambda) is between -0.2764138907713 and -0.2764138907709.')
    print('Exactly one conjugate pair of phases attains W_0.')
    print('Positive minimising argument: 2.39324227413489 < theta < 2.39324227413490 radians.')


def upper_bound():
    theta = iv.mpf('2.3932422741348965725917153882812424')
    multiplier = iv.mpc(iv.cos(theta), -iv.sin(theta))
    z = iv.mpc('-0.733052523296829803468871187856', '0.083127481904870509963603546822')
    eps = iv.mpf('1e-23')
    require(abs(z)+eps < iv.mpf(UPPER))
    f = sum(multiplier**BINS[k] * z**k / factorial(k) for k in range(DEGREE+1))
    # F'(z)=F(multiplier*z), so the same tail bound applies to F'.
    fp = sum(multiplier**BINS[k] * (multiplier*z)**k/factorial(k)
             for k in range(DEGREE+1))
    # |F''| <= e on the disk |z| <= 1. Linear Taylor approximation and Rouche.
    require((abs(f) + TAIL + 3*eps**2/2 < (abs(fp)-TAIL)*eps))
    print('A simple zero is within 10^-23 of -0.733052523296829803468871187856 + 0.083127481904870509963603546822i.')


if __name__ == '__main__':
    upper_bound()
    cover()
    cover(band=True)
    minimiser_uniqueness()
    print('Certified enclosure for W_0:', LOWER, '< W_0 <', UPPER)

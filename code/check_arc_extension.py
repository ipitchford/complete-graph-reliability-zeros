#!/usr/bin/env python3
"""Certified contour family for the deformed exponential on the phase arc [pi/30, pi/7].

F_w(z) = sum_{k>=0} w^{-k(k-1)/2} z^k / k!,  w = exp(i*theta).

For each sub-arc [a*pi, b*pi] with radius R the script proves, in Arb ball arithmetic
(python-flint), that F_w has no zero on |z| = R for any theta in the sub-arc, using a
first-order Taylor test on (theta, phi)-boxes with an explicit second-order remainder.
At each junction between neighbouring sub-arcs it proves that F_w has no zero in
the closed annulus between the two radii at the junction phase (first-order test in log z).
Option --negative-controls checks that deliberately wrong inputs are rejected.  The junction at
theta = pi/7 is matched to the contour |z| = 3/2 certified on [pi/7, pi] by
check_arc.py.  Since the zero count inside a contour is locally constant in theta
while no zero meets the contour, every phase in [pi/30, pi] has exactly one zero,
counted with multiplicity, inside its contour (conjugation gives the mirror arc).

The script fails (non-zero exit, no PASS line) unless every inequality holds.
It does not use assert statements, so python -O does not weaken it.
"""
import json
import sys
import time
from fractions import Fraction as Q

from flint import acb, arb, ctx

ctx.prec = 128
DEGREE = 60
PI = arb.pi()

# (theta_low/pi, theta_high/pi, radius); ordered from theta = pi/7 downwards.
SUBARCS = [
    (Q(1, 9), Q(1, 7), Q(37, 20)),
    (Q(1, 11), Q(1, 9), Q(21, 10)),
    (Q(1, 13), Q(1, 11), Q(47, 20)),
    (Q(1, 15), Q(1, 13), Q(13, 5)),
    (Q(1, 17), Q(1, 15), Q(57, 20)),
    (Q(1, 19), Q(1, 17), Q(63, 20)),
    (Q(1, 20), Q(1, 19), Q(33, 10)),
    (Q(1, 22), Q(1, 20), Q(7, 2)),
    (Q(1, 24), Q(1, 22), Q(15, 4)),
    (Q(1, 26), Q(1, 24), Q(4, 1)),
    (Q(1, 28), Q(1, 26), Q(17, 4)),
    (Q(1, 30), Q(1, 28), Q(91, 20)),
]
TOP_RADIUS = Q(3, 2)     # contour of check_arc.py on [pi/7, pi]
MAX_DEPTH = 80
BINOM2 = [k * (k - 1) // 2 for k in range(DEGREE + 1)]
INV_FACT = [arb(1) / arb.fac_ui(k) for k in range(DEGREE + 1)]


class CertificateFailure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise CertificateFailure(message)


def tail_bound(radius):
    r = arb(radius.numerator) / radius.denominator
    require(2 * float(radius) <= DEGREE + 2, "tail ratio")
    # sum_{k>D} r^k/k! <= 2 r^(D+1)/(D+1)! when r/(D+2) <= 1/2
    return 2 * r ** (DEGREE + 1) / arb.fac_ui(DEGREE + 1)


def first_order_excludes(theta_c, phi_c, h_theta, h_phi, powers, tail):
    """Exclude zeros of F on the box |theta-theta_c| <= h_theta, |phi-phi_c| <= h_phi, |z| = R.

    Each term is c_k e^{i alpha_k} with alpha_k = k*phi - N_k*theta.  Around the centre,
    e^{i(alpha+d)} = e^{i alpha}(1 + i d + E) with |E| <= d^2/2 and |d| <= N_k h_theta + k h_phi.
    So F = F_c + F_theta*dtheta + F_phi*dphi + E_total with |E_total| <= rem + tail.
    Three sufficient tests for F != 0 on the box: the modulus, each real component, and the
    cross product with F_phi (which removes the first-order phi term).
    """
    f = acb(0)
    f_t = acb(0)
    f_p = acb(0)
    rem = arb(0)
    for k in range(DEGREE + 1):
        b = BINOM2[k]
        term = powers[k] * (acb(0, 1) * (k * phi_c - b * theta_c)).exp()
        f += term
        f_t += acb(0, -b) * term
        f_p += acb(0, k) * term
        rem += powers[k] * (b * h_theta + k * h_phi) ** 2 / 2
    err = rem + tail
    if bool((abs(f) - h_theta * abs(f_t) - h_phi * abs(f_p)).lower() > err.upper()):
        return True
    for part in ("real", "imag"):
        value = abs(getattr(f, part)) - h_theta * abs(getattr(f_t, part)) - h_phi * abs(getattr(f_p, part))
        if bool(value.lower() > err.upper()):
            return True
    cross_c = f.real * f_p.imag - f.imag * f_p.real
    cross_t = f_t.real * f_p.imag - f_t.imag * f_p.real
    value = abs(cross_c) - h_theta * abs(cross_t)
    return bool(value.lower() > (err * abs(f_p)).upper())


def cover_contour(a, b, radius):
    """No zero of F_w on |z| = radius for theta in [a*pi, b*pi] (first-order box test)."""
    r = arb(radius.numerator) / radius.denominator
    tail = tail_bound(radius)
    powers = [r ** k * INV_FACT[k] for k in range(DEGREE + 1)]
    lo = arb(a.numerator) / a.denominator * PI
    hi = arb(b.numerator) / b.denominator * PI
    # boxes are (theta_c, h_theta, phi_c, h_phi); centres and half-widths are exact arb values
    # heuristic weights only (do not affect rigour): typical size of the theta- and phi-derivatives
    sens_theta = sum(BINOM2[k] * float(powers[k].mid()) for k in range(DEGREE + 1))
    sens_phi = sum(k * float(powers[k].mid()) for k in range(DEGREE + 1))
    pieces = 1
    width = (hi - lo) / pieces
    stack = [(lo + (2 * j + 1) * width / 2, width / 2, arb(0), PI, 0) for j in range(pieces)]
    boxes = deepest = 0
    while stack:
        tc, ht, pc, hp, depth = stack.pop()
        if first_order_excludes(tc, pc, ht, hp, powers, tail):
            boxes += 1
            deepest = max(deepest, depth)
            continue
        require(depth < MAX_DEPTH, f"contour |z|={radius} not certified near theta={float(tc.mid())}")
        # split the direction with the larger first-order sensitivity (two children cover the box)
        if float(ht.mid()) * sens_theta > float(hp.mid()) * sens_phi:
            for st in (-1, 1):
                stack.append((tc + st * ht / 2, ht / 2, pc, hp, depth + 1))
        else:
            for sp in (-1, 1):
                stack.append((tc, ht, pc + sp * hp / 2, hp / 2, depth + 1))
    return boxes, deepest


def cover_annulus(theta_frac, r_low, r_high):
    """No zero of F_w with r_low <= |z| <= r_high at theta = theta_frac*pi.

    Write z = exp(u + i*phi).  Term k is e^{-i N_k theta} e^{k(u + i phi)}/k!, so on a box
    |u-u_c| <= h_u, |phi-phi_c| <= h_phi with x = k((u-u_c) + i(phi-phi_c)), |x| <= k(h_u+h_phi):
    e^{x} = 1 + x + E with |E| <= |x|^2 e^{|x|}/2.  Hence F = F_c + G*zeta + E_total with
    G = sum k c_k, |zeta| <= sqrt(h_u^2 + h_phi^2) and |E_total| <= rem + tail.
    """
    theta = arb(theta_frac.numerator) / theta_frac.denominator * PI
    lo, hi = min(r_low, r_high), max(r_low, r_high)
    tail = tail_bound(hi)
    u0 = (arb(lo.numerator) / lo.denominator).log()
    u1 = (arb(hi.numerator) / hi.denominator).log()
    phases = [(acb(0, -BINOM2[k]) * theta).exp() * INV_FACT[k] for k in range(DEGREE + 1)]
    stack = [((u0 + u1) / 2, (u1 - u0) / 2, arb(0), PI, 0)]
    boxes = deepest = 0
    while stack:
        uc, hu, pc, hp, depth = stack.pop()
        zc = (acb(uc, pc)).exp()
        f = acb(0)
        g = acb(0)
        rem = arb(0)
        zk = acb(1)
        for k in range(DEGREE + 1):
            ck = phases[k] * zk
            f += ck
            g += k * ck
            x = k * (hu + hp)
            rem += abs(ck) * x * x * x.exp() / 2
            zk *= zc
        bound = abs(g) * (hu * hu + hp * hp).sqrt() + rem + tail
        if bool(abs(f).lower() > bound.upper()):
            boxes += 1
            deepest = max(deepest, depth)
            continue
        require(depth < MAX_DEPTH, f"annulus at theta={theta_frac}*pi not certified")
        if float(hu.mid()) > float(hp.mid()):
            for su in (-1, 1):
                stack.append((uc + su * hu / 2, hu / 2, pc, hp, depth + 1))
        else:
            for sp in (-1, 1):
                stack.append((uc, hu, pc + sp * hp / 2, hp / 2, depth + 1))
    return boxes, deepest


def main():
    started = time.time()
    record = {"degree": DEGREE, "precision_bits": ctx.prec, "subarcs": [], "junctions": []}
    # Sub-arcs must tile [THETA0, pi/7] without gaps.
    require(SUBARCS[0][1] == Q(1, 7), "top of the family must be pi/7")
    for upper, lower in zip(SUBARCS, SUBARCS[1:]):
        require(upper[0] == lower[1], "sub-arcs must abut")
    previous_radius = TOP_RADIUS
    for a, b, radius in SUBARCS:
        jb, jd = cover_annulus(b, previous_radius, radius) if radius != previous_radius else (0, 0)
        record["junctions"].append({"theta_over_pi": str(b), "radii": [str(previous_radius), str(radius)],
                                    "boxes": jb, "max_depth": jd})
        boxes, deepest = cover_contour(a, b, radius)
        record["subarcs"].append({"theta_over_pi": [str(a), str(b)], "radius": str(radius),
                                  "boxes": boxes, "max_depth": deepest})
        print(f"no zero on |z|={radius} for {a}*pi <= theta <= {b}*pi: {boxes} boxes, depth {deepest}", flush=True)
        previous_radius = radius
    theta0 = SUBARCS[-1][0]
    record["theta0_over_pi"] = str(theta0)
    record["arc_mass"] = str(1 - theta0)
    record["seconds"] = round(time.time() - started, 1)
    print(f"PASS: one simple zero inside the certified contour for every theta in [{theta0}*pi, pi];"
          f" arc mass {1 - theta0}")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as handle:
            json.dump(record, handle, indent=1)


def negative_controls():
    """Both tests must FAIL on inputs where a zero lies on the circle or in the annulus."""
    failures = 0
    try:
        cover_contour(Q(1, 9), Q(1, 7), Q(8, 5))      # the zero rho(theta) crosses |z| = 8/5 on this arc
    except CertificateFailure:
        failures += 1
    try:
        cover_annulus(Q(1, 7), Q(13, 10), Q(3, 2))    # rho(pi/7) = 1.398... lies in this annulus
    except CertificateFailure:
        failures += 1
    if failures != 2:
        print("FAIL: a negative control was certified")
        sys.exit(1)
    print("NEGATIVE CONTROLS PASS: both deliberately wrong inputs were rejected")


if __name__ == "__main__":
    if "--negative-controls" in sys.argv:
        negative_controls()
        sys.exit(0)
    try:
        main()
    except CertificateFailure as failure:
        print("FAIL:", failure)
        sys.exit(1)

#!/usr/bin/env python3
"""Figure: zeros of R_n for K_n (n = 24 and 40) against the predicted inner curve.

Floating illustration only.  Zeros of the nontrivial factor H_n are found with
Arb (python-flint) polynomial root isolation, which returns certified enclosures;
the predicted curve is |q| = exp(-b_n(w)) with b_n(w) = (log (n-1)! - n log|rho|)/(N + n Re lambda),
where rho(theta) is the zero of F_w continued from -pi/4 at theta = pi.
Writes figures/profile.pdf, figures/profile.png and figures/profile-data.json.
"""
import json
import math
import sys
from pathlib import Path

from flint import acb, arb, acb_poly, fmpz_poly, ctx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("figures")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 60


def reliability(n):
    """Integer coefficients of R_n in ascending powers of q (component recurrence)."""
    values = [None, [1]]
    for m in range(2, n + 1):
        coefficients = [1] + [0] * (m * (m - 1) // 2)
        for k in range(1, m):
            cut, ways = k * (m - k), math.comb(m - 1, k - 1)
            for degree, value in enumerate(values[k]):
                coefficients[cut + degree] -= ways * value
        values.append(coefficients)
    return values[n]


def nontrivial_zeros(n):
    poly = fmpz_poly(reliability(n))
    one_minus_q = fmpz_poly([1, -1])
    for _ in range(n - 1):
        poly, remainder = divmod(poly, one_minus_q)
        if remainder != 0:
            raise ArithmeticError("expected (1-q)^(n-1) factor")
    ctx.prec = 2000
    roots = acb_poly(poly).roots(maxprec=40000)
    return [complex(float(r.real.mid()), float(r.imag.mid())) for r in roots]


def rho_and_lambda(theta, guess):
    """Newton on the degree-DEG truncation of F_w; returns rho, lambda = -rho F''/(2F')."""
    w_inv = complex(math.cos(theta), -math.sin(theta))
    coeff = [w_inv ** (k * (k - 1) // 2) / math.factorial(k) for k in range(DEG + 1)]
    z = guess
    for _ in range(60):
        f = sum(c * z ** k for k, c in enumerate(coeff))
        fp = sum(k * c * z ** (k - 1) for k, c in enumerate(coeff) if k)
        step = f / fp
        z -= step
        if abs(step) < 1e-15:
            break
    fp = sum(k * c * z ** (k - 1) for k, c in enumerate(coeff) if k)
    fpp = sum(k * (k - 1) * c * z ** (k - 2) for k, c in enumerate(coeff) if k > 1)
    return z, -z * fpp / (2 * fp)


def main():
    thetas = [math.pi - j * (math.pi - math.pi / 30) / 400 for j in range(401)]
    branch, guess = [], complex(-math.pi / 4, 0)
    for theta in thetas:
        rho, lam = rho_and_lambda(theta, guess)
        branch.append((theta, rho, lam))
        guess = rho
    data = {"branch": [[t, abs(r), l.real] for t, r, l in branch], "orders": {}}
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), sharey=False)
    teal, berry, muted = "#00796b", "#99375c", "#77717c"
    for ax, n in zip(axes, (24, 40)):
        N = n * (n - 1) // 2
        zeros = nontrivial_zeros(n)
        upper = [z for z in zeros if z.imag > 1e-12]
        predicted = []
        for theta, rho, lam in branch:
            b = (math.lgamma(n) - n * math.log(abs(rho))) / (N + n * lam.real)
            predicted.append((theta, math.exp(-b)))
        data["orders"][str(n)] = {"zeros_upper_half": [[math.atan2(z.imag, z.real), abs(z)] for z in upper],
                                  "predicted": predicted}
        ax.scatter([math.atan2(z.imag, z.real) for z in upper], [abs(z) for z in upper],
                   s=7, color=teal, label="zeros of $R_n$")
        ax.plot([p[0] for p in predicted], [p[1] for p in predicted], color=berry, lw=1.1,
                label="predicted inner curve")
        ax.axvline(math.pi / 30, color=muted, lw=0.7, ls=":")
        ax.set_title(f"$K_{{{n}}}$", fontsize=10)
        ax.set_xlabel("argument of $q$ (radians)", fontsize=9)
        ax.set_xlim(0, math.pi)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("modulus of $q$", fontsize=9)
    axes[0].legend(fontsize=7, frameon=False, loc="lower left")
    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "profile.pdf", bbox_inches="tight")
    fig.savefig(OUT / "profile.png", dpi=200, bbox_inches="tight")
    with open(OUT / "profile-data.json", "w") as handle:
        json.dump(data, handle)
    print("profile figure written")


if __name__ == "__main__":
    main()

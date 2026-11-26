#!/usr/bin/env python3
"""Finite-size check of the local zero expansion near the negative real axis.

For q = -exp(-b_n - t/N) the zeros of R_n near the model zeros t0 of 1 + (-1)^d e^{-t}
are computed by a secant iteration on the exact component recurrence at high
precision, then compared with

  first order : t0 - (a^2/2) n b_n^2 + (-1)^(n-1) n e^{-(n-1) b_n}
  refined     : t0 - (a^2/2) n (b_n + t0/N)^2 - (2a^2/3) n b_n^3 + (-1)^(n-1) n e^{-(n-1) b_n}

with a = pi/4.  The refined remainder is O(n^-2).  Writes JSON with every value.
Floating diagnostics only (mpmath, 120 digits); residuals of the computed roots are
recorded and must be below 1e-60 or the script exits with status 1.
"""
import json
import sys
import mpmath as mp

mp.mp.dps = 120


def evaluate(n, q):
    powers = [mp.mpf(1)]
    for _ in range(n * n // 4 + 1):
        powers.append(powers[-1] * q)
    values = [None, mp.mpf(1)]
    for i in range(2, n + 1):
        ways, terms = 1, []
        for k in range(1, i):
            terms.append(ways * powers[k * (i - k)] * values[k])
            ways = ways * (i - k) // k
        values.append(1 - mp.fsum(terms))
    return values[-1]


def main(out):
    a = mp.pi / 4
    rows, failures = [], []
    cases = [(n, 0) for n in (15, 16, 31, 32, 63, 64, 127, 128, 255, 256)] + [(64, 1), (128, 1), (256, 1)]
    for n, k in cases:
        N = n * (n - 1) // 2
        d = N - n + 1
        b = (mp.loggamma(n) - n * mp.log(a)) / N
        t0 = (2 * k * mp.pi * 1j) if d % 2 == 1 else ((2 * k + 1) * mp.pi * 1j)
        if d % 2 == 0 and k == 0:
            continue
        f = lambda t: evaluate(n, -mp.exp(-b - t / N))
        start = t0 if k else mp.mpf(0)
        root = mp.findroot(f, (start, start + mp.mpf('0.05')), solver='secant', tol=mp.mpf('1e-90'))
        residual = abs(f(root))
        if residual > mp.mpf('1e-60'):
            failures.append((n, k))
        sign_term = (-1) ** (n - 1) * n * mp.exp(-(n - 1) * b)
        first = t0 - a**2 / 2 * n * b**2 + sign_term
        refined = t0 - a**2 / 2 * n * (b + t0 / N) ** 2 - 2 * a**2 / 3 * n * b**3 + sign_term
        rows.append({
            "n": n, "n_mod_4": n % 4, "model_zero_t0": mp.nstr(t0, 12), "zero_t": mp.nstr(root, 15),
            "q": mp.nstr(-mp.exp(-b - root / N), 15),
            "first_order_prediction": mp.nstr(first, 12), "first_order_error": mp.nstr(abs(root - first), 6),
            "refined_prediction": mp.nstr(refined, 12), "refined_error": mp.nstr(abs(root - refined), 6),
            "n2_times_refined_error": mp.nstr(n * n * abs(root - refined), 6),
            "recurrence_residual": mp.nstr(residual, 3)})
        print(n, k, mp.nstr(root, 12), mp.nstr(abs(root - first), 4), mp.nstr(abs(root - refined), 4),
              mp.nstr(n * n * abs(root - refined), 4))
    with open(out, "w") as handle:
        json.dump({"a": "pi/4", "digits": mp.mp.dps, "rows": rows}, handle, indent=1)
    if failures:
        print("FAIL: residual too large for", failures)
        sys.exit(1)
    print("PASS: all secant roots have recurrence residual below 1e-60")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "finite_size.json")

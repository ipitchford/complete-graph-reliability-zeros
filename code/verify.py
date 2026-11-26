#!/usr/bin/env python3
"""Independent finite and asymptotic checks for complete-graph reliability.

Run with Python 3, mpmath, and sympy (only used on the Sturm fallback).
All polynomial coefficients in the exact section are integers in ascending order.
"""
from functools import lru_cache
from itertools import combinations
from math import comb
import mpmath as mp


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)



def edge_count(n):
    return n * (n - 1) // 2


def reliability_polynomials(limit):
    """Condition on the component of a fixed vertex, with q the failure probability."""
    values = [None, [1]]
    for n in range(2, limit + 1):
        coefficients = [1] + [0] * edge_count(n)
        for k in range(1, n):
            cut = k * (n - k)
            ways = comb(n - 1, k - 1)
            for degree, value in enumerate(values[k]):
                coefficients[cut + degree] -= ways * value
        values.append(coefficients)
    return values


def enumerate_graphs(n):
    edges = list(combinations(range(n), 2))
    connected_by_size = [0] * (len(edges) + 1)
    for edge_mask in range(1 << len(edges)):
        parent = list(range(n))
        def representative(i):
            while i != parent[i]:
                i = parent[i]
            return i
        for i, (first, second) in enumerate(edges):
            if edge_mask & (1 << i):
                parent[representative(first)] = representative(second)
        if all(representative(i) == representative(0) for i in range(n)):
            connected_by_size[edge_mask.bit_count()] += 1
    # Each graph with m edges has probability (1-q)^m q^(N-m).
    coefficients = [0] * (len(edges) + 1)
    for m, count in enumerate(connected_by_size):
        for j in range(m + 1):
            coefficients[len(edges) - m + j] += (-1)**j * comb(m, j) * count
    return coefficients, sum(connected_by_size)


def variations(poly):
    signs = [c > 0 for c in poly if c]
    return sum(a != b for a, b in zip(signs, signs[1:]))


@lru_cache(maxsize=24)
def pascal(degree):
    return [[comb(n, k) for k in range(n + 1)] for n in range(degree + 1)]


def shift_one(poly):
    # p(t+1) parametrises t_original > 1.
    d = len(poly) - 1
    nums = pascal(d)
    return [sum(poly[j] * nums[j][i] for j in range(i, d + 1)) for i in range(d + 1)]


def interval_zero_one(poly):
    # (1+t)^d p(t/(1+t)) parametrises 0 < t_original < 1.
    return shift_one(poly[::-1])[::-1]


def sturm_positive(poly):
    import sympy as sp
    x = sp.Symbol('x')
    p = sp.Poly.from_list(poly[::-1], x, domain=sp.ZZ)
    return sum(mult * factor.count_roots(0, sp.oo) for factor, mult in p.sqf_list()[1])


def exact_positive_count(poly, stats, depth=0):
    """Descartes certificates count multiplicities. Split at 1 if inconclusive."""
    while poly[0] == 0:
        poly = poly[1:]  # endpoint 0 is excluded
    while poly[-1] == 0:
        poly = poly[:-1]
    stats[0] += 1
    var = variations(poly)
    if var <= 1:
        return var
    if depth == 12:
        stats[1] += 1
        return sturm_positive(poly)
    middle = shift_one(poly)
    mult = next(k for k, c in enumerate(middle) if c)
    return (exact_positive_count(interval_zero_one(poly), stats, depth + 1)
            + mult + exact_positive_count(middle, stats, depth + 1))


def evaluate(n, q):
    """Analytic recurrence at a high-precision complex argument."""
    powers = [mp.mpf(1)]
    for j in range(n*n//4):
        powers.append(powers[-1] * q)
    values = [None, mp.mpf(1)]
    for i in range(2, n + 1):
        # The varying binomial coefficient is updated with exact integer arithmetic.
        ways = 1
        terms = []
        for k in range(1, i):
            terms.append(ways * powers[k*(i-k)] * values[k])
            ways = ways * (i-k) // k
        values.append(1 - mp.fsum(terms))
    return values[-1]


def centre(n):
    a = mp.pi / 4
    return (mp.loggamma(n) - n * mp.log(a)) / edge_count(n)


def local(n, t):
    return evaluate(n, -mp.exp(-centre(n) - t/edge_count(n)))


def fmt(x):
    return mp.nstr(x, 11)


def numerical_sections():
    print('\nTransition window (t as in the paper)')
    print('n  sign   R_n(t=0)      |R_n(1+i) - (1+sign*exp(-1-i))|  cross-precision difference')
    for n in (15, 16, 31, 32, 33, 64, 65, 128, 129, 256, 257):
        sign = (-1)**edge_count(n-1)
        t = mp.mpc(1, 1)
        with mp.workdps(80):
            low = [local(n, v) for v in (0, t)]
        with mp.workdps(150):
            high = [local(n, v) for v in (0, t)]
            discrepancy = max(abs(a-b) for a,b in zip(high, low))
            require(discrepancy < mp.mpf('1e-60'))
            print(n, sign, fmt(high[0]), fmt(abs(high[1] - (1 + sign * mp.exp(-t)))), fmt(discrepancy))
    print('\nNegative-real zero near -1 (high precision secant iteration)')
    print('n  t=N(s-s*)   q                 absolute recurrence residual')
    with mp.workdps(90):
        for n in (15, 16, 32, 64, 128, 256):
            root = mp.findroot(lambda t: local(n, t), (0, mp.mpf('0.1')),
                               solver='secant', tol=mp.mpf('1e-55'))
            q = -mp.exp(-centre(n) - root/edge_count(n))
            with mp.workdps(150):
                residual = abs(evaluate(n, q))
            require(residual < mp.mpf('1e-50'))
            print(n, fmt(root), fmt(q), fmt(residual))


def main():
    limit = 24
    polys = reliability_polynomials(limit)
    print('Complete graph: exact enumeration of labelled connected spanning subgraphs')
    print('n    possible graphs    connected graphs       recurrence equality')
    for n in range(1, 7):
        coefficients, connected = enumerate_graphs(n)
        require(coefficients == polys[n])
        print(n, 2**edge_count(n), connected, 'exact')
    print('\nExact root counts, with multiplicity, certified by integer Descartes subdivisions')
    print('n   (-infinity,-1)    at -1     (-1,0)   nodes   Sturm fallbacks')
    for n in range(2, limit+1):
        poly = [c * (-1)**j for j,c in enumerate(polys[n])]
        stats = [0,0]
        outside = exact_positive_count(shift_one(poly), stats)
        inside = exact_positive_count(interval_zero_one(poly), stats)
        endpoint = int(sum(poly) == 0)
        print(n, outside, endpoint, inside, stats[0], stats[1])
        require(outside == endpoint == 0)
        require(inside == int(n % 4 in (0, 3)))
    numerical_sections()


if __name__ == '__main__':
    main()

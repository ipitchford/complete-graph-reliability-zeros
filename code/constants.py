#!/usr/bin/env python3
"""Derived constants quoted in the text (floating, 50 digits) from the certified enclosures."""
import json, sys
import mpmath as mp
mp.mp.dps = 50
a = mp.pi / 4
A_lo, A_hi = mp.mpf('0.73775075747822433247'), mp.mpf('0.73775075747822433248')
th_lo, th_hi = mp.mpf('2.39324227413489'), mp.mpf('2.39324227413490')
out = {
    "a=pi/4": mp.nstr(a, 20),
    "H_negative_axis=-2(1+log a)": mp.nstr(-2 * (1 + mp.log(a)), 12),
    "H_inner=-2(1+log A) lower": mp.nstr(-2 * (1 + mp.log(A_hi)), 15),
    "H_inner=-2(1+log A) upper": mp.nstr(-2 * (1 + mp.log(A_lo)), 15),
    "theta_star_degrees_lo": mp.nstr(th_lo * 180 / mp.pi, 10),
    "theta_star_degrees_hi": mp.nstr(th_hi * 180 / mp.pi, 10),
    "Jensen_sum_I0(4)=sum 4^k/(k!)^2": mp.nstr(mp.besseli(0, 4), 10),
    "(2/a)^4": mp.nstr((2 / a) ** 4, 10),
    "2/a": mp.nstr(2 / a, 10),
    "log 2": mp.nstr(mp.log(2), 10),
}
for n in (100, 1000, 10**6):
    hn = mp.exp(-(2 * mp.log(n) - 2 * (1 + mp.log(A_lo))) / n)
    out[f"leading-order h_n at n={n}"] = mp.nstr(hn, 8)
for k, v in out.items():
    print(k, v)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "constants.json", "w"), indent=1)

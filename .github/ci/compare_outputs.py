# SPDX-License-Identifier: MIT
"""Compare a replay copy with the recorded outputs under the declared policy:
- results/*.txt identical byte for byte (except status.txt and environment.txt, which record run times and the interpreter);
- results/finite_size.json and results/constants.json identical;
- results/arc_extension.json identical apart from the key "seconds" (run time);
- figures/profile-data.json: same structure; integers, strings and booleans identical; float64 plotting values agree to
  relative 1e-9 (last-digit differences across platforms allowed). Figure images are regenerated but not compared.
usage: compare_outputs.py REPLAY_ROOT RECORDED_ROOT"""
import json, math, sys
from pathlib import Path
new, old = Path(sys.argv[1]), Path(sys.argv[2]); bad = []
for p in sorted((old / 'results').glob('*.txt')):
    if p.name in ('status.txt', 'environment.txt'): continue
    q = new / 'results' / p.name
    if not q.is_file() or q.read_bytes() != p.read_bytes(): bad.append(p.name)
for n in ('finite_size.json', 'constants.json'):
    if (new / 'results' / n).read_bytes() != (old / 'results' / n).read_bytes(): bad.append(n)
def drop(x):
    if isinstance(x, dict): return {k: drop(v) for k, v in x.items() if k != 'seconds'}
    if isinstance(x, list): return [drop(v) for v in x]
    return x
if drop(json.loads((new / 'results/arc_extension.json').read_text())) != drop(json.loads((old / 'results/arc_extension.json').read_text())): bad.append('arc_extension.json')
def close(a, b, path='$'):
    if type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool)):
        return [f'{path}: type']
    if isinstance(a, dict):
        if a.keys() != b.keys(): return [f'{path}: keys']
        return [e for k in a for e in close(a[k], b[k], f'{path}.{k}')]
    if isinstance(a, list):
        if len(a) != len(b): return [f'{path}: length']
        return [e for i, (x, y) in enumerate(zip(a, b)) for e in close(x, y, f'{path}[{i}]')]
    if isinstance(a, float) or isinstance(b, float):
        return [] if math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-300) else [f'{path}: {a} vs {b}']
    return [] if a == b else [f'{path}: {a!r} vs {b!r}']
e = close(json.loads((new / 'figures/profile-data.json').read_text()), json.loads((old / 'figures/profile-data.json').read_text()))
if e: bad.append(f'profile-data.json {e[:5]}')
if bad: print('DIFFERENCES', bad); sys.exit(1)
print('OK all recorded outputs reproduced under the declared policy')

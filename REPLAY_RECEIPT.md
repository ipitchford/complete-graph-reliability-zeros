# Replay receipt (v0.1.0-candidate, 26 November 2026)

Release replay by the publisher from a fresh copy of the approved candidate's code/ (outputs regenerated into empty
results/ and figures/ folders), working directory code/, macOS arm64. Interpreter as printed inside the run: CPython
3.12.11, python-flint 0.7.1, mpmath 1.3.0 (see [ENVIRONMENT.txt](ENVIRONMENT.txt)).

| Program | Flags | Exit | Output vs recorded |
|---|---|---|---|
| verify.py | -E -s -B -O | 0 | identical |
| whittaker_certificate.py | -E -s -B -O | 0 | identical |
| whittaker_certificate.py --mpmath | -E -s -B -O | 0 | identical |
| check_intervals_independent.py | -E -s -B -O | 0 | identical |
| check_arc.py | -E -s -B -O | 0 | identical |
| check_intervals_independent.py --arc | -E -s -B -O | 0 | identical |
| check_arc_extension.py ../results/arc_extension.json | -E -s -B -O | 0 | text identical; JSON identical apart from "seconds" |
| check_arc_extension.py --negative-controls | -E -s -B -O | 0 | identical (both deliberately wrong inputs rejected) |
| finite_size.py ../results/finite_size.json | -E -s -B -O | 0 | text and JSON identical |
| constants.py ../results/constants.json | -E -s -B | 0 | text and JSON identical |
| profile_figure.py ../figures | -E -s -B | 0 | text, profile-data.json and profile.png identical; PDF differs only by its embedded date |
| plot_roots.py | -E -s -B | 0 | text and clique-roots.png identical; PDF differs only by its embedded date |
| Negative control | whittaker_certificate.py with the claimed enclosure of W_0 moved to [0.738, 0.738 + 1e-20], -E -s -B -O | 1 | `ArithmeticError: insufficient depth for claimed lower endpoint` |

stderr was empty for every passing run. The lead's earlier clean replay (24-25 November 2026, same interpreter)
produced the same output. After the replay only the title-page DOI line of manuscript.tex, the rebuilt PDF and the
release records changed; code/, results/ and figures/ are byte-identical to the approved candidate.

Comparison policy for other platforms (declared before tagging; also in CI, README.md and ENVIRONMENT.txt): printed
outputs identical (except status.txt and environment.txt); finite_size.json and constants.json identical;
arc_extension.json identical apart from "seconds"; profile-data.json floats to relative 1e-9; images not compared.

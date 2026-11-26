# The innermost reliability zeros of complete graphs and Waldvogel's constant

**Version 0.1.0-candidate. DOI [10.5281/zenodo.24546697](https://doi.org/10.5281/zenodo.24546697)** (all versions:
[10.5281/zenodo.24546696](https://doi.org/10.5281/zenodo.24546696)). Author: Anonymous. Text CC0-1.0; code MIT.
Unrefereed research candidate: the assurance is internal AI-assisted review and computer checks only (see
[REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)).

## What this is about

Take a network in which every link can fail, independently, with the same probability q. The *reliability polynomial*
gives the chance that the surviving links still connect every node. For the complete network K_n, where every pair of
the n nodes is linked, it is the chance that a random graph is connected.

Because it is a polynomial in q, one can plug in complex values of q and ask where it vanishes. These "reliability
zeros" have been studied since Brown and Colbourn conjectured in 1992 that they all lie in the closed unit disk. Royle
and Sokal showed that this is false in general, but only just. For the complete graphs, pictures show the zeros
forming an oval that creeps out towards the unit circle as n grows. Buys proved that the reliability tends to 1 inside
the disk and asked how the zeros are spread round the circle.

## What is new

1. **How close the zeros come to the origin.** The innermost zeros sit at distance about (2 log n - 1.392)/n inside the
   unit circle. The constant -1.392 is -2(1 + log W_0), where W_0 = 0.7377... is a constant Waldvogel studied in
   1983: the smallest possible modulus of a zero of a family of entire functions (the "deformed exponential").
2. **In which direction.** The innermost zeros approach two mirror-image directions, at about plus and minus 137
   degrees (2.393 radians), not the negative real axis. This part is a computer-assisted proof with interval
   arithmetic.
3. **The shape of the curve.** On the arc from pi/30 to 59pi/30, which carries 29/30 of the zeros, the zeros follow an
   explicit curve, one zero per gap of size 2pi/N, where N = n(n-1)/2 is the degree. Figure 1 shows the zeros of K_24
   and K_40 on the predicted curve.
4. **Negative real zeros.** For all large n there is exactly one negative real zero when n leaves remainder 0 or 3 on
   division by 4, and none otherwise. It is located very precisely: its distance from the unit circle is about
   (2 log n - 1.517)/n, so it lies farther out than the innermost complex zeros. Exact computation confirms the count
   for every n up to 24.
5. **A short criterion for matroids.** The zeros of the spanning reliability of loopless matroids spread evenly round
   the unit circle exactly when rank divided by size tends to zero. This answers a question of Buys; the proof rests on
   a theorem of Hughes and Nikeghbali.

**Idea of the proof.** The generating function for connected graphs is cut at half its length. A graph has at most one
component with more than n/2 vertices, so the reliability splits exactly into a "one giant component" term, close to 1
(Gilbert 1959), and an "only small components" term governed by the smallest zero of the deformed exponential function.
A zero needs the two terms to cancel, which pins the radius. The closest approach happens where that smallest zero is
smallest, and that is Waldvogel's constant.

**Limits.** The negative-axis count is proved for all sufficiently large n, with no explicit threshold, and checked
exactly up to n = 24. The curve result leaves out an arc of 1/30 of the circle near q = 1 and says nothing about a
vanishing fraction of the zeros, which could lie anywhere; whether every zero lies inside the unit circle is still
open. The directions, the precise radius formula and the 29/30 arc rest on interval-arithmetic certificates. Classical
inputs are credited in the paper (Gilbert, Sokal, Waldvogel, Hughes and Nikeghbali, Alsmeyer et al.).

## Contents

| Path | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (15 pages). Section 1 is a plain introduction with a roadmap. |
| [figures/](figures/) | Figure 1 (profile.pdf/.png, data profile-data.json) and Figure 2 (clique-roots.pdf/.png, zeros of K_14) |
| [code/verify.py](code/verify.py) | Exact checks: enumeration through K_6, exact real-zero counts for 2 <= n <= 24 |
| [code/whittaker_certificate.py](code/whittaker_certificate.py), [code/certified_intervals.py](code/certified_intervals.py) | Interval certificate for W_0, the direction and the curvature (integer-interval backend; `--mpmath` for the mpmath.iv backend) |
| [code/check_intervals_independent.py](code/check_intervals_independent.py) | Independent re-implementation of the certificates (python-flint Arb balls); `--arc` for the arc check |
| [code/check_arc.py](code/check_arc.py), [code/check_arc_extension.py](code/check_arc_extension.py) | Arc certificates: [pi/7, pi] and the twelve sub-arcs down to pi/30; `--negative-controls` runs two deliberately wrong inputs |
| [code/finite_size.py](code/finite_size.py) | Table 2: zeros from the recurrence at 120 digits against the asymptotic formula |
| [code/constants.py](code/constants.py) | Derived constants (1.39170..., 1.51687..., 137.12 degrees) |
| [code/profile_figure.py](code/profile_figure.py), [code/plot_roots.py](code/plot_roots.py) | The two figures |
| [results/](results/) | Recorded outputs of every program |

Machine-readable records: [AI_INDEX.md](AI_INDEX.md), [RESEARCH_GATES.json](RESEARCH_GATES.json),
[refs.bib](refs.bib), [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json), [MANIFEST.sha256](MANIFEST.sha256).

## Replay

Python 3.12 with the pinned packages in [requirements.txt](requirements.txt). From `code/` (the programs import
their siblings, so use `-E -s` rather than `-I`):

```sh
python3 -E -s -B -O verify.py
python3 -E -s -B -O whittaker_certificate.py            # and with --mpmath
python3 -E -s -B -O check_intervals_independent.py      # and with --arc
python3 -E -s -B -O check_arc.py
python3 -E -s -B -O check_arc_extension.py ../results/arc_extension.json
python3 -E -s -B -O check_arc_extension.py --negative-controls
python3 -E -s -B -O finite_size.py ../results/finite_size.json
python3 -E -s -B constants.py ../results/constants.json
python3 -E -s -B profile_figure.py ../figures
python3 -E -s -B plot_roots.py
```

Each run takes seconds (about one minute in all). Certificate programs raise `CertificateFailure` (or another
`ArithmeticError`) on any failed inequality and do not rely on `assert`, so they also fail correctly under `-O`.

**Comparison policy** (declared before tagging; the same in CI, [ENVIRONMENT.txt](ENVIRONMENT.txt) and
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)): printed outputs `results/*.txt` identical (except `status.txt` and
`environment.txt`, which record run times and the interpreter); `finite_size.json` and `constants.json` identical;
`arc_extension.json` identical apart from the key `seconds`; the float64 plotting values in
`figures/profile-data.json` agree to relative 1e-9. Figure images are regenerated but not compared byte for byte.
GitHub CI ([.github/workflows/replay.yml](.github/workflows/replay.yml)) runs every program on Linux, applies this
policy with [.github/ci/compare_outputs.py](.github/ci/compare_outputs.py) and runs a negative control.

## Citation

Anonymous (2026). *The innermost reliability zeros of complete graphs and Waldvogel's constant*, version
0.1.0-candidate. Zenodo. https://doi.org/10.5281/zenodo.24546697

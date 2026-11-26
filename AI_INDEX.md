# AI index: The innermost reliability zeros of complete graphs and Waldvogel's constant (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with the problem, what was known, what is new, the idea of the proof and a roadmap).
DOI 10.5281/zenodo.24546697 (concept 10.5281/zenodo.24546696). Status: unrefereed candidate; internal AI-assisted
review and computer checks only.

## Notation

R_n(q) = all-terminal reliability polynomial of K_n (probability that K_n stays connected when each edge fails
independently with probability q); degree N = n(n-1)/2. h_n = smallest modulus of a zero of R_n. F_w(z) = deformed
exponential sum_k w^{k(k-1)/2} z^k / k! with |w| = 1 (Macintyre's function); rho(w) its zero of smallest modulus.
A = W_0 = min over unimodular w of |rho(w)| = 0.73775075747822433247... (Waldvogel 1983, rotation-restricted Whittaker
constant); a = pi/4 = |rho(-1)|. H_n(q) = -n log|q| - 2 log n (rescaled radial coordinate).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A | Zero-counting measures of the spanning reliability of loopless matroids M_j converge to the uniform measure on the unit circle iff rank/size -> 0 (answers Buys) | Proposition 1 | written proof; rests on the deterministic theorem of Hughes and Nikeghbali |
| B | Negative real axis: for all sufficiently large n exactly one zero in (-1,0) when n = 0, 3 (mod 4) and none otherwise; located with error O(n^-4) in q (O(n^-2) on the 1/N scale); no zeros on q <= -1 or q > 0 other than q = 1 | Theorem 7; eq. (9); Table 2 | written proof; uses Alsmeyer et al. Remark 2.2; no explicit threshold; exact count checked for 2 <= n <= 24 (verify.py) |
| C | Inner radius: -n log h_n - 2 log n -> -2(1 + log W_0) = -1.39170152... | Theorem 8 | written proof (Lemmas 2-6: zero-free disk, half-degree split, phase-uniform approximation, isolated smallest zero, Jensen bound) |
| D | Extremal directions: the minimum defining W_0 is attained at exactly one conjugate pair theta = +-theta*, theta* = 2.39324227413489... (about 137.12 degrees); innermost zeros approach these directions; h_n given with absolute error o(n^-2) | Theorem 9 (computer-assisted) | interval certificates (two backends plus an independent Arb re-implementation) |
| E | Joint limit law of (arg q, H_n(q)) on every closed I in (pi/30, 59pi/30) and on [pi/30, 59pi/30] (mass 29/30); one zero per angular gap 2pi/N on the explicit curve | Theorem 10 (computer-assisted); Figure 1 | written proof plus the arc certificates (check_arc.py for [pi/7, pi]; check_arc_extension.py for twelve sub-arcs pi/30..pi/7) |
| F | Negative axis is not extremal: constant -2(1 + log(pi/4)) = -1.51687104946, about 0.125/n farther in modulus | Section 6 | constants.py |

## Evidence

| Check | Program | Output | Arithmetic |
|---|---|---|---|
| Enumeration through K_6 (26704 connected spanning subgraphs); exact real-zero counts for 2 <= n <= 24 | [code/verify.py](code/verify.py) | [results/verify.txt](results/verify.txt) | exact integers; Descartes/Sturm |
| Enclosure of W_0 in (0.73775075747822433247, 0.73775075747822433248); theta* and lambda enclosures; curvature > 0.31 (958 boxes at depth 37, 355 at depth 14) | [code/whittaker_certificate.py](code/whittaker_certificate.py) (+ `--mpmath`) | [results/whittaker_verification.txt](results/whittaker_verification.txt), [results/whittaker_verification_mpmath.txt](results/whittaker_verification_mpmath.txt) | outward-rounded integer intervals; mpmath.iv |
| Independent re-implementation | [code/check_intervals_independent.py](code/check_intervals_independent.py) | [results/independent_interval_verification.txt](results/independent_interval_verification.txt), [results/independent_arc_verification.txt](results/independent_arc_verification.txt) | Arb balls (python-flint) |
| Arc [pi/7, pi], radius 3/2 (365 boxes) | [code/check_arc.py](code/check_arc.py) | [results/arc_verification.txt](results/arc_verification.txt) | integer intervals |
| Twelve sub-arcs pi/30..pi/7, radii 37/20..91/20, junction annuli; PASS, arc mass 29/30; two negative controls rejected | [code/check_arc_extension.py](code/check_arc_extension.py) | [results/arc_extension.txt](results/arc_extension.txt), [.json](results/arc_extension.json), [results/arc_extension_negative_controls.txt](results/arc_extension_negative_controls.txt) | Arb balls |
| Table 2: zeros from the recurrence at 120 digits vs eq. (9); residuals < 1e-60 | [code/finite_size.py](code/finite_size.py) | [results/finite_size.txt](results/finite_size.txt), [.json](results/finite_size.json) | mpmath 120 digits |
| Derived constants (1.39170152316291, 1.51687104946, 137.1226817 degrees, Jensen sums) | [code/constants.py](code/constants.py) | [results/constants.txt](results/constants.txt), [.json](results/constants.json) | mpmath (floating; not a certificate) |
| Figure 1: zeros of K_24 and K_40 against the predicted curve | [code/profile_figure.py](code/profile_figure.py) | [figures/profile.pdf](figures/profile.pdf), [figures/profile-data.json](figures/profile-data.json) | Arb root isolation, float64 plotting |
| Figure 2: 78 nontrivial zeros of K_14 in disjoint certified disks of radius 1e-22 | [code/plot_roots.py](code/plot_roots.py) | [figures/clique-roots.pdf](figures/clique-roots.pdf), [results/plot_roots.txt](results/plot_roots.txt) | Arb root isolation |

Numbers quoted on the release page: W_0 = 0.7377...; -1.392 (= -2(1 + log W_0)); -1.517 (negative axis); about 137
degrees / 2.393 radians; 29/30 (arc pi/30 to 59pi/30); exact check up to n = 24; Figure 1 orders 24 and 40. All are
proved in the paper and reproduced by the programs above.

## Limits

The negative-axis count is eventual (no explicit threshold); checked exactly for 2 <= n <= 24. Theorem 10 omits an arc
of 1/30 of the circle around q = 1 and allows o(N) exceptional zeros anywhere; whether every nontrivial zero lies in
the open unit disk (Buys, Question 4.1; Sokal's conjecture) is open. Theorems 9 and 10 are computer-assisted. Table 2
and constants.py values are high-precision floating evaluations, not certificates.

## Replay and records

Replay commands and comparison policy: [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt),
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI [.github/workflows/replay.yml](.github/workflows/replay.yml) with
[.github/ci/compare_outputs.py](.github/ci/compare_outputs.py). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior
art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Licences:
[LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).

# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (19-25 November 2026)

1. **One combined five-perspective Claude review panel** (editor, methodology and proof checking, domain and
   literature, perspective, and devil's advocate) on two earlier drafts about the same problem, with the cited sources
   read. The panel found the analytic proofs correct as far as line-by-line checking goes and reproduced the key
   constants independently; the work needed was presentation and merging, not mathematics.
2. **One revision round with extensions.** The two drafts were merged into one paper with one notation and a
   context-first introduction; the abstract states exactly what is proved for all large n and what is checked; the
   computer-assisted theorems are labelled; the finite-size formula gained its next term (Table 2); the argument
   excluding zeros with q <= -1 now rests on Alsmeyer et al. instead of a source that was not cited; attributions to
   Gilbert, Chang and Shrock, Brown and DeGagné, Royle and Sokal, and Bencs, Piombi and Regts were added. Extensions: a
   new ball-arithmetic certificate extends the curve result from 6/7 to 29/30 of the circle; a new Figure 1 compares
   the zeros of K_24 and K_40 with the predicted curve. Certificate programs no longer rely on `assert`.
3. **Re-review** of the revision and a **proof audit** of every new or rewritten proof. Two local defects were found
   and fixed (an error order in Theorem 7(2) corrected to O((log n)^2/n); a restored step in the proof of Theorem
   7(3)); the new arc certificate was required to pass in a clean replay before release, and did; wording and citation
   items were applied.
4. **Clean replays** of every program from fresh copies, with negative controls, by the lead and again by the
   publisher (see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)).

No Astra/Codex second-model review took place in this batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| "Two manuscripts with clashing notation; the headline is buried." | One paper; the inner radius leads; Table 1 separates classical inputs from new statements. |
| "The abstract overclaims the negative-axis count." | Stated as "for all sufficiently large n", with the exact check for 2 <= n <= 24; an explicit threshold is an open problem. |
| "The finite-size formula hides a term of the same order." | The cubic term is now explicit (Waldvogel's expansion credited); remainder O(n^-2) on the 1/N scale; Table 2 shows the errors. |
| "Which parts are computer-assisted?" | Theorems 9 and 10 are labelled; Section 8 says which computations are exact, interval-certified or floating, and which backend gives which bound. |
| "The arc of 6/7 is a certificate limit, not a natural one." | Extended to 29/30 by a new certificate; the remaining arc near q = 1 is stated as open with the reason. |
| Not addressed | An explicit threshold for the negative-axis count, the arc near q = 1, and whether every zero lies in the unit disk are listed as open problems in Section 9. |

## Scope of assurance

Internal AI-assisted review only: written proofs checked by AI reviewers, computer checks replayed, sources read. No
human specialist review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
Literature searches were of limited scope (see [RESEARCH_GATES.json](RESEARCH_GATES.json)).

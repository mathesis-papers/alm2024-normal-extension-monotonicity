# Common transversal probabilities

This note answers Problem 5.2 of Aivazidis, Loukaki and Müller,
arXiv:2209.00986v2, p. 37 (journal version, p. 1280): for every finite group
G, subgroup H and normal subgroup N, P_G(H) ≤ P_G(HN). Equality holds
exactly when N ≤ H or H is normal in G.

It also characterizes equality in the homomorphism inequality (5.1), proves
the coprime-order assertion in Problem 5.3, and gives a finite-quotient
probability model for the finite-index setting in Problem 5.4. The S3 chain
with values 1, 1/2, 1 illustrates the limits of arbitrary containment in
Problem 5.7; the source does not assert monotonicity for arbitrary containment.
The general equality classification in Problem 5.3 and Problems 5.1, 5.5,
5.6 and 5.8 are outside the note's coverage.

## Files

- `note.tex`: complete proofs and bibliography.
- `note.pdf`: mathematical note.
- `metadata.json`: identity, creators, scope, license and declaration link.
- `verify.py`: standard-library Python 3 exact finite checks.
- `checks.json`: deterministic output for the explicitly listed groups.
- `SHA256SUMS`: hashes of the five source-package files.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a TeX distribution with the packages named in `note.tex`:

```sh
pdflatex -halt-on-error note.tex
pdflatex -halt-on-error note.tex
```

The bibliography is included in the source; BibTeX is unnecessary.

## Exact finite checks

```sh
python3 verify.py > checks-rebuilt.json
cmp checks.json checks-rebuilt.json
```

The output has `status: PASS` when all assertions hold. The original 24-group
family covers 238 subgroups and 1,628 subgroup/normal pairs, with actual
transversal enumeration, component surjections, both equality criteria,
homomorphism inequalities and quotient identities. Additional checks cover
576 grid-weight pairs, 12 coprime direct products, and the rectangular bound
in S3 × S3. The group list and all finite ranges appear in `checks.json`.
These are specified examples, not an enumeration of all groups below a
given order. Integer multiplication tables and rational arithmetic are used;
there are no additional Python dependencies.

## Scope of the proofs

The declaration `D5.S3.FiniteGroups.TransversalProbabilityMonotone.result`
proves the universal finite-group normal-extension inequality, with actual
two-sided transversal subsets and denominator |H|^[G:H]. Its exact source is:

https://github.com/the-omega-institute/trureturing/blob/13f4288c2940bd2aa4020cec789956b9265769f5/D5/S3/FiniteGroups/TransversalProbabilityMonotone.lean

The structural equality, homomorphism, direct-product and finite-index
statements have complete written proofs in the note. The infinite-group
statement is a reduction through the core of a finite-index subgroup to a
finite quotient. It does not put a uniform distribution on the elements of
an infinite coset. Finite checks support the specified examples, rather than
establishing the universal assertions.

The counting formula is credited to Aivazidis, Loukaki and Müller,
Proposition 3.2. The coset-intersection geometry and the use of finite core
quotients are credited to Button, Chiodo and Zeron-Medina Laris (2013/2014).

## Creators and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. CC BY 4.0 (`cc-by-4.0`).

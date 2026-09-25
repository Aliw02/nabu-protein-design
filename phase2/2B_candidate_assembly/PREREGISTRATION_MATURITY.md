# Phase 2B.1 — Assembly Evidence-Maturity Ablation

## Motivation

Assembly V1 passed engineering constraints at the 5% RhlA bootstrap, but the subset of proposals with historical truth did not contain a Top-1% hit.

This may be caused by insufficient model evidence rather than by the assembly mechanism itself.

At the 5% RhlA bootstrap, the frozen V8.3 hierarchy is still evidence-limited.

## Single changed variable

Keep fixed:

- ReferenceProtein rules;
- mutation vocabulary;
- CandidateAssembler V1;
- target order = 3;
- beam width = 256;
- proposal count = 64;
- frozen V8.3 mathematics;
- post-freeze evaluation protocol.

Change only the amount of measured evidence.

## Evidence levels

Use deterministic identity-only nested seeds:

- 5%;
- 10%;
- 20%.

No label is used to choose seed membership.

## Metrics

For each evidence level:

- proposal scoreable count;
- assembly trace;
- proposal SHA before truth lookup;
- historical oracle match count;
- best true fitness among matched proposals;
- Top-1% hits among matched proposals;
- mean and median historical fitness percentile among matched proposals;
- router mode;
- main/pair/triplet/quartet memory sizes.

## Interpretation

If assembly quality improves materially with evidence maturity, 2C should use a knowledge-maturity gate before enabling assembly.

If assembly quality does not improve, do not hide the result by increasing the budget again. Diagnose or redesign CandidateAssembler before combining it with acquisition.

This is RhlA development evidence only.

# 2B — Candidate Assembly

2B begins after the 2A gate opened with an explicit acquisition limitation.

This stage is intentionally isolated from closed-loop acquisition.

## V1 layers

- `ReferenceProtein` — validates mutation source residue against the reference sequence.
- `MutationVocabulary` — frozen allowed substitutions.
- `MeasurementAggregator` — canonical replicate aggregation before V8.3.
- `CandidateAssembler` — bounded beam expansion from mutation vocabulary to higher-order candidate proposals.

## V1 development POC

The initial RhlA POC:

- uses the deterministic 5% measured bootstrap;
- derives the mutation vocabulary from identities only;
- infers a development reference sequence from consistent source residues in that vocabulary;
- fits frozen V8.3 only on measured labels;
- assembles triple-mutant proposals without reading unmeasured truth;
- freezes and hashes proposals;
- only then checks whether proposals exist in the historical RhlA assay table.

This reference inference is a development convenience.

A real biological deployment must receive an independently supplied reference sequence.

## Boundary

2B does not solve or modify the PhoQ acquisition limitation.

Assembly + acquisition is deferred to 2C.

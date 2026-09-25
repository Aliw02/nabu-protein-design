# Phase 2D-A AAV — Adapter Scope Result

**Status: ADAPTER_SCOPE_FAIL (label-blind).**

The identity-only V2 preflight read no target values.

Pinned AAV source:
- archive SHA256: `ad91ba8d5b390d793fc9393f8003ff7b0290fbe2e13db58cb6ad72bc981a99fd`
- reference SHA256: `9b5e572c2a18d27482b629efeb48f573e866fe61bb6aea58bb8c087e4177fbcb`

Published BO candidate pool:
- total: **59,459**
- representable by frozen NABU token grammar: **27,018 (45.44%)**
- unrepresentable: **32,441 (54.56%)**

The pool contains length-changing variants. Frozen V8.3 mutation identities encode
position-indexed substitutions and deletion targets (`A10C`, `A10*`) but do not
define insertion identities.

The preregistered rule forbids silently dropping unsupported identities.
Therefore the AAV head-to-head is not run on a reduced subset and no AAV
fitness targets are revealed under Blind V1.

This is an external benchmark scope limitation, not a negative fitness result.
A future architecture with frozen indel semantics would require a new phase/version,
not a post-hoc Phase-2D retune.

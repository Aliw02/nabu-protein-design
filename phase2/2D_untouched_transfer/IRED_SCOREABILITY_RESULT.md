# Phase 2D-B IRED — Frozen scoreability handling

This rule was frozen **before target reveal**.

Identity-only result:

- fit/train identities: 3,746
- validation identities held out: 662
- official test identities: 4,178
- test identities with all main mutations supported: 3,648
- test identities with at least one supported internal pair: 133
- frozen V8.3 scoreable identities: **131 / 4,178 = 3.135%**

The primary full-test IRED result is therefore classified as:

> **MODEL_SCOPE_FAIL**

No B2 fallback, score imputation, pair-threshold relaxation, or silent test subsetting
is allowed. Those would change the frozen system after blind identity contact.

A one-time target reveal is still allowed for a clearly labeled diagnostic evaluation
on the exact 131 identity-only scoreable variants. That diagnostic does not count as
a Phase-2D full-test PASS and is not directly compared against published full-test
FLIP2 scores.

Source SHA256:
`aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74`

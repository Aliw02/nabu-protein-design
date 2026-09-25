# Phase 2D Identity Preflight V1 — Audit Note

## Preserved failed preflight

Workflow run: `36163646395`  
Artifact: `10876508881`  
Artifact digest: `sha256:be990151bfd621015f2668f0ac67bd07ed3024cae1798f19ff4e80d52b696891`

V1 read no target values (`target_values_read=false`) and recorded:

- AAV archive SHA256: `ad91ba8d5b390d793fc9393f8003ff7b0290fbe2e13db58cb6ad72bc981a99fd`
- AAV reference SHA256: `9b5e572c2a18d27482b629efeb48f573e866fe61bb6aea58bb8c087e4177fbcb`
- IRED source SHA256: `aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74`

## Label-independent bug

V1 computed the AAV representation gate across every row in `sampled.csv`,
including rows outside the published BO training pool.

The released protein-uq loader defines the BO training pool as:

`set == "train" AND validation is NaN`

and holds `validation == True` rows separately.

Therefore V1's AAV gate did not test the preregistered candidate pool.

## Fix

V2 changes only identity-pool scoping:

- AAV representation coverage is evaluated on the exact published BO train pool.
- Source bytes, reference, frozen NABU token grammar, scientific metrics and
  winner rules are unchanged.
- No target values were read to identify or fix this bug.

The V1 artifact is preserved and is not overwritten.

## V2 contract-test rerun audit

Two automatic V2 workflow attempts failed before any blind source download:

- run `36163922402` at head `794bd03f9cfa67efe635749bb4960d878363e485`;
- run `36163945248` at head `2c92fcf5e24ead81b22d9eb758952f0ce81e9274`.

Failure occurred in the synthetic contract fixture. The fixture encoded normal AAV
train/test rows with `validation=False`, while the released protein-uq loader
uses `validation.isna()` for the BO training pool. Consequently the synthetic
published pool was empty and the representability assertion failed.

This is label-independent and occurred before download steps. The fix changes
only the synthetic fixture to use missing validation values for ordinary rows.
No scientific rule, source byte, model rule, metric, or blind label is changed.

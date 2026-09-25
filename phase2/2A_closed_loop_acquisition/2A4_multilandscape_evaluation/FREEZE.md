# Phase 2A.4 Freeze

## Status

**PASS WITH LIMITATION**

The preregistered multi-landscape evaluation completed successfully on:

- GB1;
- TrpB;
- PhoQ.

Protocol was fixed before evaluation:

- deterministic identity-only pool = 2048;
- bootstrap = 64;
- final budget = 256;
- primary cadence = micro_8;
- sensitivity cadence = micro_16;
- no per-landscape retuning.

## Generalization result

### GB1

At micro_8:

- random AUC = 0.492134;
- historical 50/50 AUC = 0.576258;
- V2 AUC = 0.506032.

Historical 50/50 improved over random.

V2 degraded relative to the historical controller at micro_8.

At micro_16, historical and V2 tied on the primary discovery metrics.

### TrpB

At micro_8:

- random AUC = 0.659724;
- historical 50/50 AUC = 0.882817;
- V2 AUC = 0.882817.

Historical 50/50 and V2 both materially improved over random and tied each other.

The same controller pair also tied at micro_16.

### PhoQ

At micro_8:

- random AUC = 0.565568;
- historical 50/50 AUC = 0.400004;
- V2 AUC = 0.400004.

Random acquisition outperformed both active controllers.

At micro_16:

- historical 50/50 AUC = 0.447095;
- V2 AUC = 0.474971.

V2 improved over historical 50/50, but the landscape-level acquisition failure was not resolved.

## What is established

- closed-loop active acquisition gains reproduce on more than one landscape;
- historical 50/50 is strong on RhlA, GB1 and TrpB, but not robust on PhoQ;
- Relative-Evidence V2 is behaviorally valid but not generally superior;
- controller performance is landscape-sensitive;
- no per-landscape policy switch is allowed.

## Important failure diagnostic

The PhoQ result is a real failure case and must remain visible.

Final scoreable coverage was much lower under active controllers than under random acquisition, but this is a correlation with the failure, not proof of causation.

Do not tune a new controller on PhoQ before 2B merely to erase this failure.

## Gate decision

The 2A gate to **2B Candidate Assembly is OPEN WITH LIMITATION**.

Reason:

- the closed-loop mechanism is leak-safe and reproducible;
- batch staleness is understood;
- gains reproduce on GB1 and TrpB;
- the PhoQ failure is explicit and bounded;
- controller limitations are now known well enough to keep assembly development separate from acquisition.

2B must test candidate assembly as its own mechanism.

Do not claim that acquisition is universally solved.

The assembly + acquisition combination remains deferred to 2C.

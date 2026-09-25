# NABU Phase 2 Workspace

This directory is the single working area for NABU Phase 2.

Start with [ROADMAP.md](ROADMAP.md). Do not skip ahead between sections unless the roadmap explicitly allows it.

## Structure

- `2A_closed_loop_acquisition/`
  - `2A0_campaign_simulator/`
  - `2A1_baselines/`
  - `2A2_micro_batch_loop/`
  - `2A3_adaptive_controller/`
  - `2A4_multilandscape_evaluation/`
- `2B_candidate_assembly/`
- `2C_assembly_plus_acquisition/`
- `2D_untouched_transfer/`

Each section owns its own `results/` directory. Results from one section must not be mixed into another section.

## Frozen dependency

Phase 2 builds above the frozen Phase-1 V8.3 core. The Phase-1 core is a reference baseline and must remain reproducible and untouched.

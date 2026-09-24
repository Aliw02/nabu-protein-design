import argparse
import hashlib
import importlib.util
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_BUDGETS = [5, 10, 20, 40]
HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]

PASSIVE_STAGE_A = (
    ROOT
    / "nabu_v8_3_rhla_sample_efficiency_stage_a"
    / "STAGE_A_MANIFEST.json"
)
PASSIVE_MEMBERSHIP = (
    ROOT
    / "nabu_v8_3_rhla_sample_efficiency_stage_a"
    / "BUDGET_MEMBERSHIP.csv"
)
PASSIVE_PREDICTIONS = (
    ROOT
    / "nabu_v8_3_rhla_sample_efficiency_stage_a"
    / "ALL_BUDGET_PREDICTIONS.csv"
)

HELPER_PATH = (
    ROOT
    / "experiments"
    / "v8_3_rhla_sample_efficiency"
    / "01_stage_a_freeze_budgets.py"
)

spec = importlib.util.spec_from_file_location("rhla_passive_helper", HELPER_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not import helper from {HELPER_PATH}")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)

v81 = helper.v81


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_prior_freeze():
    manifest = json.loads(PASSIVE_STAGE_A.read_text(encoding="utf-8"))

    checks = {
        "budget_membership": (
            PASSIVE_MEMBERSHIP,
            manifest["hashes"]["budget_membership_sha256"],
        ),
        "all_budget_predictions": (
            PASSIVE_PREDICTIONS,
            manifest["hashes"]["all_budget_predictions_sha256"],
        ),
    }

    verified = {}
    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Prior frozen RhlA artifact mismatch for {name}: "
                f"{actual} != {expected}"
            )
        verified[name] = actual

    return manifest, verified


def support_deficit(mutation_set, model):
    supports = []

    for mutation in mutation_set:
        entry = model["base"]["main"].get(mutation)
        supports.append(0 if entry is None else int(entry["support"]))

    for pair in combinations(mutation_set, 2):
        entry = model["base"]["pair"].get(pair)
        supports.append(0 if entry is None else int(entry["support"]))

    if not supports:
        return 1.0

    values = [
        1.0 / math.sqrt(1.0 + float(support))
        for support in supports
    ]
    return float(np.mean(values))


def router_mode_only(model, visible_ids, visible_labels):
    empty = pd.DataFrame(
        {
            "candidate_id": ["DUMMY"],
            "B3_RAW_PAIR": [0.0],
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": [0.0],
        }
    )
    decision = helper.apply_router(
        empty,
        model,
        visible_ids,
        visible_labels,
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_ADAPTIVE_ROUTER",
    )
    return decision


def score_acquisition_candidates(identity_pool, model, router_mode):
    work = identity_pool.copy()
    work["mutation_set"] = work["mutant"].map(v81.parse_mutations)

    work["structural_exploration"] = work["mutation_set"].map(
        lambda mutation_set: support_deficit(mutation_set, model)
    )

    scoreable = helper.eligible_mask(work, model)
    work["scoreable"] = scoreable.astype(bool)

    scored = helper.score_candidates(work, model)

    if router_mode == "GLOBAL_HIGHER_ORDER":
        scored["exploitation_raw"] = scored[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ].astype(float)
    elif router_mode == "B3_PROTECTED_NO_HIGHER_ORDER":
        scored["exploitation_raw"] = scored["B3_RAW_PAIR"].astype(float)
    elif router_mode == "RANK_PRESERVING_B3_TOP20_B5_RERANK":
        routed = helper.rank_preserving_elite_rerank(
            scored,
            "B3_RAW_PAIR",
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
            "V8_3_ACTIVE_ACQ_SCORE",
            fraction=0.20,
        )
        scored["exploitation_raw"] = scored[
            "V8_3_ACTIVE_ACQ_SCORE"
        ].astype(float)
    else:
        raise RuntimeError(f"Unknown router mode: {router_mode}")

    scored.loc[~scored["scoreable"], "exploitation_raw"] = np.nan

    scored["exploration_rank"] = scored[
        "structural_exploration"
    ].rank(method="average", pct=True)

    exploitation_rank = scored["exploitation_raw"].rank(
        method="average",
        pct=True,
    )
    scored["exploitation_rank"] = exploitation_rank.fillna(0.0)

    scored["acquisition_score"] = (
        0.5 * scored["exploration_rank"]
        + 0.5 * scored["exploitation_rank"]
    )

    return scored


def fit_selected(selected):
    model, visible_ids, labels = helper.fit_model(selected)
    decision = router_mode_only(model, visible_ids, labels)
    return model, visible_ids, labels, decision


def score_fixed_hidden(fixed_hidden, model, decision, prefix):
    scored = helper.score_candidates(fixed_hidden, model)

    if not bool(helper.eligible_mask(fixed_hidden, model).all()):
        raise RuntimeError(
            f"Fixed hidden scoreability invariant failed at {prefix}."
        )

    b3_col = f"{prefix}_B3"
    b5_col = f"{prefix}_B5"
    v83_col = f"{prefix}_V83"

    scored[b3_col] = scored["B3_RAW_PAIR"].astype(float)
    scored[b5_col] = scored[
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
    ].astype(float)

    helper.apply_router(
        scored,
        model,
        [],
        np.array([], dtype=float),
        b3_col,
        b5_col,
        v83_col,
    ) if False else None

    mode = decision["mode"]

    if mode == "GLOBAL_HIGHER_ORDER":
        scored[v83_col] = scored[b5_col].astype(float)
    elif mode == "B3_PROTECTED_NO_HIGHER_ORDER":
        scored[v83_col] = scored[b3_col].astype(float)
    elif mode == "RANK_PRESERVING_B3_TOP20_B5_RERANK":
        helper.rank_preserving_elite_rerank(
            scored,
            b3_col,
            b5_col,
            v83_col,
            fraction=0.20,
        )
    else:
        raise RuntimeError(f"Unknown router mode: {mode}")

    return scored[
        ["candidate_id", b3_col, b5_col, v83_col]
    ].copy()


def main(input_dir, out_dir):
    input_root = Path(input_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    passive_manifest, prior_verified = verify_prior_freeze()

    training_pool = pd.read_csv(input_root / "TRAINING_POOL.csv")
    hidden_ids = pd.read_csv(input_root / "HIDDEN_TRIPLE_IDS.csv")
    membership = pd.read_csv(PASSIVE_MEMBERSHIP)
    passive_predictions = pd.read_csv(PASSIVE_PREDICTIONS)

    fixed_ids = set(passive_predictions["candidate_id"].astype(str))
    fixed_hidden = hidden_ids[
        hidden_ids["candidate_id"].astype(str).isin(fixed_ids)
    ].copy()

    if len(fixed_hidden) != int(passive_manifest["fixed_hidden_rows"]):
        raise RuntimeError(
            "Fixed hidden pool no longer matches prior frozen RhlA test."
        )

    initial_ids = set(
        membership.loc[
            membership["in_05pct"].astype(bool),
            "candidate_id",
        ].astype(str)
    )

    selected = training_pool[
        training_pool["candidate_id"].astype(str).isin(initial_ids)
    ].copy()

    expected_5 = int(passive_manifest["budget_counts"]["5"])
    if len(selected) != expected_5:
        raise RuntimeError(
            f"Initial 5% membership mismatch: {len(selected)} != {expected_5}"
        )

    budget_predictions = fixed_hidden[
        ["candidate_id", "mutant", "mutation_count"]
    ].copy()

    acquisition_log_rows = []
    round_manifests = {}
    router_by_budget = {}
    diagnostics_by_budget = {}

    for budget in TARGET_BUDGETS:
        target_count = int(passive_manifest["budget_counts"][str(budget)])

        if len(selected) > target_count:
            raise RuntimeError(
                f"Selected rows exceed target at {budget}%."
            )

        if len(selected) < target_count:
            model, visible_ids, labels, decision = fit_selected(selected)

            selected_ids = set(selected["candidate_id"].astype(str))
            identity_pool = training_pool[
                ~training_pool["candidate_id"].astype(str).isin(selected_ids)
            ][
                ["candidate_id", "mutant", "mutation_count"]
            ].copy()

            acquisition = score_acquisition_candidates(
                identity_pool,
                model,
                decision["mode"],
            )

            n_add = target_count - len(selected)
            chosen = acquisition.sort_values(
                ["acquisition_score", "candidate_id"],
                ascending=[False, True],
            ).head(n_add).copy()

            selection_path = out / (
                f"ROUND_TO_{budget:02d}PCT_SELECTION.csv"
            )
            chosen[
                [
                    "candidate_id",
                    "mutant",
                    "mutation_count",
                    "scoreable",
                    "structural_exploration",
                    "exploration_rank",
                    "exploitation_raw",
                    "exploitation_rank",
                    "acquisition_score",
                ]
            ].to_csv(selection_path, index=False)

            selection_sha = sha256_file(selection_path)

            # Reveal only after the identity-only acquisition list is frozen.
            reveal = training_pool[
                training_pool["candidate_id"].astype(str).isin(
                    set(chosen["candidate_id"].astype(str))
                )
            ].copy()

            if len(reveal) != n_add:
                raise RuntimeError(
                    f"Round {budget}% reveal size mismatch."
                )

            selected = pd.concat(
                [selected, reveal],
                ignore_index=True,
            ).drop_duplicates("candidate_id", keep="first")

            round_manifests[str(budget)] = {
                "added_rows": int(n_add),
                "selection_sha256_before_label_reveal": selection_sha,
                "selection_file": selection_path.name,
            }

            for rank, row in enumerate(
                chosen.itertuples(),
                start=1,
            ):
                acquisition_log_rows.append(
                    {
                        "target_budget_percent": budget,
                        "selection_rank": rank,
                        "candidate_id": row.candidate_id,
                        "scoreable_before_reveal": bool(row.scoreable),
                        "structural_exploration": float(
                            row.structural_exploration
                        ),
                        "exploration_rank": float(row.exploration_rank),
                        "exploitation_raw": (
                            None
                            if pd.isna(row.exploitation_raw)
                            else float(row.exploitation_raw)
                        ),
                        "exploitation_rank": float(row.exploitation_rank),
                        "acquisition_score": float(row.acquisition_score),
                    }
                )

        if len(selected) != target_count:
            raise RuntimeError(
                f"Budget {budget}% count mismatch after acquisition."
            )

        model, visible_ids, labels, decision = fit_selected(selected)
        diagnostics = v81.model_diagnostics(model)

        prefix = f"A{budget:02d}"
        scored_hidden = score_fixed_hidden(
            fixed_hidden,
            model,
            decision,
            prefix,
        )

        budget_predictions = budget_predictions.merge(
            scored_hidden,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )

        router_by_budget[str(budget)] = decision
        diagnostics_by_budget[str(budget)] = {
            "visible_rows": int(len(selected)),
            "visible_by_mutation_count": {
                str(k): int(v)
                for k, v in selected[
                    "mutation_count"
                ].value_counts().sort_index().items()
            },
            "model_diagnostics": diagnostics,
        }

    acquisition_log = pd.DataFrame(acquisition_log_rows)
    acquisition_log_path = out / "ACQUISITION_LOG.csv"
    acquisition_log.to_csv(acquisition_log_path, index=False)

    predictions_path = out / "ACTIVE_BUDGET_PREDICTIONS.csv"
    budget_predictions.to_csv(predictions_path, index=False)

    freeze_rows = []
    for budget in TARGET_BUDGETS:
        for family in ["B3", "B5", "V83"]:
            arm = f"A{budget:02d}_{family}"
            ranked = budget_predictions.sort_values(
                [arm, "candidate_id"],
                ascending=[False, True],
            ).head(50)

            for rank, row in enumerate(ranked.itertuples(), start=1):
                freeze_rows.append(
                    {
                        "arm": arm,
                        "rank": rank,
                        "candidate_id": row.candidate_id,
                        "predicted_score": float(getattr(row, arm)),
                        "budget10": bool(rank <= 10),
                        "budget50": True,
                    }
                )

    freeze = pd.DataFrame(freeze_rows)
    freeze_path = out / "CANDIDATE_FREEZE.csv"
    freeze.to_csv(freeze_path, index=False)

    prereg_path = HERE.parent / "PREREGISTRATION.json"

    manifest = {
        "version": "NABU_V8_3_RHLA_ACTIVE_ACQUISITION_STAGE_A_V1",
        "status": "ACTIVE_ACQUISITION_AND_HIDDEN_PREDICTIONS_FROZEN",
        "development_only": True,
        "base_architecture_commit": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "prior_rhla_freeze_commit": (
            "7ca494746bd6ec1a5de5e61247aa85cd396a5349"
        ),
        "prior_artifact_hash_verification": prior_verified,
        "fixed_hidden_rows": int(len(fixed_hidden)),
        "initial_5pct_rows": int(expected_5),
        "router_by_budget": router_by_budget,
        "diagnostics_by_budget": diagnostics_by_budget,
        "rounds": round_manifests,
        "hashes": {
            "preregistration_sha256": sha256_file(prereg_path),
            "acquisition_log_sha256": sha256_file(acquisition_log_path),
            "active_budget_predictions_sha256": sha256_file(
                predictions_path
            ),
            "candidate_freeze_sha256": sha256_file(freeze_path),
            "sealed_hidden_truth_sha256_expected": passive_manifest[
                "input_seal"
            ]["hidden_truth_sha256"],
        },
        "hidden_truth_loaded": False,
    }

    (out / "STAGE_A_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))
    print("ACTIVE STAGE A COMPLETE: hidden truth was not loaded.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="rhla_sample_efficiency_sealed_input",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_rhla_active_acquisition_stage_a",
    )
    args = parser.parse_args()
    main(args.input, args.out)

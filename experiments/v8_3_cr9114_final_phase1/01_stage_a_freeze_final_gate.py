import argparse
import hashlib
import importlib.util
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


BUDGETS = [5, 10, 20, 40]
SEED = 161
HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]

HELPER_PATH = (
    ROOT
    / "experiments"
    / "v8_3_rhla_sample_efficiency"
    / "01_stage_a_freeze_budgets.py"
)

spec = importlib.util.spec_from_file_location(
    "nabu_frozen_helper",
    HELPER_PATH,
)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not import frozen helper from {HELPER_PATH}")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)

v81 = helper.v81


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fnv1a32(text):
    h = 2166136261
    for byte in str(text).encode("utf-8"):
        h ^= byte
        h = (h * 16777619) & 0xFFFFFFFF
    return h


def verify_seal(input_root):
    manifest = json.loads(
        (input_root / "SEAL_MANIFEST.json").read_text(
            encoding="utf-8"
        )
    )

    checks = {
        "training_ids": (
            input_root / "TRAINING_IDENTITIES.csv",
            manifest["training_ids_sha256"],
        ),
        "training_labels": (
            input_root / ".TRAINING_LABELS.sealed.csv",
            manifest["training_labels_sha256"],
        ),
        "hidden_ids": (
            input_root / "HIDDEN_4_5_IDS.csv",
            manifest["hidden_ids_sha256"],
        ),
        "hidden_truth": (
            input_root / ".HIDDEN_4_5_TRUTH.sealed.csv",
            manifest["hidden_truth_sha256"],
        ),
    }

    verified = {}
    for name, (path, expected) in checks.items():
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(
                f"Seal hash mismatch for {name}: "
                f"{actual} != {expected}"
            )
        verified[name] = actual

    return manifest, verified


def deterministic_order(training_ids):
    ordered = training_ids.copy()
    ordered["_budget_hash"] = ordered["candidate_id"].map(
        lambda cid: fnv1a32("BUDGET|" + str(cid))
    )
    ordered = ordered.sort_values(
        ["_budget_hash", "candidate_id"],
        ascending=[True, True],
    ).reset_index(drop=True)
    return ordered


def reveal_labels(training_ids, selected_ids, labels_path):
    selected_ids = set(str(x) for x in selected_ids)

    identities = training_ids[
        training_ids["candidate_id"].astype(str).isin(selected_ids)
    ].copy()

    labels = pd.read_csv(labels_path)
    labels = labels[
        labels["candidate_id"].astype(str).isin(selected_ids)
    ][["candidate_id", "DMS_score"]].copy()

    visible = identities.merge(
        labels,
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    if len(visible) != len(selected_ids):
        raise RuntimeError(
            "Training-label reveal mismatch: "
            f"requested={len(selected_ids)} revealed={len(visible)}"
        )

    return visible


def fit_visible(visible):
    model, visible_ids, labels = helper.fit_model(visible)

    dummy = pd.DataFrame({
        "candidate_id": ["DUMMY"],
        "B3_RAW_PAIR": [0.0],
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER": [0.0],
    })

    decision = helper.apply_router(
        dummy,
        model,
        visible_ids,
        labels,
        "B3_RAW_PAIR",
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER",
        "V8_3_DUMMY",
    )

    diagnostics = v81.model_diagnostics(model)

    return {
        "model": model,
        "visible": visible,
        "visible_ids": visible_ids,
        "labels": labels,
        "decision": decision,
        "diagnostics": diagnostics,
    }


def support_deficit(mutation_set, model):
    supports = []

    for mutation in mutation_set:
        entry = model["base"]["main"].get(mutation)
        supports.append(
            0 if entry is None else int(entry["support"])
        )

    for pair in combinations(mutation_set, 2):
        entry = model["base"]["pair"].get(pair)
        supports.append(
            0 if entry is None else int(entry["support"])
        )

    if not supports:
        return 1.0

    return float(
        np.mean(
            [
                1.0 / math.sqrt(1.0 + float(support))
                for support in supports
            ]
        )
    )


def acquisition_scores(identity_pool, fitted):
    work = identity_pool.copy()
    work["mutation_set"] = work["mutant"].map(
        v81.parse_mutations
    )

    work["structural_exploration"] = work[
        "mutation_set"
    ].map(
        lambda mutation_set: support_deficit(
            mutation_set,
            fitted["model"],
        )
    )

    scoreable_mask = helper.eligible_mask(
        work,
        fitted["model"],
    )
    work["scoreable"] = scoreable_mask.astype(bool)
    work["exploitation_raw"] = np.nan
    work["exploitation_rank"] = 0.0

    scoreable = work[work["scoreable"]].copy()

    if not scoreable.empty:
        scored = helper.score_candidates(
            scoreable,
            fitted["model"],
        )

        scored["ACQ_B3"] = scored[
            "B3_RAW_PAIR"
        ].astype(float)
        scored["ACQ_B5"] = scored[
            "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
        ].astype(float)

        helper.apply_router(
            scored,
            fitted["model"],
            fitted["visible_ids"],
            fitted["labels"],
            "ACQ_B3",
            "ACQ_B5",
            "ACQ_V83",
        )

        exploit = scored[
            ["candidate_id", "ACQ_V83"]
        ].rename(
            columns={"ACQ_V83": "exploitation_raw"}
        )

        work = work.drop(
            columns=["exploitation_raw"]
        ).merge(
            exploit,
            on="candidate_id",
            how="left",
            validate="one_to_one",
        )

        scoreable_idx = work["exploitation_raw"].notna()
        work.loc[
            scoreable_idx,
            "exploitation_rank",
        ] = work.loc[
            scoreable_idx,
            "exploitation_raw",
        ].rank(
            method="average",
            pct=True,
        )

    work["exploration_rank"] = work[
        "structural_exploration"
    ].rank(
        method="average",
        pct=True,
    )

    work["acquisition_score"] = (
        0.5 * work["exploration_rank"]
        + 0.5 * work["exploitation_rank"]
    )

    return work


def score_common(common_hidden, fitted, prefix):
    if not bool(
        helper.eligible_mask(
            common_hidden,
            fitted["model"],
        ).all()
    ):
        raise RuntimeError(
            f"Common candidate pool not scoreable for {prefix}."
        )

    scored = helper.score_candidates(
        common_hidden,
        fitted["model"],
    )

    b3_col = f"{prefix}_B3"
    b5_col = f"{prefix}_B5"
    v83_col = f"{prefix}_V83"

    scored[b3_col] = scored[
        "B3_RAW_PAIR"
    ].astype(float)
    scored[b5_col] = scored[
        "B5_CROSSFIT_ADAPTIVE_HIGHER_ORDER"
    ].astype(float)

    helper.apply_router(
        scored,
        fitted["model"],
        fitted["visible_ids"],
        fitted["labels"],
        b3_col,
        b5_col,
        v83_col,
    )

    return scored[
        ["candidate_id", b3_col, b5_col, v83_col]
    ].copy()


def fit_shuffled(visible, common_hidden):
    shuffled_visible = visible.copy()
    shuffled_visible["DMS_score"] = np.random.default_rng(
        SEED
    ).permutation(
        shuffled_visible["DMS_score"].to_numpy(dtype=float)
    )

    fitted = fit_visible(shuffled_visible)

    scored = score_common(
        common_hidden,
        fitted,
        "SHUFFLE40",
    )

    return fitted, scored


def main(input_dir, out_dir):
    input_root = Path(input_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seal, verified_seal = verify_seal(input_root)

    training_ids = pd.read_csv(
        input_root / "TRAINING_IDENTITIES.csv"
    )
    hidden_ids = pd.read_csv(
        input_root / "HIDDEN_4_5_IDS.csv"
    )
    labels_path = (
        input_root / ".TRAINING_LABELS.sealed.csv"
    )

    ordered = deterministic_order(training_ids)

    counts = {
        budget: int(
            seal["budget_target_counts"][str(budget)]
        )
        for budget in BUDGETS
    }

    if max(counts.values()) > len(ordered):
        raise RuntimeError(
            "CR9114 budget exceeds training identity pool."
        )

    # Freeze every passive membership before any training label is revealed.
    membership_rows = []
    passive_ids = {}

    for budget in BUDGETS:
        ids = (
            ordered.head(counts[budget])[
                "candidate_id"
            ].astype(str).tolist()
        )
        passive_ids[budget] = ids

    for row in ordered.itertuples():
        record = {
            "candidate_id": row.candidate_id,
            "budget_hash": int(row._budget_hash),
        }
        for budget in BUDGETS:
            record[f"in_{budget:02d}pct"] = bool(
                str(row.candidate_id)
                in set(passive_ids[budget])
            )
        membership_rows.append(record)

    passive_membership = pd.DataFrame(
        membership_rows
    )
    passive_membership_path = (
        out / "PASSIVE_BUDGET_MEMBERSHIP.csv"
    )
    passive_membership.to_csv(
        passive_membership_path,
        index=False,
    )

    passive = {}

    for budget in BUDGETS:
        visible = reveal_labels(
            training_ids,
            passive_ids[budget],
            labels_path,
        )
        passive[budget] = fit_visible(visible)

    # Active controller begins from the exact same passive 5% identities.
    active_selected_ids = list(passive_ids[5])
    active = {
        5: passive[5],
    }
    active_rounds = {}
    acquisition_log_rows = []

    for target_budget in [10, 20, 40]:
        current_budget = max(
            budget
            for budget in active
            if budget < target_budget
        )
        current_fitted = active[current_budget]

        selected_set = set(
            str(x) for x in active_selected_ids
        )

        identity_pool = training_ids[
            ~training_ids[
                "candidate_id"
            ].astype(str).isin(selected_set)
        ].copy()

        acquisition = acquisition_scores(
            identity_pool,
            current_fitted,
        )

        n_add = (
            counts[target_budget]
            - len(active_selected_ids)
        )
        if n_add <= 0:
            raise RuntimeError(
                f"Invalid active acquisition size at {target_budget}%."
            )

        chosen = acquisition.sort_values(
            ["acquisition_score", "candidate_id"],
            ascending=[False, True],
        ).head(n_add).copy()

        selection_path = (
            out
            / f"ACTIVE_ROUND_TO_{target_budget:02d}PCT_SELECTION.csv"
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
        ].to_csv(
            selection_path,
            index=False,
        )

        selection_sha = sha256_file(
            selection_path
        )

        chosen_ids = chosen[
            "candidate_id"
        ].astype(str).tolist()

        # The selected identity list is now frozen. Only now are labels revealed.
        active_selected_ids.extend(chosen_ids)
        active_selected_ids = list(
            dict.fromkeys(active_selected_ids)
        )

        visible = reveal_labels(
            training_ids,
            active_selected_ids,
            labels_path,
        )

        if len(visible) != counts[target_budget]:
            raise RuntimeError(
                f"Active {target_budget}% row-count mismatch: "
                f"{len(visible)} != {counts[target_budget]}"
            )

        active[target_budget] = fit_visible(
            visible
        )

        active_rounds[str(target_budget)] = {
            "source_budget_percent": int(current_budget),
            "added_rows": int(n_add),
            "selection_file": selection_path.name,
            "selection_sha256_before_label_reveal": (
                selection_sha
            ),
        }

        for rank, row in enumerate(
            chosen.itertuples(),
            start=1,
        ):
            acquisition_log_rows.append({
                "target_budget_percent": int(target_budget),
                "selection_rank": int(rank),
                "candidate_id": row.candidate_id,
                "mutation_count": int(row.mutation_count),
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
                "exploitation_rank": float(
                    row.exploitation_rank
                ),
                "acquisition_score": float(
                    row.acquisition_score
                ),
            })

    acquisition_log = pd.DataFrame(
        acquisition_log_rows
    )
    acquisition_log_path = (
        out / "ACTIVE_ACQUISITION_LOG.csv"
    )
    acquisition_log.to_csv(
        acquisition_log_path,
        index=False,
    )

    # Higher-order qualification is checked before hidden truth is touched.
    p40_diag = passive[40]["diagnostics"]
    a20_diag = active[20]["diagnostics"]

    qualification = {
        "passive_40_triplet_entries": int(
            p40_diag["triplet_crossfit"]["entries"]
        ),
        "passive_40_quartet_entries": int(
            p40_diag["quartet_crossfit"]["entries"]
        ),
        "active_20_triplet_entries": int(
            a20_diag["triplet_crossfit"]["entries"]
        ),
        "active_20_quartet_entries": int(
            a20_diag["quartet_crossfit"]["entries"]
        ),
    }

    qualification["passed"] = bool(
        qualification["passive_40_triplet_entries"] >= 1
        and qualification["passive_40_quartet_entries"] >= 1
        and qualification["active_20_triplet_entries"] >= 1
        and qualification["active_20_quartet_entries"] >= 1
    )

    if not qualification["passed"]:
        abort = {
            "status": "ABORT_BEFORE_REVEAL",
            "reason": (
                "Final CR9114 higher-order qualification failed."
            ),
            "higher_order_qualification": qualification,
        }
        (out / "STAGE_A_QUALIFICATION_ABORT.json").write_text(
            json.dumps(abort, indent=2),
            encoding="utf-8",
        )
        print(json.dumps(abort, indent=2))
        raise SystemExit(2)

    common_ids = set(
        hidden_ids["candidate_id"].astype(str)
    )

    all_fitted = []
    for budget in BUDGETS:
        all_fitted.append(
            ("passive", budget, passive[budget])
        )
        all_fitted.append(
            ("active", budget, active[budget])
        )

    for _, _, fitted in all_fitted:
        mask = helper.eligible_mask(
            hidden_ids,
            fitted["model"],
        )
        common_ids &= set(
            hidden_ids.loc[
                mask,
                "candidate_id",
            ].astype(str)
        )

    common_ids = sorted(common_ids)
    common_hidden = hidden_ids[
        hidden_ids["candidate_id"].astype(str).isin(
            set(common_ids)
        )
    ].copy()

    coverage = float(
        len(common_hidden) / len(hidden_ids)
    )

    if (
        len(common_hidden) < 500
        or coverage < 0.50
    ):
        abort = {
            "status": "ABORT_BEFORE_REVEAL",
            "reason": (
                "Common hidden 4/5-mutant scoreability "
                "qualification failed."
            ),
            "hidden_rows": int(len(hidden_ids)),
            "common_rows": int(len(common_hidden)),
            "coverage_fraction": coverage,
            "required_min_rows": 500,
            "required_min_fraction": 0.50,
            "higher_order_qualification": qualification,
        }
        (out / "STAGE_A_QUALIFICATION_ABORT.json").write_text(
            json.dumps(abort, indent=2),
            encoding="utf-8",
        )
        print(json.dumps(abort, indent=2))
        raise SystemExit(2)

    predictions = common_hidden[
        ["candidate_id", "mutant", "mutation_count"]
    ].copy()

    router_summary = {
        "passive": {},
        "active": {},
    }
    diagnostics_summary = {
        "passive": {},
        "active": {},
    }

    for family, budget, fitted in all_fitted:
        prefix = (
            f"P{budget:02d}"
            if family == "passive"
            else f"A{budget:02d}"
        )

        scored = score_common(
            common_hidden,
            fitted,
            prefix,
        )

        predictions = predictions.merge(
            scored,
            on="candidate_id",
            how="inner",
            validate="one_to_one",
        )

        router_summary[family][str(budget)] = (
            fitted["decision"]
        )
        diagnostics_summary[family][str(budget)] = {
            "visible_rows": int(
                len(fitted["visible"])
            ),
            "visible_by_mutation_count": {
                str(k): int(v)
                for k, v in fitted[
                    "visible"
                ]["mutation_count"].value_counts().sort_index().items()
            },
            "model_diagnostics": fitted[
                "diagnostics"
            ],
        }

    predictions["RANDOM_HASH"] = predictions[
        "candidate_id"
    ].map(
        lambda cid: float(
            fnv1a32("RANDOM|" + str(cid))
            / 4294967296.0
        )
    )

    shuffle_fitted, shuffle_scored = fit_shuffled(
        passive[40]["visible"],
        common_hidden,
    )

    predictions = predictions.merge(
        shuffle_scored[
            ["candidate_id", "SHUFFLE40_V83"]
        ],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )

    predictions_path = (
        out / "ALL_FINAL_GATE_PREDICTIONS.csv"
    )
    predictions.to_csv(
        predictions_path,
        index=False,
    )

    freeze_rows = []
    arms = []

    for budget in BUDGETS:
        for family_prefix in ["P", "A"]:
            for model_name in ["B3", "B5", "V83"]:
                arms.append(
                    f"{family_prefix}{budget:02d}_{model_name}"
                )

    arms.extend(
        ["RANDOM_HASH", "SHUFFLE40_V83"]
    )

    for arm in arms:
        ranked = predictions.sort_values(
            [arm, "candidate_id"],
            ascending=[False, True],
        ).head(50)

        for rank, row in enumerate(
            ranked.itertuples(),
            start=1,
        ):
            freeze_rows.append({
                "arm": arm,
                "rank": int(rank),
                "candidate_id": row.candidate_id,
                "mutation_count": int(row.mutation_count),
                "predicted_score": float(
                    getattr(row, arm)
                ),
                "budget10": bool(rank <= 10),
                "budget50": True,
            })

    freeze = pd.DataFrame(freeze_rows)
    freeze_path = out / "CANDIDATE_FREEZE.csv"
    freeze.to_csv(
        freeze_path,
        index=False,
    )

    prereg_path = (
        HERE.parent / "PREREGISTRATION.json"
    )

    manifest = {
        "version": "NABU_V8_3_CR9114_H1_FINAL_GATE_STAGE_A_V1",
        "status": "FINAL_PHASE1_STAGE_A_FROZEN_BEFORE_HIDDEN_REVEAL",
        "seal_hash_verification": verified_seal,
        "base_architecture_commit": (
            "c8afdcd7698231d95a39ef13a3fe22b6e3f507c5"
        ),
        "active_controller_commit": (
            "749b5e1332e5cf22805624f016e8eeab46fcb346"
        ),
        "budget_counts": {
            str(k): int(v)
            for k, v in counts.items()
        },
        "higher_order_qualification": qualification,
        "hidden_4_5_rows_total": int(
            len(hidden_ids)
        ),
        "fixed_common_hidden_rows": int(
            len(common_hidden)
        ),
        "fixed_common_hidden_coverage_fraction": coverage,
        "router_by_budget": router_summary,
        "diagnostics_by_budget": diagnostics_summary,
        "active_rounds": active_rounds,
        "shuffle_40": {
            "router": shuffle_fitted[
                "decision"
            ],
            "model_diagnostics": shuffle_fitted[
                "diagnostics"
            ],
        },
        "hashes": {
            "preregistration_sha256": sha256_file(
                prereg_path
            ),
            "passive_membership_sha256": sha256_file(
                passive_membership_path
            ),
            "active_acquisition_log_sha256": sha256_file(
                acquisition_log_path
            ),
            "all_predictions_sha256": sha256_file(
                predictions_path
            ),
            "candidate_freeze_sha256": sha256_file(
                freeze_path
            ),
            "sealed_training_labels_sha256_expected": seal[
                "training_labels_sha256"
            ],
            "sealed_hidden_truth_sha256_expected": seal[
                "hidden_truth_sha256"
            ],
        },
        "hidden_truth_loaded": False,
    }

    (out / "STAGE_A_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2))
    print(
        "FINAL PHASE-1 STAGE A COMPLETE: "
        "passive, active, higher-order and controls are frozen; "
        "hidden 4/5 truth was not loaded."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="cr9114_h1_final_gate_sealed_input",
    )
    parser.add_argument(
        "--out",
        default="nabu_v8_3_cr9114_final_gate_stage_a",
    )
    args = parser.parse_args()
    main(args.input, args.out)

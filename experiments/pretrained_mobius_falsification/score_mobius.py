from __future__ import annotations

import argparse
import hashlib
import json
import platform
from itertools import combinations
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from huggingface_hub import model_info, snapshot_download
from transformers import AutoModelForMaskedLM, AutoTokenizer

SEED = 161
MODEL_ID = "facebook/esm2_t6_8M_UR50D"
MODEL_REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"
EXPECTED_SAFETENSORS_SHA256 = (
    "24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189"
)
EXPECTED_SAFETENSORS_SIZE = 31384292


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_token(token: str):
    token = str(token)
    source = token[0]
    target = token[-1]
    position = int(token[1:-1])
    return source, position, target


def mutate(reference: str, tokens: tuple[str, ...]) -> str:
    chars = list(reference)
    seen = set()
    for token in tokens:
        source, position, target = parse_token(token)
        idx = position - 1
        if idx in seen:
            raise RuntimeError(f"Duplicate mutation position: {tokens}")
        seen.add(idx)
        if chars[idx] != source:
            raise RuntimeError(
                f"Reference mismatch for {token}: {chars[idx]} != {source}"
            )
        chars[idx] = target
    return "".join(chars)


def subset_label(tokens: tuple[str, ...]) -> str:
    return ";".join(tokens)


def request_key(dataset: str, subset: tuple[str, ...], token: str) -> str:
    return f"{dataset}|{subset_label(subset)}|MASK={token}"


def required_subsets(tokens: tuple[str, ...]):
    out = []
    for token in tokens:
        out.append((token,))
    for pair in combinations(tokens, 2):
        out.append(tuple(pair))
    out.append(tuple(tokens))
    # Preserve order while removing duplicates.
    seen = set()
    unique = []
    for subset in out:
        if subset not in seen:
            seen.add(subset)
            unique.append(subset)
    return unique


def build_plan(selected: pd.DataFrame):
    references = {}
    requests = {}
    candidates = []

    for row in selected.itertuples(index=False):
        dataset = str(row.dataset)
        cid = str(row.candidate_id)
        reference = str(row.reference_sequence)
        tokens = tuple(json.loads(row.mutation_tokens))

        if dataset in references and references[dataset] != reference:
            raise RuntimeError(f"{dataset}: inconsistent reference sequence.")
        references[dataset] = reference

        full_sequence = mutate(reference, tokens)
        if full_sequence != str(row.candidate_sequence):
            raise RuntimeError(
                f"{dataset}/{cid}: candidate reconstruction mismatch."
            )

        subsets = required_subsets(tokens)
        subset_keys = []

        for subset in subsets:
            seq = mutate(reference, subset)
            skey = subset_label(subset)
            subset_keys.append(skey)

            for token in subset:
                source, position, target = parse_token(token)
                key = request_key(dataset, subset, token)
                request = {
                    "request_key": key,
                    "dataset": dataset,
                    "subset_tokens": skey,
                    "subset_size": len(subset),
                    "mutation_token": token,
                    "sequence": seq,
                    "position": position,
                    "wildtype_aa": source,
                    "mutant_aa": target,
                }
                previous = requests.get(key)
                if previous is not None and previous != request:
                    raise RuntimeError(f"Request collision: {key}")
                requests[key] = request

        candidates.append(
            {
                "dataset": dataset,
                "order": int(row.order),
                "candidate_id": cid,
                "selection_rank": int(row.selection_rank),
                "selection_hash": str(row.selection_hash),
                "tokens": tokens,
                "subset_keys": subset_keys,
            }
        )

    return references, list(requests.values()), candidates


def resolve_model(model_cache: Path):
    model_cache.mkdir(parents=True, exist_ok=True)

    resolved = model_info(MODEL_ID, revision=MODEL_REVISION).sha
    if resolved != MODEL_REVISION:
        raise RuntimeError(
            f"Resolved model revision mismatch: {resolved}"
        )

    snapshot = Path(
        snapshot_download(
            repo_id=MODEL_ID,
            revision=MODEL_REVISION,
            cache_dir=str(model_cache),
            allow_patterns=[
                "config.json",
                "model.safetensors",
                "special_tokens_map.json",
                "tokenizer_config.json",
                "vocab.txt",
            ],
        )
    )

    weights = snapshot / "model.safetensors"
    if not weights.exists():
        raise RuntimeError("Pinned model.safetensors missing.")
    if weights.stat().st_size != EXPECTED_SAFETENSORS_SIZE:
        raise RuntimeError(
            f"Model size mismatch: {weights.stat().st_size}"
        )
    observed = sha256_file(weights)
    if observed != EXPECTED_SAFETENSORS_SHA256:
        raise RuntimeError(
            f"Model weights SHA256 mismatch: {observed}"
        )

    files = {}
    for path in sorted(snapshot.iterdir()):
        if path.is_file():
            files[path.name] = {
                "size": int(path.stat().st_size),
                "sha256": sha256_file(path),
            }

    return snapshot, {
        "model_id": MODEL_ID,
        "requested_revision": MODEL_REVISION,
        "resolved_revision": resolved,
        "expected_model_safetensors_sha256": EXPECTED_SAFETENSORS_SHA256,
        "expected_model_safetensors_size": EXPECTED_SAFETENSORS_SIZE,
        "files": files,
    }


def score_requests(requests, tokenizer, model, batch_size):
    out = {}
    mask_id = tokenizer.mask_token_id
    if mask_id is None:
        raise RuntimeError("Tokenizer has no mask token.")

    aa_ids = {
        aa: tokenizer.convert_tokens_to_ids(aa)
        for aa in "ACDEFGHIKLMNPQRSTVWY"
    }
    if any(v is None or int(v) < 0 for v in aa_ids.values()):
        raise RuntimeError("Amino-acid token mapping failed.")

    for start in range(0, len(requests), batch_size):
        batch = requests[start : start + batch_size]
        sequences = [x["sequence"] for x in batch]

        enc = tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            add_special_tokens=True,
        )
        input_ids = enc["input_ids"].clone()

        for row_idx, req in enumerate(batch):
            pos = int(req["position"])
            seq = req["sequence"]
            token_pos = pos  # ESM tokenizer includes BOS at token index 0.
            expected_aa = seq[pos - 1]
            expected_id = int(aa_ids[expected_aa])
            observed_id = int(input_ids[row_idx, token_pos].item())
            if observed_id != expected_id:
                raise RuntimeError(
                    f"Tokenizer alignment failed for {req['request_key']}: "
                    f"{observed_id} != {expected_id}"
                )
            input_ids[row_idx, token_pos] = mask_id

        enc["input_ids"] = input_ids
        with torch.no_grad():
            logits = model(**enc).logits

        for row_idx, req in enumerate(batch):
            token_pos = int(req["position"])
            lp = torch.log_softmax(
                logits[row_idx, token_pos].float(),
                dim=-1,
            )
            wt_id = int(aa_ids[req["wildtype_aa"]])
            mut_id = int(aa_ids[req["mutant_aa"]])
            wt_logp = float(lp[wt_id].item())
            mut_logp = float(lp[mut_id].item())
            preference = float(mut_logp - wt_logp)

            out[req["request_key"]] = {
                "wildtype_logp": wt_logp,
                "mutant_logp": mut_logp,
                "preference": preference,
            }

    return out


def f_value(dataset, subset, scores):
    subset = tuple(subset)
    return float(
        sum(
            scores[request_key(dataset, subset, token)]["preference"]
            for token in subset
        )
    )


def candidate_metrics(candidate, scores):
    dataset = candidate["dataset"]
    tokens = tuple(candidate["tokens"])

    singles = {
        token: f_value(dataset, (token,), scores)
        for token in tokens
    }

    pair_f = {}
    pair_e2 = {}
    for pair in combinations(tokens, 2):
        pair = tuple(pair)
        pair_f[pair] = f_value(dataset, pair, scores)
        pair_e2[pair] = float(
            pair_f[pair] - singles[pair[0]] - singles[pair[1]]
        )

    full_f = f_value(dataset, tokens, scores)
    model_o2 = float(
        sum(singles.values()) + sum(pair_e2.values())
    )
    model_ho = float(full_f - model_o2)

    full_preferences = np.asarray(
        [
            scores[request_key(dataset, tokens, token)]["preference"]
            for token in tokens
        ],
        dtype=float,
    )
    single_preferences = np.asarray(
        [singles[token] for token in tokens],
        dtype=float,
    )
    context_shift = full_preferences - single_preferences

    return {
        "model_f_full": full_f,
        "model_o2": model_o2,
        "model_ho": model_ho,
        "context_mean_control": float(np.mean(context_shift)),
        "context_sum_control": float(np.sum(context_shift)),
        "single_mean_control": float(np.mean(single_preferences)),
        "single_sum_control": float(np.sum(single_preferences)),
    }


def aggregate(selected, requests, candidates, scores):
    request_by_key = {x["request_key"]: x for x in requests}

    masked_rows = []
    for key in sorted(scores):
        req = request_by_key[key]
        val = scores[key]
        masked_rows.append(
            {
                "request_key": key,
                "dataset": req["dataset"],
                "subset_tokens": req["subset_tokens"],
                "subset_size": req["subset_size"],
                "mutation_token": req["mutation_token"],
                "position": req["position"],
                "wildtype_aa": req["wildtype_aa"],
                "mutant_aa": req["mutant_aa"],
                "wildtype_logp": val["wildtype_logp"],
                "mutant_logp": val["mutant_logp"],
                "preference": val["preference"],
            }
        )

    subset_rows = []
    candidate_rows = []

    for candidate in candidates:
        dataset = candidate["dataset"]
        cid = candidate["candidate_id"]
        tokens = tuple(candidate["tokens"])

        for subset in required_subsets(tokens):
            subset_rows.append(
                {
                    "dataset": dataset,
                    "order": candidate["order"],
                    "candidate_id": cid,
                    "selection_rank": candidate["selection_rank"],
                    "subset_tokens": subset_label(subset),
                    "subset_size": len(subset),
                    "model_f": f_value(dataset, subset, scores),
                }
            )

        metrics = candidate_metrics(candidate, scores)
        candidate_rows.append(
            {
                "dataset": dataset,
                "order": candidate["order"],
                "candidate_id": cid,
                "selection_rank": candidate["selection_rank"],
                "selection_hash": candidate["selection_hash"],
                **metrics,
            }
        )

    return (
        pd.DataFrame(masked_rows),
        pd.DataFrame(subset_rows),
        pd.DataFrame(candidate_rows),
    )


def determinism_check(
    selected,
    requests,
    candidates,
    tokenizer,
    model,
    batch_size,
    first_scores,
    first_candidate_df,
):
    chosen = set()
    for _, cell in selected.groupby(["dataset", "order"], sort=True):
        for row in (
            cell.sort_values("selection_rank")
            .head(8)
            .itertuples(index=False)
        ):
            chosen.add((str(row.dataset), str(row.candidate_id)))

    candidate_lookup = {
        (x["dataset"], x["candidate_id"]): x
        for x in candidates
    }
    chosen_candidates = [candidate_lookup[key] for key in sorted(chosen)]

    needed_keys = set()
    for candidate in chosen_candidates:
        dataset = candidate["dataset"]
        tokens = tuple(candidate["tokens"])
        for subset in required_subsets(tokens):
            for token in subset:
                needed_keys.add(request_key(dataset, subset, token))

    # Replay the exact full request stream with identical ordering and
    # batch boundaries. Comparing a reduced subset would change padding and
    # matrix shapes, which can introduce tiny deterministic floating-point
    # differences even when the model/runtime is otherwise reproducible.
    replay_requests = list(requests)
    replay_scores = score_requests(
        replay_requests,
        tokenizer,
        model,
        batch_size,
    )

    masked_diffs = []
    for key in sorted(needed_keys):
        for field in ("wildtype_logp", "mutant_logp", "preference"):
            masked_diffs.append(
                abs(
                    float(first_scores[key][field])
                    - float(replay_scores[key][field])
                )
            )

    first_by_candidate = {
        (str(r.dataset), str(r.candidate_id)): r
        for r in first_candidate_df.itertuples(index=False)
    }

    derived_diffs = []
    for candidate in chosen_candidates:
        first = first_by_candidate[
            (candidate["dataset"], candidate["candidate_id"])
        ]
        replay = candidate_metrics(candidate, replay_scores)
        for field in (
            "model_f_full",
            "model_o2",
            "model_ho",
            "context_mean_control",
            "context_sum_control",
            "single_mean_control",
            "single_sum_control",
        ):
            derived_diffs.append(
                abs(float(getattr(first, field)) - float(replay[field]))
            )

    masked_max = 0.0 if not masked_diffs else float(max(masked_diffs))
    derived_max = 0.0 if not derived_diffs else float(max(derived_diffs))
    overall = max(masked_max, derived_max)

    return {
        "candidate_count": int(len(chosen_candidates)),
        "request_count": int(len(replay_requests)),
        "masked_score_max_abs_difference": masked_max,
        "derived_score_max_abs_difference": derived_max,
        "overall_max_abs_difference": overall,
        "exact": bool(overall == 0.0),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--visible", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--model-cache", required=True)
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()

    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.set_num_threads(2)
    torch.use_deterministic_algorithms(True)

    visible = Path(args.visible)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(visible / "SELECTED_CANDIDATES.csv")
    expected = 2 * 2 * 192
    if len(selected) != expected:
        raise RuntimeError(
            f"Visible candidate count mismatch: {len(selected)}"
        )

    if int(selected["selection_rank"].min()) != 192:
        raise RuntimeError("Fresh-holdout selection rank must start at 192.")
    if int(selected["selection_rank"].max()) != 383:
        raise RuntimeError("Fresh-holdout selection rank must end at 383.")

    references, requests, candidates = build_plan(selected)

    snapshot, model_manifest = resolve_model(Path(args.model_cache))
    tokenizer = AutoTokenizer.from_pretrained(
        snapshot,
        local_files_only=True,
    )
    model = AutoModelForMaskedLM.from_pretrained(
        snapshot,
        local_files_only=True,
        use_safetensors=True,
    )
    model.eval()

    max_position_embeddings = int(
        getattr(model.config, "max_position_embeddings", 0)
    )
    max_sequence_length = int(
        max(selected["candidate_sequence"].astype(str).str.len())
    )
    if (
        max_position_embeddings
        and max_sequence_length + 2 > max_position_embeddings
    ):
        raise RuntimeError("Sequence exceeds ESM context length.")

    scores = score_requests(
        requests,
        tokenizer,
        model,
        args.batch_size,
    )

    masked_df, subset_df, candidate_df = aggregate(
        selected,
        requests,
        candidates,
        scores,
    )

    det = determinism_check(
        selected,
        requests,
        candidates,
        tokenizer,
        model,
        args.batch_size,
        scores,
        candidate_df,
    )
    if not det["exact"]:
        raise RuntimeError(f"Determinism contract failed: {det}")

    masked_df.to_csv(out / "MASKED_LOG_ODDS.csv", index=False)
    subset_df.to_csv(out / "SUBSET_MODEL_F.csv", index=False)
    candidate_df.to_csv(
        out / "CANDIDATE_MODEL_SCORES.csv",
        index=False,
    )

    model_manifest.update(
        {
            "architecture": "EsmForMaskedLM",
            "eval_mode": True,
            "fine_tuned": False,
            "gradients_enabled": False,
            "max_position_embeddings": max_position_embeddings,
            "maximum_sequence_length_observed": max_sequence_length,
        }
    )
    (out / "MODEL_MANIFEST.json").write_text(
        json.dumps(model_manifest, indent=2),
        encoding="utf-8",
    )

    runtime = {
        "version": "NABU_PRETRAINED_MOBIUS_SCORING_V1",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "seed": SEED,
        "torch_num_threads": int(torch.get_num_threads()),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "torch": version("torch"),
        "transformers": version("transformers"),
        "huggingface_hub": version("huggingface_hub"),
        "safetensors": version("safetensors"),
        "batch_size": int(args.batch_size),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUNTIME_MANIFEST.json").write_text(
        json.dumps(runtime, indent=2),
        encoding="utf-8",
    )

    diagnostics = {
        "selected_candidates": int(len(selected)),
        "unique_masked_requests": int(len(requests)),
        "masked_score_rows": int(len(masked_df)),
        "subset_model_f_rows": int(len(subset_df)),
        "candidate_score_rows": int(len(candidate_df)),
        "request_counts_by_dataset": {
            str(k): int(v)
            for k, v in masked_df.groupby("dataset").size().items()
        },
        "candidate_counts_by_cell": {
            f"{dataset}|{int(order)}": int(count)
            for (dataset, order), count in (
                candidate_df.groupby(["dataset", "order"]).size().items()
            )
        },
        "determinism": det,
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "SCORING_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    hashes = {
        path.name: sha256_file(path)
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "SCORE_HASHES.json"
    }
    (out / "SCORE_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

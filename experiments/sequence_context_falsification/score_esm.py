from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from collections import defaultdict
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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_token(token: str):
    token = str(token)
    if len(token) < 3:
        raise RuntimeError(f"Invalid mutation token: {token}")
    source = token[0]
    target = token[-1]
    position = int(token[1:-1])
    return source, position, target


def single_mutant_sequence(reference: str, token: str) -> str:
    source, position, target = parse_token(token)
    chars = list(reference)
    idx = position - 1
    if chars[idx] != source:
        raise RuntimeError(
            f"Reference/token mismatch for {token}: {chars[idx]}"
        )
    chars[idx] = target
    return "".join(chars)


def request_key(dataset, candidate_id, token, context_type):
    return f"{dataset}|{candidate_id}|{token}|{context_type}"


def build_requests(selected: pd.DataFrame):
    candidate_requests = []
    single_requests_by_key = {}

    for row in selected.itertuples(index=False):
        mutations = json.loads(row.mutation_tokens)
        reference = str(row.reference_sequence)
        candidate = str(row.candidate_sequence)

        for token in mutations:
            source, position, target = parse_token(token)
            if reference[position - 1] != source:
                raise RuntimeError(
                    f"{row.dataset}/{row.candidate_id}: "
                    f"reference mismatch for {token}"
                )
            if candidate[position - 1] != target:
                raise RuntimeError(
                    f"{row.dataset}/{row.candidate_id}: "
                    f"candidate mismatch for {token}"
                )

            candidate_requests.append(
                {
                    "request_key": request_key(
                        row.dataset,
                        row.candidate_id,
                        token,
                        "candidate",
                    ),
                    "dataset": row.dataset,
                    "candidate_id": str(row.candidate_id),
                    "order": int(row.order),
                    "mutation_token": token,
                    "context_type": "candidate",
                    "sequence": candidate,
                    "position": position,
                    "wildtype_aa": source,
                    "mutant_aa": target,
                }
            )

            single_key = f"{row.dataset}|{token}"
            if single_key not in single_requests_by_key:
                single_requests_by_key[single_key] = {
                    "request_key": single_key,
                    "dataset": row.dataset,
                    "candidate_id": "__SINGLE__",
                    "order": 1,
                    "mutation_token": token,
                    "context_type": "single",
                    "sequence": single_mutant_sequence(reference, token),
                    "position": position,
                    "wildtype_aa": source,
                    "mutant_aa": target,
                }

    return candidate_requests, list(single_requests_by_key.values())


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
        raise RuntimeError("Pinned safetensors weights missing.")
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
        "files": files,
    }


def score_requests(requests, tokenizer, model, batch_size):
    out = {}
    mask_id = tokenizer.mask_token_id
    if mask_id is None:
        raise RuntimeError("Tokenizer mask token is unavailable.")

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
            if pos < 1 or pos > len(seq):
                raise RuntimeError("Mask position out of range.")

            token_pos = pos
            expected_aa = seq[pos - 1]
            expected_id = int(aa_ids[expected_aa])
            observed_id = int(input_ids[row_idx, token_pos].item())
            if observed_id != expected_id:
                raise RuntimeError(
                    f"Tokenizer residue alignment failed at {req['request_key']}: "
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


def aggregate(
    selected,
    candidate_requests,
    single_requests,
    candidate_scores,
    single_scores,
):
    mutation_rows = []
    candidate_acc = defaultdict(list)
    single_acc = defaultdict(list)

    single_by_dataset_token = {
        (req["dataset"], req["mutation_token"]): single_scores[
            req["request_key"]
        ]
        for req in single_requests
    }

    candidate_meta = {
        (str(r.dataset), str(r.candidate_id)): r
        for r in selected.itertuples(index=False)
    }

    for req in candidate_requests:
        cscore = candidate_scores[req["request_key"]]
        sscore = single_by_dataset_token[
            (req["dataset"], req["mutation_token"])
        ]
        shift = float(cscore["preference"] - sscore["preference"])

        key = (req["dataset"], req["candidate_id"])
        candidate_acc[key].append(shift)
        single_acc[key].append(float(sscore["preference"]))

        mutation_rows.append(
            {
                "dataset": req["dataset"],
                "order": req["order"],
                "candidate_id": req["candidate_id"],
                "mutation_token": req["mutation_token"],
                "position": req["position"],
                "wildtype_aa": req["wildtype_aa"],
                "mutant_aa": req["mutant_aa"],
                "candidate_wildtype_logp": cscore["wildtype_logp"],
                "candidate_mutant_logp": cscore["mutant_logp"],
                "candidate_preference": cscore["preference"],
                "single_wildtype_logp": sscore["wildtype_logp"],
                "single_mutant_logp": sscore["mutant_logp"],
                "single_preference": sscore["preference"],
                "context_shift": shift,
            }
        )

    candidate_rows = []
    for key, shifts in candidate_acc.items():
        dataset, cid = key
        meta = candidate_meta[key]
        shifts = np.asarray(shifts, dtype=float)
        singles = np.asarray(single_acc[key], dtype=float)
        candidate_rows.append(
            {
                "dataset": dataset,
                "order": int(meta.order),
                "candidate_id": cid,
                "selection_rank": int(meta.selection_rank),
                "selection_hash": str(meta.selection_hash),
                "esm_context_mean": float(np.mean(shifts)),
                "esm_context_sum": float(np.sum(shifts)),
                "esm_context_std": float(np.std(shifts)),
                "esm_single_mean": float(np.mean(singles)),
                "mutation_count": int(len(shifts)),
            }
        )

    return pd.DataFrame(mutation_rows), pd.DataFrame(candidate_rows)


def determinism_check(
    selected,
    candidate_requests,
    tokenizer,
    model,
    batch_size,
    first_scores,
):
    chosen_ids = set()
    for (_, _), cell in selected.groupby(["dataset", "order"], sort=True):
        for row in cell.sort_values("selection_rank").head(16).itertuples():
            chosen_ids.add((str(row.dataset), str(row.candidate_id)))

    replay_requests = [
        req
        for req in candidate_requests
        if (req["dataset"], req["candidate_id"]) in chosen_ids
    ]
    replay = score_requests(
        replay_requests,
        tokenizer,
        model,
        batch_size,
    )

    diffs = []
    for req in replay_requests:
        key = req["request_key"]
        a = first_scores[key]
        b = replay[key]
        for field in (
            "wildtype_logp",
            "mutant_logp",
            "preference",
        ):
            diffs.append(abs(float(a[field]) - float(b[field])))

    max_abs = 0.0 if not diffs else float(max(diffs))
    return {
        "candidate_count": int(len(chosen_ids)),
        "mutation_request_count": int(len(replay_requests)),
        "max_abs_probability_score_difference": max_abs,
        "exact": bool(max_abs == 0.0),
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

    candidate_requests, single_requests = build_requests(selected)

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
    if max_position_embeddings and max_sequence_length + 2 > max_position_embeddings:
        raise RuntimeError(
            "Sequence exceeds pinned ESM context length."
        )

    single_scores = score_requests(
        single_requests,
        tokenizer,
        model,
        args.batch_size,
    )
    candidate_scores = score_requests(
        candidate_requests,
        tokenizer,
        model,
        args.batch_size,
    )

    mutation_df, candidate_df = aggregate(
        selected,
        candidate_requests,
        single_requests,
        candidate_scores,
        single_scores,
    )

    det = determinism_check(
        selected,
        candidate_requests,
        tokenizer,
        model,
        args.batch_size,
        candidate_scores,
    )
    if not det["exact"]:
        raise RuntimeError(f"Determinism contract failed: {det}")

    mutation_df.to_csv(
        out / "ESM_MUTATION_SCORES.csv",
        index=False,
    )
    candidate_df.to_csv(
        out / "ESM_CANDIDATE_SCORES.csv",
        index=False,
    )

    model_manifest.update(
        {
            "architecture": "EsmForMaskedLM",
            "eval_mode": True,
            "fine_tuned": False,
            "gradient_enabled": False,
            "max_position_embeddings": max_position_embeddings,
            "max_sequence_length_observed": max_sequence_length,
        }
    )
    (out / "MODEL_MANIFEST.json").write_text(
        json.dumps(model_manifest, indent=2),
        encoding="utf-8",
    )

    request_counts = {
        "candidate_masked_requests": int(len(candidate_requests)),
        "unique_single_masked_requests": int(len(single_requests)),
        "total_masked_forward_examples": int(
            len(candidate_requests) + len(single_requests)
        ),
        "candidate_rows": int(len(candidate_df)),
        "mutation_rows": int(len(mutation_df)),
        "sequence_length_by_dataset": {
            dataset: {
                "min": int(
                    cell["candidate_sequence"].astype(str).str.len().min()
                ),
                "max": int(
                    cell["candidate_sequence"].astype(str).str.len().max()
                ),
            }
            for dataset, cell in selected.groupby("dataset")
        },
        "unique_mutations_by_dataset": {
            dataset: int(
                len(
                    {
                        token
                        for raw in cell["mutation_tokens"]
                        for token in json.loads(raw)
                    }
                )
            )
            for dataset, cell in selected.groupby("dataset")
        },
        "determinism": det,
    }
    (out / "SCORING_DIAGNOSTICS.json").write_text(
        json.dumps(request_counts, indent=2),
        encoding="utf-8",
    )

    runtime = {
        "python": sys.version,
        "platform": platform.platform(),
        "seed": SEED,
        "torch_num_threads": int(torch.get_num_threads()),
        "packages": {
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "torch": version("torch"),
            "transformers": version("transformers"),
            "huggingface_hub": version("huggingface-hub"),
            "safetensors": version("safetensors"),
        },
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out / "RUNTIME_MANIFEST.json").write_text(
        json.dumps(runtime, indent=2),
        encoding="utf-8",
    )

    hashes = {
        p.name: sha256_file(p)
        for p in sorted(out.iterdir())
        if p.is_file()
    }
    (out / "SCORE_OUTPUT_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(request_counts, indent=2))


if __name__ == "__main__":
    main()

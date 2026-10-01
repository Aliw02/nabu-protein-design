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
    return token[0], int(token[1:-1]), token[-1]


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
    return ";".join(tokens) if tokens else "__WT__"


def required_subsets(tokens: tuple[str, ...]):
    out = [tuple()]
    out.extend((t,) for t in tokens)
    out.extend(tuple(pair) for pair in combinations(tokens, 2))
    out.append(tuple(tokens))
    seen = set()
    unique = []
    for subset in out:
        if subset not in seen:
            seen.add(subset)
            unique.append(subset)
    return unique


def resolve_model(model_cache: Path):
    model_cache.mkdir(parents=True, exist_ok=True)
    resolved = model_info(MODEL_ID, revision=MODEL_REVISION).sha
    if resolved != MODEL_REVISION:
        raise RuntimeError(f"Resolved revision mismatch: {resolved}")

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


def build_items(selected: pd.DataFrame):
    items = []
    candidate_specs = []

    for row in selected.itertuples(index=False):
        dataset = str(row.dataset)
        cid = str(row.candidate_id)
        reference = str(row.reference_sequence)
        tokens = tuple(json.loads(row.mutation_tokens))
        positions = tuple(parse_token(t)[1] for t in tokens)

        reconstructed = mutate(reference, tokens)
        if reconstructed != str(row.candidate_sequence):
            raise RuntimeError(
                f"{dataset}/{cid}: candidate reconstruction mismatch."
            )

        subset_indices = []
        for subset in required_subsets(tokens):
            idx = len(items)
            subset_indices.append(idx)
            items.append(
                {
                    "dataset": dataset,
                    "order": int(row.order),
                    "candidate_id": cid,
                    "selection_rank": int(row.selection_rank),
                    "selection_hash": str(row.selection_hash),
                    "subset_tokens": tuple(subset),
                    "subset_label": subset_label(tuple(subset)),
                    "subset_size": int(len(subset)),
                    "sequence": mutate(reference, tuple(subset)),
                    "candidate_positions": positions,
                }
            )

        candidate_specs.append(
            {
                "dataset": dataset,
                "order": int(row.order),
                "candidate_id": cid,
                "selection_rank": int(row.selection_rank),
                "selection_hash": str(row.selection_hash),
                "tokens": tokens,
                "subset_indices": subset_indices,
            }
        )

    return items, candidate_specs


def extract(items, tokenizer, model, batch_size):
    vectors = []
    hidden_size = int(model.config.hidden_size)

    for start in range(0, len(items), batch_size):
        batch = items[start : start + batch_size]
        sequences = [x["sequence"] for x in batch]

        enc = tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            add_special_tokens=True,
        )
        with torch.no_grad():
            out = model(
                **enc,
                output_hidden_states=True,
                return_dict=True,
            )
        hidden = out.hidden_states[-1]

        for row_idx, item in enumerate(batch):
            token_positions = [
                int(position)
                for position in item["candidate_positions"]
            ]
            if not token_positions:
                raise RuntimeError("Candidate position set is empty.")

            vec = hidden[
                row_idx,
                token_positions,
                :,
            ].float().mean(dim=0)
            arr = vec.cpu().numpy().astype(np.float32, copy=False)
            if arr.shape != (hidden_size,):
                raise RuntimeError(
                    f"Unexpected vector shape: {arr.shape}"
                )
            if not np.isfinite(arr).all():
                raise RuntimeError("Non-finite representation.")
            vectors.append(arr.copy())

    matrix = np.stack(vectors, axis=0)
    return matrix


def candidate_r_ho(candidate, items, matrix):
    tokens = tuple(candidate["tokens"])
    index_map = {}
    for idx in candidate["subset_indices"]:
        item = items[idx]
        index_map[tuple(item["subset_tokens"])] = matrix[idx]

    empty = index_map[tuple()]
    e1 = {
        token: index_map[(token,)] - empty
        for token in tokens
    }
    e2 = {}
    for pair in combinations(tokens, 2):
        pair = tuple(pair)
        e2[pair] = (
            index_map[pair]
            - empty
            - e1[pair[0]]
            - e1[pair[1]]
        )

    e_o2 = empty.copy()
    for token in tokens:
        e_o2 = e_o2 + e1[token]
    for pair in combinations(tokens, 2):
        e_o2 = e_o2 + e2[tuple(pair)]

    full = index_map[tokens]
    r_ho = full - e_o2
    return r_ho.astype(np.float32, copy=False)


def build_candidate_matrix(candidate_specs, items, matrix):
    rows = []
    vectors = []

    for candidate in candidate_specs:
        vec = candidate_r_ho(candidate, items, matrix)
        rows.append(
            {
                "dataset": candidate["dataset"],
                "order": candidate["order"],
                "candidate_id": candidate["candidate_id"],
                "selection_rank": candidate["selection_rank"],
                "selection_hash": candidate["selection_hash"],
                "r_ho_l2": float(np.linalg.norm(vec)),
                "r_ho_mean": float(np.mean(vec)),
                "r_ho_std": float(np.std(vec)),
            }
        )
        vectors.append(vec)

    return pd.DataFrame(rows), np.stack(vectors, axis=0)


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
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(visible / "SELECTED_CANDIDATES.csv")
    if len(selected) != 768:
        raise RuntimeError(f"Selected row count mismatch: {len(selected)}")
    if int(selected["selection_rank"].min()) != 384:
        raise RuntimeError("Fresh holdout must start at rank 384.")
    if int(selected["selection_rank"].max()) != 575:
        raise RuntimeError("Fresh holdout must end at rank 575.")

    items, candidate_specs = build_items(selected)

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

    hidden_size = int(model.config.hidden_size)
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

    matrix_first = extract(
        items,
        tokenizer,
        model,
        args.batch_size,
    )
    matrix_replay = extract(
        items,
        tokenizer,
        model,
        args.batch_size,
    )

    max_abs = float(
        np.max(np.abs(matrix_first - matrix_replay))
    )
    if max_abs != 0.0:
        raise RuntimeError(
            f"Representation determinism failed: max_abs={max_abs}"
        )

    subset_meta = pd.DataFrame(
        [
            {
                "row_index": i,
                "dataset": item["dataset"],
                "order": item["order"],
                "candidate_id": item["candidate_id"],
                "selection_rank": item["selection_rank"],
                "selection_hash": item["selection_hash"],
                "subset_tokens": item["subset_label"],
                "subset_size": item["subset_size"],
                "candidate_positions": json.dumps(
                    list(item["candidate_positions"])
                ),
            }
            for i, item in enumerate(items)
        ]
    )

    candidate_meta, r_ho = build_candidate_matrix(
        candidate_specs,
        items,
        matrix_first,
    )

    np.save(
        out_dir / "SUBSET_REPRESENTATIONS.npy",
        matrix_first,
        allow_pickle=False,
    )
    subset_meta.to_csv(
        out_dir / "SUBSET_METADATA.csv",
        index=False,
    )
    np.save(
        out_dir / "CANDIDATE_R_HO.npy",
        r_ho,
        allow_pickle=False,
    )
    candidate_meta.to_csv(
        out_dir / "CANDIDATE_R_HO_METADATA.csv",
        index=False,
    )

    model_manifest.update(
        {
            "architecture": "EsmForMaskedLM",
            "representation_layer": "last_hidden_state",
            "pooling": "mean_over_full_candidate_mutation_positions",
            "hidden_size": hidden_size,
            "eval_mode": True,
            "fine_tuned": False,
            "gradients_enabled": False,
            "max_position_embeddings": max_position_embeddings,
            "maximum_sequence_length_observed": max_sequence_length,
        }
    )
    (out_dir / "MODEL_MANIFEST.json").write_text(
        json.dumps(model_manifest, indent=2),
        encoding="utf-8",
    )

    runtime = {
        "version": "NABU_FROZEN_EMBEDDING_EXTRACTION_V1",
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
    (out_dir / "RUNTIME_MANIFEST.json").write_text(
        json.dumps(runtime, indent=2),
        encoding="utf-8",
    )

    diagnostics = {
        "selected_candidates": int(len(selected)),
        "subset_sequence_rows": int(len(items)),
        "representation_hidden_size": hidden_size,
        "candidate_r_ho_rows": int(len(candidate_meta)),
        "candidate_r_ho_shape": [
            int(x) for x in r_ho.shape
        ],
        "representation_replay_max_abs_difference": max_abs,
        "deterministic_replay_exact": bool(max_abs == 0.0),
        "subset_representation_l2_mean": float(
            np.mean(np.linalg.norm(matrix_first, axis=1))
        ),
        "candidate_r_ho_l2_mean": float(
            np.mean(np.linalg.norm(r_ho, axis=1))
        ),
        "candidate_r_ho_l2_std": float(
            np.std(np.linalg.norm(r_ho, axis=1))
        ),
        "nucb_consumed": False,
        "phase3_opened": False,
    }
    (out_dir / "REPRESENTATION_DIAGNOSTICS.json").write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    hashes = {
        path.name: sha256_file(path)
        for path in sorted(out_dir.iterdir())
        if path.is_file() and path.name != "REPRESENTATION_HASHES.json"
    }
    (out_dir / "REPRESENTATION_HASHES.json").write_text(
        json.dumps(hashes, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

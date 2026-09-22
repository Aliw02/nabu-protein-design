"""
NABU Protein Design — Standalone Runner
=======================================
Ultra-Fast Epistatic Combinatorial Protein Design Pipeline.

Usage:
    python run_nabu.py PHOT_CHLRE_Chen_2023.csv --out results_nabu
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Ensure 'src' is in sys.path for direct script execution and IDE resolution
sys_path_src = str(Path(__file__).resolve().parent / "src")
if sys_path_src not in sys.path:
    sys.path.insert(0, sys_path_src)

from nabu_protein.core import (
    NabuProteinModel,
    ComponentMemoryModel,
    parse_mutations,
    evaluate_metrics,
    fnv1a_hash,
)


def main():
    parser = argparse.ArgumentParser(description="NABU Protein Design Engine")
    parser.add_argument("csv", default="PHOT_CHLRE_Chen_2023.csv", nargs="?", help="Path to DMS CSV dataset")
    parser.add_argument("--out", default="results_nabu", help="Output directory")
    parser.add_argument("--aggregation", choices=["sum", "average"], default="sum", help="Epistasis aggregation mode")
    parser.add_argument("--seed", type=int, default=161, help="Random seed")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"[+] Loading dataset: {args.csv}")
    df = pd.read_csv(args.csv)
    
    mut_col = next(c for c in ["mutant", "mutation", "mutations"] if c in df.columns)
    score_col = next(c for c in ["DMS_score", "score", "fitness", "mean"] if c in df.columns)
    
    df = df[[mut_col, score_col]].dropna()
    df["ms"] = df[mut_col].map(parse_mutations)
    df = df[df["ms"].map(lambda x: 2 <= len(x) <= 5)]
    df["b"] = df["ms"].map(lambda x: fnv1a_hash(":".join(x)) % 10)

    visible_df = df[df["b"] <= 6].copy()
    hidden_df = df[df["b"] >= 7].copy()

    print(f"[*] Filtered variants: {len(df):,} | Visible: {len(visible_df):,} | Hidden: {len(hidden_df):,}")

    # Train NABU & Baseline
    nabu = NabuProteinModel(aggregation=args.aggregation)
    nabu.fit(visible_df["ms"].tolist(), visible_df[score_col].tolist())

    base = ComponentMemoryModel()
    base.fit(visible_df["ms"].tolist(), visible_df[score_col].tolist())

    # Filter eligible
    eligible = hidden_df[hidden_df["ms"].map(lambda x: all(m in base.component_effects for m in x))].copy()
    print(f"[*] Eligible hidden variants: {len(eligible):,}")

    eligible["component_score"] = base.predict(eligible["ms"].tolist())
    eligible["nabu_score"] = nabu.predict(eligible["ms"].tolist())

    # Candidate manifest
    candidates = eligible.sort_values("nabu_score", ascending=False)[[
        mut_col, "component_score", "nabu_score"
    ]].head(50)
    manifest_path = out_dir / "candidate_manifest_pre_reveal.csv"
    candidates.to_csv(manifest_path, index=False)
    
    sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (out_dir / "candidate_manifest_sha256.txt").write_text(sha + "\n")

    # Metrics
    m_nabu = evaluate_metrics(eligible[score_col].values, eligible["nabu_score"].values)
    m_base = evaluate_metrics(eligible[score_col].values, eligible["component_score"].values)

    # Shuffled control (true random permutation)
    np.random.seed(args.seed)
    shuffled_scores = np.random.permutation(visible_df[score_col].values)
    shuffled_nabu = NabuProteinModel(aggregation=args.aggregation)
    shuffled_nabu.fit(visible_df["ms"].tolist(), shuffled_scores)
    shuf_preds = shuffled_nabu.predict(eligible["ms"].tolist())
    shuf_corr = float(spearmanr(shuf_preds, eligible[score_col].values).statistic)

    R = {
        "architecture": "NABU (Structured State Memory + Pairwise Epistatic Residuals)",
        "rows": len(df),
        "visible": len(visible_df),
        "hidden": len(hidden_df),
        "eligible_hidden": len(eligible),
        "manifest_sha256": sha,
        "nabu_metrics": m_nabu,
        "component_memory": m_base,
        "shuffled_control_spearman": shuf_corr,
        "preregistered_passed": bool(m_nabu["spearman"] > m_base["spearman"] and abs(shuf_corr) < 0.05)
    }

    (out_dir / "results.json").write_text(json.dumps(R, indent=2))
    eligible.drop(columns=["ms"]).to_csv(out_dir / "revealed_predictions.csv", index=False)

    print("\n" + "=" * 55)
    print("NABU BENCHMARK EXECUTION COMPLETE")
    print("=" * 55)
    print(f"• Spearman Correlation  : {m_nabu['spearman']:.4f} (Base: {m_base['spearman']:.4f})")
    print(f"• Top-1 Percentile      : {m_nabu['top1_percentile']*100:.2f}%")
    print(f"• Top-5 Mean Percentile : {m_nabu['top5_mean_percentile']*100:.2f}%")
    print(f"• Top-50 Enrichment     : {m_nabu['top50_enrichment']:.2f}x")
    print(f"• Shuffle Spearman      : {shuf_corr:.4f} (Threshold: < 0.05)")
    print(f"• Preregistration Gate  : {'PASSED [OK]' if R['preregistered_passed'] else 'FAILED'}")
    print("=" * 55)
    print(f"Results saved to: {out_dir.absolute()}\n")


if __name__ == "__main__":
    main()

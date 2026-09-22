"""
NABU Command-Line Interface (CLI)
=================================
CLI commands for training, scoring, and benchmarking protein variant libraries.
"""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .core import NabuProteinModel, ComponentMemoryModel, parse_mutations, evaluate_metrics, fnv1a_hash


def main():
    parser = argparse.ArgumentParser(
        description="NABU Protein Design: Ultra-fast combinatorial protein fitness prediction"
    )
    parser.add_argument("csv", help="Path to DMS CSV dataset (must contain mutant and fitness/DMS_score columns)")
    parser.add_argument("--out", default="results_nabu", help="Output directory (default: results_nabu)")
    parser.add_argument("--aggregation", choices=["sum", "average"], default="sum", help="Epistasis aggregation mode")
    parser.add_argument("--min-mutations", type=int, default=2, help="Minimum mutations per variant (default: 2)")
    parser.add_argument("--max-mutations", type=int, default=5, help="Maximum mutations per variant (default: 5)")
    parser.add_argument("--seed", type=int, default=161, help="Random seed for controls (default: 161)")

    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"🧬 Loading dataset: {args.csv}")
    df = pd.read_csv(args.csv)
    
    # Identify mutant and fitness columns
    mut_col = next(c for c in ["mutant", "mutation", "mutations"] if c in df.columns)
    score_col = next(c for c in ["DMS_score", "score", "fitness", "mean"] if c in df.columns)
    
    df = df[[mut_col, score_col]].dropna()
    df["ms"] = df[mut_col].map(parse_mutations)
    df = df[df["ms"].map(lambda x: args.min_mutations <= len(x) <= args.max_mutations)]
    df["b"] = df["ms"].map(lambda x: fnv1a_hash(":".join(x)) % 10)

    # 70% visible, 30% hidden split
    visible_df = df[df["b"] <= 6].copy()
    hidden_df = df[df["b"] >= 7].copy()

    print(f"📊 Total filtered variants: {len(df):,}")
    print(f"   ├─ Visible (Train): {len(visible_df):,}")
    print(f"   └─ Hidden (Test): {len(hidden_df):,}")

    # Train NABU
    model = NabuProteinModel(aggregation=args.aggregation)
    model.fit(visible_df["ms"].tolist(), visible_df[score_col].tolist())

    # Component Memory baseline
    comp_model = ComponentMemoryModel()
    comp_model.fit(visible_df["ms"].tolist(), visible_df[score_col].tolist())

    # Filter eligible hidden variants
    eligible_hidden = hidden_df[
        hidden_df["ms"].map(lambda x: all(m in comp_model.component_effects for m in x))
    ].copy()
    print(f"🎯 Eligible hidden variants evaluated: {len(eligible_hidden):,}")

    # Predict
    eligible_hidden["component_score"] = comp_model.predict(eligible_hidden["ms"].tolist())
    eligible_hidden["nabu_score"] = model.predict(eligible_hidden["ms"].tolist())

    # Generate Candidate Manifest
    candidates = eligible_hidden.sort_values("nabu_score", ascending=False)[[
        mut_col, "component_score", "nabu_score"
    ]].head(50)
    manifest_path = out_dir / "candidate_manifest_pre_reveal.csv"
    candidates.to_csv(manifest_path, index=False)
    
    manifest_sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (out_dir / "candidate_manifest_sha256.txt").write_text(manifest_sha + "\n")

    # Evaluate Metrics
    nabu_metrics = evaluate_metrics(eligible_hidden[score_col].values, eligible_hidden["nabu_score"].values)
    comp_metrics = evaluate_metrics(eligible_hidden[score_col].values, eligible_hidden["component_score"].values)

    # True Shuffled Control
    np.random.seed(args.seed)
    shuffled_scores = np.random.permutation(visible_df[score_col].values)
    shuffled_model = NabuProteinModel(aggregation=args.aggregation)
    shuffled_model.fit(visible_df["ms"].tolist(), shuffled_scores)
    shuffled_preds = shuffled_model.predict(eligible_hidden["ms"].tolist())
    shuffle_spearman = float(spearmanr(shuffled_preds, eligible_hidden[score_col].values).statistic)

    results = {
        "dataset": args.csv,
        "total_rows": len(df),
        "visible_rows": len(visible_df),
        "hidden_rows": len(hidden_df),
        "eligible_hidden_rows": len(eligible_hidden),
        "manifest_sha256": manifest_sha,
        "nabu_metrics": nabu_metrics,
        "component_memory_metrics": comp_metrics,
        "shuffle_control_spearman": shuffle_spearman,
        "validation_passed": bool(
            nabu_metrics["spearman"] > comp_metrics["spearman"] and
            abs(shuffle_spearman) < 0.05
        )
    }

    results_json = out_dir / "results.json"
    results_json.write_text(json.dumps(results, indent=2))
    
    revealed_path = out_dir / "revealed_predictions.csv"
    eligible_hidden.drop(columns=["ms"]).to_csv(revealed_path, index=False)

    print("\n" + "=" * 50)
    print("🏆 NABU BENCHMARK RESULTS")
    print("=" * 50)
    print(f"• Spearman Correlation : {nabu_metrics['spearman']:.4f} (Additive: {comp_metrics['spearman']:.4f})")
    print(f"• Top-1 Percentile     : {nabu_metrics['top1_percentile']*100:.2f}%")
    print(f"• Top-5 Mean Percentile: {nabu_metrics['top5_mean_percentile']*100:.2f}%")
    print(f"• Top-50 Enrichment    : {nabu_metrics['top50_enrichment']:.2f}x")
    print(f"• Shuffle Control Corr : {shuffle_spearman:.4f} (Target: < 0.05)")
    print(f"• Validation Pass Gate : {'✅ PASSED' if results['validation_passed'] else '❌ FAILED'}")
    print("=" * 50)
    print(f"📁 Output files saved in: {out_dir.absolute()}\n")


if __name__ == "__main__":
    main()

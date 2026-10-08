from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ALL_DEG_FILE = Path(
    r"D:/0work/0wholetransctiptome/2sRNAminic/17HIV1/GEO_HIV_bulk/GSE289893_HIV_vs_Control_24h_DESeq2_results.csv"
)
TARGET_FILE = Path(
    r"D:/0work/0wholetransctiptome/2sRNAminic/17HIV1/intersect_symbols.csv"
)
OUT_DIR = Path(__file__).resolve().parent
N_PERM = 1000
SEED = 20260621


def main() -> None:
    all_df = pd.read_csv(ALL_DEG_FILE)
    target_df = pd.read_csv(TARGET_FILE)

    targets = (
        target_df["x"].astype(str).str.strip().replace("", np.nan).dropna().drop_duplicates()
    )

    all_df = all_df.dropna(subset=["log2FoldChange"]).copy()
    all_df["GeneKey_clean"] = all_df["GeneKey"].astype(str).str.strip()
    all_df["up"] = pd.to_numeric(all_df["log2FoldChange"], errors="coerce") > 0

    target_set = set(targets)
    matched = all_df[all_df["GeneKey_clean"].isin(target_set)].copy()
    missing = sorted(target_set - set(matched["GeneKey_clean"]))

    observed_n = len(matched)
    observed_up = int(matched["up"].sum())
    observed_prop = observed_up / observed_n

    rng = np.random.default_rng(SEED)
    random_counts = np.empty(N_PERM, dtype=int)
    random_props = np.empty(N_PERM)
    for i in range(N_PERM):
        sample = all_df.sample(
            n=observed_n,
            replace=False,
            random_state=int(rng.integers(0, 2**32 - 1)),
        )
        random_counts[i] = int(sample["up"].sum())
        random_props[i] = random_counts[i] / observed_n

    permutation_p = (np.sum(random_props >= observed_prop) + 1) / (N_PERM + 1)
    background_up = int(all_df["up"].sum())
    hypergeom_p = hypergeom.sf(observed_up - 1, len(all_df), background_up, observed_n)

    pd.DataFrame(
        {
            "iteration": np.arange(1, N_PERM + 1),
            "random_up_count": random_counts,
            "random_up_prop": random_props,
        }
    ).to_csv(OUT_DIR / "hiv_target_up_permutation_distribution.csv", index=False)

    pd.DataFrame({"missing_target_symbol": missing}).to_csv(
        OUT_DIR / "hiv_target_missing_symbols.csv", index=False
    )

    summary_lines = [
        f"Input target symbols: {len(targets)}",
        f"Matched target genes in DESeq2 table: {observed_n}",
        f"Missing target symbols: {len(missing)}",
        f"Observed up genes: {observed_up}/{observed_n}",
        f"Observed up proportion: {observed_prop:.6f}",
        f"Background up genes: {background_up}/{len(all_df)}",
        f"Background up proportion: {background_up / len(all_df):.6f}",
        f"Random mean up proportion ({N_PERM} permutations): {random_props.mean():.6f}",
        f"Random SD up proportion: {random_props.std(ddof=1):.6f}",
        f"Random 2.5%-97.5% interval: {np.quantile(random_props, 0.025):.6f}, {np.quantile(random_props, 0.975):.6f}",
        f"One-sided permutation p-value P(random >= observed): {permutation_p:.6f}",
        f"One-sided hypergeometric p-value: {hypergeom_p:.6g}",
    ]
    (OUT_DIR / "hiv_target_up_permutation_summary.txt").write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=160)
    ax.hist(random_props, bins=30, color="#8fb7d9", edgecolor="white")
    ax.axvline(observed_prop, color="#c43b3b", linewidth=2, label="HIV targets")
    ax.set_xlabel("Up-regulated proportion in random gene sets")
    ax.set_ylabel("Permutation count")
    ax.set_title("HIV target genes are enriched for up-regulation")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "hiv_target_up_permutation_hist.pdf")


if __name__ == "__main__":
    main()

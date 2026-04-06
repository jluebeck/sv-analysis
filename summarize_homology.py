#!/usr/bin/env python3
"""
summarize_homology.py — Summarize homology-detection improvement over AA baseline.

Loads final_augmented.tsv files produced by SVRecalibrator (refine.py --mode both)
from a batch output directory and reports how many SVs gained microhomology by
split-read and/or scaffold approaches that were not present in the AA input calls.

Usage:
    python summarize_homology.py <outdir> [--sample SAMPLE ...] [--min-hom-len N]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_sp_hom_len(val):
    """Convert sp_hom_len (stored as string) to int or NaN."""
    if pd.isna(val) or val == "" or val == "nan":
        return np.nan
    if str(val).strip() == "N/A":
        return 0
    try:
        return int(val)
    except (ValueError, TypeError):
        return np.nan


def load_augmented(outdir, samples=None):
    """Load all final_augmented.tsv files under outdir/*/final_augmented.tsv."""
    frames = []
    missing = []

    candidate_dirs = (
        [os.path.join(outdir, s) for s in samples]
        if samples
        else [
            os.path.join(outdir, d)
            for d in sorted(os.listdir(outdir))
            if os.path.isdir(os.path.join(outdir, d))
        ]
    )

    for sample_dir in candidate_dirs:
        sample = os.path.basename(sample_dir)
        tsv = os.path.join(sample_dir, "final_augmented.tsv")
        if not os.path.isfile(tsv):
            missing.append(sample)
            continue
        try:
            df = pd.read_csv(tsv, sep="\t", low_memory=False)
        except EmptyDataError:
            continue  # empty file — sample produced no SVs
        except Exception as e:
            print(f"Warning: could not load {tsv}: {e}", file=sys.stderr)
            missing.append(sample)
            continue
        if df.empty:
            continue
        df["_sample"] = sample
        frames.append(df)

    if missing:
        print(f"Warning: no final_augmented.tsv for {len(missing)} sample(s): {', '.join(missing)}",
              file=sys.stderr)
    if not frames:
        sys.exit("No data loaded — check outdir and sample names.")
    return pd.concat(frames, ignore_index=True)


def classify(df, min_hom_len):
    """Add boolean columns for AA / split-read / scaffold homology presence."""
    # AA homology (optional column)
    if "homology_len" in df.columns:
        df["aa_hom"] = pd.to_numeric(df["homology_len"], errors="coerce").fillna(0) >= min_hom_len
    else:
        df["aa_hom"] = False

    # Split-read homology: sp_hom_len is a string; positive = homology
    df["_sp_hom_len_int"] = df["sp_hom_len"].apply(parse_sp_hom_len)
    df["sp_hom"] = df["_sp_hom_len_int"] >= min_hom_len

    # Scaffold homology: sc_hom_len is numeric; positive = homology
    df["_sc_hom_len_num"] = pd.to_numeric(df["sc_hom_len"], errors="coerce")
    df["sc_hom"] = df["_sc_hom_len_num"] >= min_hom_len

    df["either_hom"] = df["sp_hom"] | df["sc_hom"]
    return df


def pct(n, total):
    return f"{n:>5d}  ({100 * n / total:.1f}%)" if total else f"{n:>5d}  (  n/a)"


def print_section(title):
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


def summarize(df, min_hom_len):
    total = len(df)
    has_aa = "homology_len" in df.columns

    print_section("Dataset overview")
    samples = df["_sample"].nunique()
    print(f"  Samples loaded           : {samples}")
    print(f"  Total SVs                : {total}")
    if has_aa:
        n_aa = df["aa_hom"].sum()
        print(f"  SVs with AA homology     : {pct(n_aa, total)}")
        print(f"  SVs without AA homology  : {pct(total - n_aa, total)}")
    else:
        print("  (homology_len column absent — AA homology comparison unavailable)")

    # ── SVRecalibrator detection across ALL SVs ────────────────────────────────
    print_section("SVRecalibrator homology detection (all SVs)")
    n_sp   = df["sp_hom"].sum()
    n_sc   = df["sc_hom"].sum()
    n_both = (df["sp_hom"] & df["sc_hom"]).sum()
    n_any  = df["either_hom"].sum()
    n_sp_only = (df["sp_hom"] & ~df["sc_hom"]).sum()
    n_sc_only = (~df["sp_hom"] & df["sc_hom"]).sum()

    print(f"  Detected by split reads  : {pct(n_sp,   total)}")
    print(f"  Detected by scaffold     : {pct(n_sc,   total)}")
    print(f"  Detected by both         : {pct(n_both, total)}")
    print(f"  Detected by either       : {pct(n_any,  total)}")
    print(f"    split only             : {pct(n_sp_only, total)}")
    print(f"    scaffold only          : {pct(n_sc_only, total)}")

    # ── Gain over AA ──────────────────────────────────────────────────────────
    if has_aa:
        print_section(f"Homology gain over AA (SVs without AA homology, min_hom_len={min_hom_len})")
        no_aa = df[~df["aa_hom"]]
        N = len(no_aa)
        print(f"  Universe (no AA hom)     : {N}")
        if N:
            g_sp      = no_aa["sp_hom"].sum()
            g_sc      = no_aa["sc_hom"].sum()
            g_both    = (no_aa["sp_hom"] & no_aa["sc_hom"]).sum()
            g_any     = no_aa["either_hom"].sum()
            g_sp_only = (no_aa["sp_hom"] & ~no_aa["sc_hom"]).sum()
            g_sc_only = (~no_aa["sp_hom"] & no_aa["sc_hom"]).sum()
            g_none    = (~no_aa["either_hom"]).sum()

            print(f"  Gained by either         : {pct(g_any,     N)}")
            print(f"    split reads only       : {pct(g_sp_only, N)}")
            print(f"    scaffold only          : {pct(g_sc_only, N)}")
            print(f"    both approaches        : {pct(g_both,    N)}")
            print(f"  Still no homology        : {pct(g_none,    N)}")

        # Concordance on SVs that DO have AA homology
        print_section("Concordance on SVs with AA homology")
        with_aa = df[df["aa_hom"]]
        M = len(with_aa)
        print(f"  Universe (AA hom)        : {M}")
        if M:
            c_any  = with_aa["either_hom"].sum()
            c_none = M - c_any
            print(f"  Confirmed by SVRecalib   : {pct(c_any,  M)}")
            print(f"  Not confirmed            : {pct(c_none, M)}")

    # ── Homology length distributions ─────────────────────────────────────────
    print_section("Homology length distributions (detected SVs only)")

    def len_stats(series, label):
        s = pd.to_numeric(series, errors="coerce").dropna()
        s = s[s > 0]
        if s.empty:
            print(f"  {label}: no data")
            return
        print(f"  {label} (n={len(s)}):")
        print(f"    median={s.median():.0f}  mean={s.mean():.1f}  "
              f"min={s.min():.0f}  max={s.max():.0f}  "
              f"p25={s.quantile(0.25):.0f}  p75={s.quantile(0.75):.0f}")

    len_stats(df.loc[df["sp_hom"], "_sp_hom_len_int"], "split-read hom len")
    len_stats(df.loc[df["sc_hom"], "_sc_hom_len_num"], "scaffold   hom len")
    if has_aa:
        len_stats(df.loc[df["aa_hom"], "homology_len"], "AA         hom len")

    # ── Per-sample breakdown ───────────────────────────────────────────────────
    print_section("Per-sample breakdown")
    grp = df.groupby("_sample")
    rows = []
    for sample, g in grp:
        n = len(g)
        aa  = g["aa_hom"].sum() if has_aa else np.nan
        sp  = g["sp_hom"].sum()
        sc  = g["sc_hom"].sum()
        any_ = g["either_hom"].sum()
        no_aa_g = g[~g["aa_hom"]] if has_aa else g
        gain = no_aa_g["either_hom"].sum() if has_aa else np.nan
        rows.append({"sample": sample, "SVs": n, "AA_hom": aa,
                     "sp_hom": sp, "sc_hom": sc, "either": any_, "new_gain": gain})
    tbl = pd.DataFrame(rows).set_index("sample")
    if not has_aa:
        tbl = tbl.drop(columns=["AA_hom", "new_gain"])
    print(tbl.to_string())
    print()
    return tbl


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="batch output directory (contains one subdirectory per sample)")
    p.add_argument("--sample", nargs="+", metavar="NAME", help="restrict to these sample(s)")
    p.add_argument(
        "--min-hom-len",
        type=int,
        default=1,
        metavar="N",
        help="minimum homology length to count as detected (default: 1)",
    )
    p.add_argument(
        "--out-csv",
        metavar="FILE",
        help="write per-sample summary table to this CSV file",
    )
    args = p.parse_args()

    if not os.path.isdir(args.outdir):
        sys.exit(f"outdir not found: {args.outdir}")

    df = load_augmented(args.outdir, args.sample)
    df = classify(df, args.min_hom_len)
    tbl = summarize(df, args.min_hom_len)
    if args.out_csv:
        tbl.to_csv(args.out_csv)
        print(f"Per-sample table written to {args.out_csv}")


if __name__ == "__main__":
    main()

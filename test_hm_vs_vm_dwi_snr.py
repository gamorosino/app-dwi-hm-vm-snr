#!/usr/bin/env python3
"""
Paired statistical test of DWI SNR: Horizontal Meridian (HM) vs Vertical
Meridian (UVM+LVM), using the CSV produced by compute_dwi_snr_hm_vm.py
(columns: subject, meridian, mean_snr, voxel_count).

Unlike voxel-count/volume metrics (which are summed across UVM+LVM), SNR is
a ratio and is combined across UVM+LVM with a voxel-count-weighted average.
This is exact (not just a convenient approximation) as long as HM/UVM/LVM
SNR values share the same noise estimate, which compute_dwi_snr_hm_vm.py
guarantees by reusing one noise mask across all three ROIs per subject.

Usage:
    python test_hm_vs_vm_dwi_snr.py hm_vm_dwi_snr.csv
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats

HM_LABELS = ["HM"]
VM_LABELS = ["UVM", "LVM"]


def weighted_group_snr(df, meridians):
    sub = df[df["meridian"].isin(meridians)].dropna(subset=["mean_snr"])

    def combine(g):
        weights = g["voxel_count"].to_numpy()
        values = g["mean_snr"].to_numpy()
        if weights.sum() == 0:
            return np.nan
        return np.average(values, weights=weights)

    return sub.groupby("subject")[["mean_snr", "voxel_count"]].apply(combine)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    required = {"subject", "meridian", "mean_snr", "voxel_count"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"CSV missing columns: {missing}")

    hm = weighted_group_snr(df, HM_LABELS)
    vm = weighted_group_snr(df, VM_LABELS)

    per_subj = pd.DataFrame({"HM": hm, "VM": vm}).dropna()
    n = len(per_subj)
    if n < 2:
        sys.exit(f"Not enough subjects with both HM and VM SNR ({n}) to run a paired test.")

    hm_vals = per_subj["HM"].to_numpy()
    vm_vals = per_subj["VM"].to_numpy()
    diff = hm_vals - vm_vals

    t_stat, t_p = stats.ttest_rel(hm_vals, vm_vals)
    w_stat, w_p = stats.wilcoxon(hm_vals, vm_vals)
    cohens_d = diff.mean() / diff.std(ddof=1)

    print(f"Subjects with both HM and VM SNR: {n}")
    print()
    print(f"HM SNR:   mean={hm_vals.mean():.3f}  sd={hm_vals.std(ddof=1):.3f}  median={np.median(hm_vals):.3f}")
    print(f"VM SNR:   mean={vm_vals.mean():.3f}  sd={vm_vals.std(ddof=1):.3f}  median={np.median(vm_vals):.3f}")
    print(f"Diff (HM-VM): mean={diff.mean():.3f}  sd={diff.std(ddof=1):.3f}  median={np.median(diff):.3f}")
    print()
    print("Paired t-test:        t = {:.4f}, p = {:.3e}".format(t_stat, t_p))
    print("Wilcoxon signed-rank: W = {:.4f}, p = {:.3e}".format(w_stat, w_p))
    print(f"Cohen's d (paired):   {cohens_d:.4f}")
    print()

    alpha = 0.05
    sig = "significantly" if t_p < alpha else "not significantly"
    direction = "greater than" if diff.mean() > 0 else "less than"
    print(f"=> HM DWI SNR is {sig} different from VM DWI SNR (paired t-test p={t_p:.3e}, alpha={alpha}).")
    print(f"   On average HM SNR is {direction} VM SNR.")


if __name__ == "__main__":
    main()

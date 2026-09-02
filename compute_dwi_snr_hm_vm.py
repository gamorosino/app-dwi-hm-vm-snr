#!/usr/bin/env python3
"""
Compute DWI signal-to-noise ratio (SNR) within the Horizontal Meridian (HM)
and Vertical Meridian (UVM/LVM) cortical ROIs, per subject.

Meridian ROI definition mirrors compute_hm_vm_from_raw_retinotopy.py in
VISCONTI_analysis (both hemispheres combined, restricted to any Benson14
visual area and eccentricity 0-90 deg):
    HM  : polar angle 75-105 deg
    UVM : polar angle 0-15 deg
    LVM : polar angle 165-180 deg

SNR is computed with scilpy's compute_snr (mean DWI signal inside the ROI
mask divided by the std of a background-noise mask, estimated once via
median_otsu on the b0 volume(s) unless --noise-mask is given), averaged
across the b0 volume(s). All three ROIs share the same noise estimate, so
their SNR values are directly comparable and voxel-count-weighted averaging
across ROIs is valid.

The polarAngle/eccentricity/varea maps (from the "prf" Benson/neuropythy
output) typically live on the FreeSurfer-conformed grid (e.g. 256^3 @ 1mm),
which differs from the DWI grid. If their affine puts them in the same
physical (scanner RAS) space as the DWI volume -- true when both derive
from the same T1w, as in the VISCONTI pipeline -- this script resamples them
onto the DWI grid with nearest-neighbor interpolation before building the
ROI masks, so no separate registration step is required. If your data is in
a genuinely different space (a different subject/session, or unregistered),
warp the retinotopy maps into DWI space yourself first.

Single subject:
    python compute_dwi_snr_hm_vm.py \\
        --subject 100206 \\
        --dwi dwi.nii.gz --bval dwi.bvals --bvec dwi.bvecs \\
        --polar-angle polarAngle.nii.gz \\
        --eccentricity eccentricity.nii.gz \\
        --varea varea.nii.gz \\
        --out-csv snr_100206.csv --out-dir ./work_100206

Batch over subjects under a base directory (mirrors
compute_hm_vm_from_raw_retinotopy.py --base-dir usage):
    python compute_dwi_snr_hm_vm.py \\
        --base-dir ~/data/VISCONTI \\
        --retino-subdir prf \\
        --dwi-rel dwi/dwi.nii.gz --bval-rel dwi/dwi.bvals --bvec-rel dwi/dwi.bvecs \\
        --output-csv hm_vm_dwi_snr.csv --out-dir ./work
"""

import argparse
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to

try:
    from scilpy.utils.image import compute_snr
except ImportError:
    from scilpy.image.utils import compute_snr  # newer scilpy layout

MERIDIAN_ANGLE_RANGES = {
    "HM": (75, 105),
    "UVM": (0, 15),
    "LVM": (165, 180),
}


def _load_resampled(path, target_img):
    """Load a (possibly 4D-with-singleton) map and resample (nearest
    neighbor) onto target_img's grid if the grids differ."""
    img = nib.load(str(path))
    data = np.asarray(img.dataobj).squeeze()
    img = nib.Nifti1Image(data, img.affine)  # drop any singleton 4th dim
    target_shape_affine = (target_img.shape[:3], target_img.affine)
    if img.shape[:3] != target_img.shape[:3] or not np.allclose(img.affine, target_img.affine):
        img = resample_from_to(img, target_shape_affine, order=0)
    return img.get_fdata().squeeze()


def compute_subject_snr(subject, dwi, bval, bvec, polar_angle, eccentricity,
                         varea, out_dir, b0_thr=50.0, noise_mask=None,
                         varea_labels=None):
    """Compute per-meridian mean DWI SNR for one subject.

    varea_labels: optional iterable of integer varea labels to restrict the
        ROI to (e.g. [1] for V1 only). Default (None) uses any label > 0
        (all Benson14 visual areas), matching compute_hm_vm_from_raw_retinotopy.py.

    Returns a list of dicts with keys: subject, meridian, mean_snr,
    n_b0_used, voxel_count.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dwi_img = nib.load(str(dwi))
    angle = _load_resampled(polar_angle, dwi_img)
    ecc = _load_resampled(eccentricity, dwi_img)
    varea_data = _load_resampled(varea, dwi_img)
    affine = dwi_img.affine

    if varea_labels:
        in_varea = np.isin(varea_data, list(varea_labels))
    else:
        in_varea = varea_data > 0
    in_ecc = (ecc >= 0) & (ecc <= 90)

    basename = str(out_dir / f"{subject}_snr")
    shared_noise_mask = noise_mask

    rows = []
    for meridian, (lo, hi) in MERIDIAN_ANGLE_RANGES.items():
        roi_mask = in_varea & in_ecc & (angle >= lo) & (angle <= hi)
        voxel_count = int(roi_mask.sum())

        if voxel_count == 0:
            print(f"[WARN] {subject}: no voxels for meridian {meridian}, skipping")
            rows.append(dict(subject=subject, meridian=meridian,
                              mean_snr=np.nan, n_b0_used=0, voxel_count=0))
            continue

        roi_mask_path = out_dir / f"{subject}_{meridian}_mask.nii.gz"
        nib.save(nib.Nifti1Image(roi_mask.astype(np.uint8), affine), roi_mask_path)

        values = compute_snr(
            str(dwi), str(bval), str(bvec), b0_thr, str(roi_mask_path),
            noise_mask=shared_noise_mask, basename=basename,
        )

        if shared_noise_mask is None:
            # compute_snr wrote its auto-derived background-noise mask here,
            # as float32. Its own loader (get_data_as_mask) rejects float
            # masks, so re-cast to uint8 before reusing it for the remaining
            # meridians (this keeps the same noise floor shared across all
            # three ROIs).
            shared_noise_mask = basename + "_noise_mask.nii.gz"
            noise_img = nib.load(shared_noise_mask)
            nib.save(
                nib.Nifti1Image((noise_img.get_fdata() > 0).astype(np.uint8), noise_img.affine),
                shared_noise_mask,
            )

        b0_snrs = [v["snr"] for v in values.values() if v["bval"] <= b0_thr]
        rows.append(dict(
            subject=subject,
            meridian=meridian,
            mean_snr=float(np.mean(b0_snrs)),
            n_b0_used=len(b0_snrs),
            voxel_count=voxel_count,
        ))

    return rows


def find_subjects(base_dir, retino_subdir, dwi_rel, bval_rel, bvec_rel):
    subjects = []
    for subj_dir in sorted(Path(base_dir).iterdir()):
        if not subj_dir.is_dir():
            continue
        retino_dir = subj_dir / retino_subdir
        required = [
            retino_dir / "polarAngle.nii.gz",
            retino_dir / "eccentricity.nii.gz",
            retino_dir / "varea.nii.gz",
            subj_dir / dwi_rel,
            subj_dir / bval_rel,
            subj_dir / bvec_rel,
        ]
        if all(p.exists() for p in required):
            subjects.append(subj_dir.name)
    return subjects


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--subject", help="Single-subject mode: subject label")
    parser.add_argument("--dwi", help="Single-subject mode: path to DWI nifti")
    parser.add_argument("--bval", help="Single-subject mode: path to bvals")
    parser.add_argument("--bvec", help="Single-subject mode: path to bvecs")
    parser.add_argument("--polar-angle", help="Single-subject mode: polarAngle.nii.gz")
    parser.add_argument("--eccentricity", help="Single-subject mode: eccentricity.nii.gz")
    parser.add_argument("--varea", help="Single-subject mode: varea.nii.gz")

    parser.add_argument("--base-dir", help="Batch mode: root dir containing one subdir per subject")
    parser.add_argument("--retino-subdir", default="prf",
                         help="Batch mode: subject-relative dir with polarAngle/eccentricity/varea (default: prf)")
    parser.add_argument("--dwi-rel", default="dwi/dwi.nii.gz", help="Batch mode: subject-relative DWI path")
    parser.add_argument("--bval-rel", default="dwi/dwi.bvals", help="Batch mode: subject-relative bvals path")
    parser.add_argument("--bvec-rel", default="dwi/dwi.bvecs", help="Batch mode: subject-relative bvecs path")
    parser.add_argument("--limit", type=int, default=None, help="Batch mode: process only first N subjects")

    parser.add_argument("--varea-labels", default=None,
                         help="Comma-separated varea label(s) to restrict the ROI to (e.g. '1' for V1 only). "
                              "Default: any label > 0 (all Benson14 visual areas).")
    parser.add_argument("--b0-thr", type=float, default=50.0,
                         help="b-value threshold below which a volume is treated as a b0 (default: 50)")
    parser.add_argument("--noise-mask", default=None,
                         help="Optional precomputed background-noise mask, shared across subjects. "
                              "If omitted, a noise mask is auto-estimated per subject from the b0(s).")
    parser.add_argument("--out-dir", default="./dwi_snr_work",
                         help="Directory for intermediate ROI/noise masks")
    parser.add_argument("--out-csv", default=None, help="Single-subject mode: output CSV path")
    parser.add_argument("--output-csv", default="hm_vm_dwi_snr.csv", help="Batch mode: output CSV path")
    args = parser.parse_args()

    varea_labels = ([int(v) for v in args.varea_labels.split(",")]
                     if args.varea_labels else None)

    single_mode = args.subject is not None
    if single_mode:
        missing = [name for name in ("dwi", "bval", "bvec", "polar_angle", "eccentricity", "varea")
                   if getattr(args, name) is None]
        if missing:
            sys.exit(f"--subject requires also: {['--' + m.replace('_', '-') for m in missing]}")

        rows = compute_subject_snr(
            args.subject, args.dwi, args.bval, args.bvec,
            args.polar_angle, args.eccentricity, args.varea,
            out_dir=args.out_dir, b0_thr=args.b0_thr, noise_mask=args.noise_mask,
            varea_labels=varea_labels,
        )
        df = pd.DataFrame(rows)
        out_csv = args.out_csv or f"{args.subject}_dwi_snr.csv"
        df.to_csv(out_csv, index=False)
        print(df.to_string(index=False))
        print(f"\nSaved to {out_csv}")
        return

    if not args.base_dir:
        sys.exit("Provide either --subject (single-subject mode) or --base-dir (batch mode)")

    base_dir = Path(args.base_dir).expanduser()
    subjects = find_subjects(base_dir, args.retino_subdir, args.dwi_rel, args.bval_rel, args.bvec_rel)
    if args.limit:
        subjects = subjects[: args.limit]
    if not subjects:
        sys.exit(
            f"No subjects found under {base_dir} with retinotopy subdir "
            f"'{args.retino_subdir}' and DWI files '{args.dwi_rel}'/'{args.bval_rel}'/'{args.bvec_rel}'. "
            f"Adjust --retino-subdir/--dwi-rel/--bval-rel/--bvec-rel to match your layout."
        )

    print(f"Found {len(subjects)} subjects with retinotopy + DWI data.")

    all_rows = []
    for i, subject in enumerate(subjects, 1):
        subj_dir = base_dir / subject
        retino_dir = subj_dir / args.retino_subdir
        out_dir = Path(args.out_dir) / subject
        try:
            rows = compute_subject_snr(
                subject,
                dwi=subj_dir / args.dwi_rel,
                bval=subj_dir / args.bval_rel,
                bvec=subj_dir / args.bvec_rel,
                polar_angle=retino_dir / "polarAngle.nii.gz",
                eccentricity=retino_dir / "eccentricity.nii.gz",
                varea=retino_dir / "varea.nii.gz",
                out_dir=out_dir,
                b0_thr=args.b0_thr,
                noise_mask=args.noise_mask,
                varea_labels=varea_labels,
            )
            all_rows.extend(rows)
        except Exception as exc:  # noqa: BLE001
            print(f"[WARN] skipping {subject}: {exc}", file=sys.stderr)

        if i % 10 == 0 or i == len(subjects):
            print(f"  processed {i}/{len(subjects)} subjects")

    df = pd.DataFrame(all_rows)
    df.to_csv(args.output_csv, index=False)
    print(f"\nSaved {len(df)} rows ({df['subject'].nunique()} subjects) to {args.output_csv}")

    summary = df.groupby("meridian")[["mean_snr", "voxel_count"]].mean()
    print("\nPer-meridian group means (across subjects):")
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()

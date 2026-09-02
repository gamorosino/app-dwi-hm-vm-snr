# DWI SNR: Horizontal vs Vertical Meridian

Checks whether DWI signal-to-noise ratio differs between the Horizontal
Meridian (HM) and Vertical Meridian (UVM+LVM) cortical ROIs used in the
VISCONN retinotopy/connectivity analysis. This is a confound check for the
tractography-based connectivity claims: if HM and VM regions have
systematically different DWI SNR, that alone could produce a spurious
difference in streamline counts between them.

## Scripts

- `compute_dwi_snr_hm_vm.py` — per-subject: builds the HM/UVM/LVM ROI masks
  from `polarAngle.nii.gz` / `eccentricity.nii.gz` / `varea.nii.gz` (same
  angle-range definition as `compute_hm_vm_from_raw_retinotopy.py` in
  `VISCONN_analysis/code/scripts`: HM=75-105deg, UVM=0-15deg, LVM=165-180deg,
  restricted to varea>0 and eccentricity 0-90deg), then computes DWI SNR in
  each ROI with scilpy's `compute_snr` (mean signal in the ROI / std of a
  shared background-noise mask, auto-estimated once per subject via
  `median_otsu`). Supports both a single-subject CLI and a `--base-dir`
  batch mode. Use `--varea-labels 1` to restrict to V1 only instead of all
  visual areas.
- `test_hm_vs_vm_dwi_snr.py` — paired t-test / Wilcoxon / Cohen's d on the
  resulting CSV, combining UVM+LVM into one "VM" value per subject via a
  voxel-count-weighted average (valid because all three ROIs share the same
  noise estimate).

Run `compute_dwi_snr_hm_vm.py` in an environment with `scilpy` + `dipy`
installed (e.g. the `tract_align` conda env); `test_hm_vs_vm_dwi_snr.py`
only needs `pandas`/`scipy`/`numpy`.

## Space handling

`polarAngle`/`eccentricity`/`varea` (the "prf" outputs) typically live on
the FreeSurfer-conformed grid (e.g. 256^3 @ 1mm), which differs from the DWI
grid. `compute_dwi_snr_hm_vm.py` resamples them onto the DWI grid with
nearest-neighbor interpolation before building the ROI masks — this only
works because the conformed grid and the subject's native T1w/DWI grid share
the same physical (scanner RAS) space, which was confirmed on the test
dataset below (varea>0 voxels transform correctly into the occipital pole in
DWI-image coordinates). No separate registration step is needed in this
pipeline. If your retinotopy maps come from a genuinely different
space/session, warp them into DWI space yourself first.

## Test dataset

A one-subject example lives at `~/data/test_DWI_HM_VM_SNR/`:

```
prf/polarAngle.nii.gz, prf/eccentricity.nii.gz, prf/varea.nii.gz   (256^3 @ 1mm)
dwi/dwi.nii.gz, dwi/dwi.bvals, dwi/dwi.bvecs                        (145x174x145x288 @ 1.25mm, 18 b0 + 3 shells)
```

```bash
python compute_dwi_snr_hm_vm.py \
    --subject test001 \
    --dwi ~/data/test_DWI_HM_VM_SNR/dwi/dwi.nii.gz \
    --bval ~/data/test_DWI_HM_VM_SNR/dwi/dwi.bvals \
    --bvec ~/data/test_DWI_HM_VM_SNR/dwi/dwi.bvecs \
    --polar-angle ~/data/test_DWI_HM_VM_SNR/prf/polarAngle.nii.gz \
    --eccentricity ~/data/test_DWI_HM_VM_SNR/prf/eccentricity.nii.gz \
    --varea ~/data/test_DWI_HM_VM_SNR/prf/varea.nii.gz \
    --out-dir ./work --out-csv snr_test001.csv
```

Verified output:

```
subject meridian  mean_snr  n_b0_used  voxel_count
test001       HM 15.520436         18         7753
test001      UVM 13.006310         18         3090
test001      LVM 17.106253         18         4477
```

(~1.5-2 min/subject, dominated by scilpy reloading the full 4D DWI volume
once per ROI.)

## Usage

Single subject: see the test dataset example above.

Batch over subjects (adjust `--retino-subdir`/`--dwi-rel`/`--bval-rel`/
`--bvec-rel` if your per-subject layout differs from `<subject>/prf/...` +
`<subject>/dwi/dwi.{nii.gz,bvals,bvecs}`):

```bash
python compute_dwi_snr_hm_vm.py \
    --base-dir ~/data/VISCONN \
    --out-dir ./work --output-csv hm_vm_dwi_snr.csv

python test_hm_vs_vm_dwi_snr.py hm_vm_dwi_snr.csv
```

<h1 align="center">DWI SNR (HM vs VM)</h1>

<p align="center">
  Computes DWI b0 signal-to-noise ratio in horizontal- and vertical-meridian visual cortex ROIs built from Benson14 retinotopy maps, using scilpy compute_snr.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.956"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.956-blue.svg" alt="Run on brainlife.io"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-see_LICENSE-lightgrey.svg" alt="License"></a>
</p>

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Computes DWI b0 signal-to-noise ratio in horizontal- and vertical-meridian visual cortex ROIs built from Benson14 retinotopy maps, using scilpy compute_snr.**

Builds HM (polar angle 75-105 deg), UVM (0-15 deg) and LVM (165-180 deg) ROIs from the prf polar angle, eccentricity and visual area maps (resampled to the DWI grid), estimates a background noise mask with median_otsu and computes the mean b0 SNR in each ROI. Use it as a confound check that HM and VM regions do not differ in DWI SNR before comparing their streamline counts. The prf maps must be in the same physical space as the DWI.

| | |
|---|---|
| **Inputs** | `neuro/dwi` · `neuro/prf` |
| **Outputs** | `neuro/SNR_and_dMRI_stats` |
| **Container** | `gamorosino/tract_align:v1.1` |
| **Resources** | 4 GB memory · 1 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.956](https://doi.org/10.25663/brainlife.app.956) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `dwi` | `neuro/dwi` | yes | DWI image with bvals and bvecs. |
| `prf` | `neuro/prf` | yes | Retinotopy maps (polarAngle, eccentricity, varea), for example from Benson14, aligned to the DWI. |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `stats` | `neuro/SNR_and_dMRI_stats` | Per-subject mean DWI SNR (b0) in HM/UVM/LVM cortical ROIs, one row per meridian, with voxel counts. |

<details>
<summary><b>Parameters</b> (2)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `varea_labels` | string |  | Comma-separated varea label(s) to restrict the ROI to (e.g. '1' for V1 only). Empty = all Benson14 visual areas (varea>0). |
| `b0_thr` | number | `50` | b-value threshold below which a volume is treated as a b0. |

</details>

<!-- brainlife-app:end -->

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.956](https://doi.org/10.25663/brainlife.app.956), choose `neuro/dwi`, `neuro/prf` from your project, set the parameters and press **Execute**. The results are saved to your project as `neuro/SNR_and_dMRI_stats`.

### Locally

You need [Singularity/Apptainer](https://apptainer.org/) (the App runs in `gamorosino/tract_align:v1.1`), [`jq`](https://jqlang.github.io/jq/).

1. Clone this repository.
2. Create `config.json` with the paths to your inputs and the parameters:

   ```json
   {
       "dwi": "/path/to/dwi/dwi",
       "bvals": "/path/to/dwi/bvals",
       "bvecs": "/path/to/dwi/bvecs",
       "polarAngle": "/path/to/prf/polarAngle",
       "eccentricity": "/path/to/prf/eccentricity",
       "varea": "/path/to/prf/varea",
       "varea_labels": "",
       "b0_thr": 50
   }
   ```

3. Run `./main`.

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

## Citing

If you use this App, please cite it, the software it builds on and brainlife.io:

- DWI SNR (HM vs VM). brainlife.io App. [doi:10.25663/brainlife.app.956](https://doi.org/10.25663/brainlife.app.956)
- Garyfallidis et al. Dipy, a library for the analysis of diffusion MRI data. *Frontiers in Neuroinformatics* 8 (2014). [doi:10.3389/fninf.2014.00008](https://doi.org/10.3389/fninf.2014.00008)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

Maintained by [@gamorosino](https://github.com/gamorosino).

## License

See [LICENSE](LICENSE).

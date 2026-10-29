"""Map saved SVM-RFE selection frequencies to the local Canlab2018 atlas.

Produces a self-contained time-slider fragment, static axial montages, a
single-play GIF, and primary-contrast 4-D NIfTI maps. No decoder is refitted.
The atlas mapping follows the extraction code's ordinal convention; exact
historical atlas-version identity is not independently established.
"""
from pathlib import Path
import argparse
import base64
import gzip
import hashlib
import json

import h5py
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.processing import resample_from_to
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[4]
MODELS = ["cHRF", "cHRFderiv", "sHRF"]
JOBS = {"cHRF": "job_36477504_task_0", "cHRFderiv": "job_36477598_task_1", "sHRF": "job_36477153_task_2"}
TIMES = np.arange(41) * .72
CUTS = [-54, -30, -6, 18, 42, 66]
OVERVIEW_LAGS = [0, 6, 12, 18, 24, 30, 36, 40]


def pack(array):
    return base64.b64encode(gzip.compress(array.tobytes(order="C"), mtime=0)).decode("ascii")


def decode(h, ref):
    return "".join(chr(int(v)) for v in h[ref][()].ravel() if v)


def fingerprint(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return {"path": str(path), "sha256": h.hexdigest()}


def read_frequencies():
    numerators, denominators, all_counts, performance, labels = [], [], [], [], None
    provenance = []
    for model in MODELS:
        path = ROOT / "data/feature_selection/feature_selection" / model / JOBS[model] / "observed_wm_interaction.mat"
        print(f"Reading selection maps: {model}", flush=True)
        with h5py.File(path, "r") as h:
            assert int(h["nBoot"][()].item()) == 100
            current_labels = [decode(h, v) for v in h["contrast_labels"][()].ravel()]
            if labels is None: labels = current_labels
            assert labels == current_labels
            fs = h["FeatureSelection"]
            partial = np.zeros((6, 41, 489), dtype=np.uint8)
            den = np.zeros((6, 41), dtype=np.uint8)
            selected_counts = np.zeros_like(partial)
            for ref in fs["bootstrap"][()].ravel():
                q = h[ref]
                assert q["valid"][()].all()
                assert q["eligible"][()].all(), "Compact all-fit reconstruction requires all 489 eligible."
                selected = q["selected"][()].astype(np.uint8)
                assert np.array_equal(selected.sum(axis=-1), q["count"][()])
                reduced = q["count"][()] < 489
                partial += selected * reduced[..., None]
                den += reduced
                selected_counts += selected
            # Full-atlas fits contribute one inclusion for every ROI. Keeping
            # only partial numerators avoids storing redundant data inline.
            assert np.array_equal(partial + (100 - den)[..., None], selected_counts)
            assert np.array_equal(selected_counts, fs["selected_boot"][()])
            assert np.allclose(selected_counts / 100, fs["selection_frequency"][()])
            numerators.append(partial); denominators.append(den); all_counts.append(selected_counts)
            performance.append({"baseline": np.round(h["TA_mean"][()] * 100, 2).tolist(),
                                "selected": np.round(fs["accuracy_mean"][()] * 100, 2).tolist()})
        provenance.append(fingerprint(path))
    return np.stack(numerators), np.stack(denominators), np.stack(all_counts), labels, performance, provenance


def static_montage(values, denominator, atlas, anatomy, affine, model, label, lags, path, threshold=.8):
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "savefig.dpi": 160})
    fig, axes = plt.subplots(len(lags), len(CUTS), figsize=(12.4, 1.65 * len(lags) + 1),
                             squeeze=False, layout="constrained")
    cmap = plt.get_cmap("viridis").copy(); cmap.set_bad((0, 0, 0, 0))
    for row, ti in enumerate(lags):
        adequate = denominator[ti] >= 20
        freq = np.r_[np.nan, values[ti]]
        n = int(np.sum(values[ti] >= threshold)) if adequate else 0
        for col, z in enumerate(CUTS):
            zi = int(round((z - affine[2, 3]) / affine[2, 2]))
            plane = atlas[:, :, zi].T
            ax = axes[row, col]
            ax.imshow(anatomy[:, :, zi].T, origin="lower", cmap="gray", vmin=0, vmax=1)
            overlay = freq[plane]
            overlay[(overlay < threshold) | (plane == 0) | (not adequate)] = np.nan
            ax.imshow(overlay, origin="lower", cmap=cmap, vmin=0, vmax=1, alpha=.95)
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0: ax.set_title(f"MNI z = {z:+d} mm")
            if col == 0:
                ax.set_ylabel(f"{TIMES[ti]:.2f} s\n{n} ROIs; {int(denominator[ti])} reduced fits", fontsize=9)
                ax.text(.02, .05, "L", transform=ax.transAxes, color="white", fontsize=9)
                ax.text(.94, .05, "R", transform=ax.transAxes, color="white", fontsize=9)
            for spine in ax.spines.values(): spine.set_visible(False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=cmap), ax=axes.ravel().tolist(),
                       location="bottom", fraction=.025, pad=.01, shrink=.7)
    bar.set_label("Selection frequency among reduced-subset fits (only >=80% shown)")
    title = label.replace("_LoadDiff", "").replace("_vs_", " versus ")
    fig.suptitle(f"{model}: {title} load interaction\nExploratory stability display; >=20 reduced fits; neurological orientation", fontsize=12)
    fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def animation(values, denominator, atlas, anatomy, affine, path):
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 17)
    small = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 13)
    frames = []
    cmap = plt.get_cmap("viridis")
    for ti in range(41):
        frame = Image.new("RGB", (960, 580), "white")
        draw = ImageDraw.Draw(frame)
        count = int(np.sum(values[ti] >= .8)) if denominator[ti] >= 20 else 0
        draw.text((18, 8), f"cHRF | Body versus Face | lag {TIMES[ti]:.2f} s", fill="black", font=font)
        draw.text((18, 35), f"Selection >=80%; {count} ROIs across brain; {int(denominator[ti])} reduced fits", fill="black", font=small)
        lookup = np.r_[0, values[ti]]
        for col, z in enumerate(CUTS):
            zi = int(round((z-affine[2, 3])/affine[2, 2]))
            plane = atlas[:, :, zi].T
            gray = anatomy[:, :, zi].T
            rgb = np.repeat(gray[..., None], 3, axis=-1)
            active = (plane > 0) & (lookup[plane] >= .8) & (denominator[ti] >= 20)
            rgb[active] = cmap(lookup[plane[active]])[:, :3]
            img = Image.fromarray(np.uint8(np.clip(rgb[::-1], 0, 1) * 255)).resize((290, 234))
            x, y = 18 + (col % 3) * 315, 78 + (col // 3) * 250
            frame.paste(img, (x, y))
            draw.text((x + 5, y + 5), f"z={z:+d} mm", fill="white", font=small)
            draw.text((x + 5, y + 212), "L", fill="white", font=small)
            draw.text((x + 275, y + 212), "R", fill="white", font=small)
        frames.append(frame)
    # No loop extension: a single temporal pass rather than an endless animation.
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=400, optimize=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/svm_brain_maps_20261026")
    parser.add_argument("--template", type=Path, default=Path(__file__).with_name("predictive-brain-lags.template.html"))
    args = parser.parse_args(); out = args.out; out.mkdir(parents=True, exist_ok=True)
    atlas_dir = ROOT / "outputs/wm_scb_downstream_20261026/atlas"
    atlas_path = atlas_dir / "canlab2018_labels.nii"
    labels_path = atlas_dir / "atlas_regions.csv"
    template_path = Path("C:/Users/Ryan_/Documents/MATLAB/spm/canonical/avg152T1.nii")
    atlas_img = nib.as_closest_canonical(nib.load(atlas_path))
    atlas = np.asarray(atlas_img.dataobj).astype(np.uint16)
    assert np.array_equal(np.unique(atlas), np.arange(490))
    assert nib.aff2axcodes(atlas_img.affine) == ("R", "A", "S")
    assert np.allclose(atlas_img.affine[:3, :3], np.eye(3))
    region_labels = pd.read_csv(labels_path).sort_values("atlas_roi_index")
    assert region_labels.atlas_roi_index.tolist() == list(range(1, 490))
    t1 = resample_from_to(nib.load(template_path), atlas_img, order=1)
    anatomy = np.asarray(t1.dataobj, dtype=np.float32)
    anatomy = np.clip(anatomy / np.percentile(anatomy[anatomy > 0], 99.5), 0, 1)
    numerators, denominators, counts, labels, performance, sources = read_frequencies()
    partial = np.divide(numerators, denominators[..., None], out=np.full(numerators.shape, np.nan, dtype=np.float32),
                        where=denominators[..., None] > 0)
    region_labels.to_csv(out / "atlas_regions.csv", index=False)
    # Two-mm nearest-neighbor parcel preview preserves all 489 labels.
    preview = atlas[::2, ::2, ::2].astype("<u2")
    anatomy_preview = np.rint(anatomy[::2, ::2, ::2] * 31).astype(np.uint8)
    assert len(np.unique(preview)) == 490
    affine = atlas_img.affine.copy(); affine[:3, :3] *= 2
    payload = dict(models=MODELS, contrasts=labels, times=TIMES.round(2).tolist(),
        labels=region_labels.atlas_label.tolist(), networks=region_labels.atlas_network.tolist(),
        shape=list(preview.shape), affine=affine.tolist(), numerator=pack(numerators),
        denominator=pack(denominators), atlas=pack(preview), anatomy=pack(anatomy_preview), performance=performance)
    html = args.template.read_text(encoding="utf-8").replace("__PAYLOAD__", json.dumps(payload, separators=(",", ":")))
    assert "__PAYLOAD__" not in html
    destination = out / "predictive-brain-lags.html"
    destination.write_text(html, encoding="utf-8")
    assert destination.stat().st_size < 1_000_000, destination.stat().st_size
    rows = []
    for mi, model in enumerate(MODELS):
        for ci, label in enumerate(labels):
            for ti, sec in enumerate(TIMES):
                rows.append(dict(model=model, contrast=label, lag_index=ti+1, time_sec=sec,
                    reduced_fits=int(denominators[mi, ci, ti]),
                    stable_regions_80pct=int(np.sum(partial[mi, ci, ti] >= .8)) if denominators[mi, ci, ti] >= 20 else 0,
                    display_adequate=bool(denominators[mi, ci, ti] >= 20)))
            filename = out / f"{model}_{label}_brain_lags.png"
            static_montage(partial[mi, ci], denominators[mi, ci], atlas, anatomy, atlas_img.affine,
                           model, label, OVERVIEW_LAGS, filename)
        # Export primary contrast without thresholding: frequencies, not binary
        # significance masks. Inadequate-denominator lag volumes are NaN.
        lut = np.concatenate([np.zeros((1, 41), np.float32), partial[mi, 0].T.astype(np.float32)], axis=0)
        image4d = lut[preview]
        image4d[..., denominators[mi, 0] < 20] = np.nan
        nii = nib.Nifti1Image(image4d, affine)
        nii.header.set_xyzt_units("mm", "sec"); nii.header.set_zooms((2, 2, 2, .72))
        nii.header["descrip"] = b"SVM reduced-fit selection frequency; NOT significance; lag 0:0.72:28.8 s"
        nib.save(nii, out / f"{model}_Body_vs_Face_selection_frequency_4d.nii.gz")
        del image4d, nii
        print(f"Saved montages and primary 4D map: {model}", flush=True)
    for page, start in enumerate(range(0, 41, 7), 1):
        static_montage(partial[0, 0], denominators[0, 0], atlas, anatomy, atlas_img.affine,
            "cHRF", labels[0], list(range(start, min(start+7, 41))), out / f"cHRF_Body_vs_Face_all_lags_page_{page}.png")
    animation(partial[0, 0], denominators[0, 0], atlas, anatomy, atlas_img.affine, out / "cHRF_Body_vs_Face_all_lags.gif")
    pd.DataFrame(rows).to_csv(out / "display_counts_by_lag.csv", index=False)
    sources.extend(fingerprint(p) for p in [atlas_path, labels_path, template_path])
    metadata = dict(atlas_provenance=(atlas_dir / "atlas_provenance.txt").read_text(), sources=sources,
        dimensions=list(preview.shape), preview_resolution_mm=2, atlas_parcels_represented=489,
        neurological_orientation="left hemisphere on image left in axial and coronal views",
        default_threshold=.8, minimum_reduced_fits=20,
        frequency_definition="selected count among fits with fewer than 489 retained, divided by number of such fits",
        historical_atlas_identity="ordinal convention consistent with extraction code; historical version unverified",
        no_model_refitting=True, fragment_bytes=destination.stat().st_size)
    (out / "mapping_provenance.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    note = """# Predictive regional maps across lags

The maps show nested SVM-RFE selection frequency, not activation, coefficient
sign, p-values or a confidence-band result. Every model has 100 outer fits.
The default is frequency among reduced-subset fits only, with an exploratory
80% display threshold and at least 20 such fits. All-fit frequency is available
in the interactive viewer. Thresholding does not create statistical significance.

The viewer includes three HRF models, six load-interaction contrasts and all
41 lags from 0 to 28.80 seconds at 0.72-second steps. Move the lag slider or
play one temporal pass. Axial and coronal views use neurological orientation
(left hemisphere on the left). Slice sliders use MNI millimeters. Hover or tap
a parcel to inspect its label and current selection frequency. Gray tissue
is context; absent color can mean frequency below the display threshold.

Static overview montages show 0, 4.32, 8.64, 12.96, 17.28, 21.60, 25.92 and
28.80 seconds. There are six all-lag pages and a single-play GIF for cHRF
Body versus Face. These are displays of the same underlying data, not new
analyses chosen for a favorable result. Six fixed axial slices cannot display
every parcel; the viewer's three movable planes and 4D files support inspection.

The three primary-contrast NIfTI files contain unthresholded frequencies at
2 mm with the fourth dimension ordered by lag. Inadequate-denominator volumes
are NaN. Other atlas-background voxels are zero. Their time step is 0.72 seconds.

Mapping uses the local Canlab2018 export and its ordinal 1–489 region mapping,
matching the convention in Fit_HRF_Residuals_HCP.m. The exact atlas version
used in historical extraction cannot independently be established. Anatomical
names are therefore conditional on that ordinal mapping. The viewer and NIfTI
preview use a 2 mm nearest-neighbor atlas and retain all 489 labels. Static
montages use the native 1 mm atlas on the resampled SPM average MNI152 T1.
No anatomical or temporal smoothing is applied to selection frequencies.

Reproduce from the repository root:

```powershell
python code/residual_contrast/analysis/svm_downstream/plot_svm_brain_lags.py
```

Original results are unchanged. Input hashes and mapping assumptions are in
mapping_provenance.json. Display counts and reduced-fit denominators are in
display_counts_by_lag.csv. Atlas names and networks are in atlas_regions.csv.
"""
    (out / "README.md").write_text(note, encoding="utf-8")
    print(f"Finished: {destination} ({destination.stat().st_size:,} bytes)", flush=True)


if __name__ == "__main__":
    main()

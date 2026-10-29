"""Descriptive postprocessing of saved nested SVM-RFE results (no refitting).

Run from the repository root:
    python code/residual_contrast/analysis/svm_downstream/analyze_svm_results.py
HDF5 arrays use the reversed MATLAB dimension order, retained explicitly here:
bootstrap x contrast x lag x ROI. Raw condition arithmetic preserves single
precision to match the original MATLAB WM wrapper. No new p-values are made.
"""
from pathlib import Path
import argparse
import hashlib
import json
import warnings

import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[4]
MODELS = ["cHRF", "cHRFderiv", "sHRF"]
JOBS = {"cHRF": "job_36477504_task_0", "cHRFderiv": "job_36477598_task_1",
        "sHRF": "job_36477153_task_2"}
POP_JOBS = {"cHRF": "job_35821286", "cHRFderiv": "job_35821287", "sHRF": "job_35821288"}
COLORS = {"cHRF": "#146c94", "cHRFderiv": "#cc6b2c", "sHRF": "#6b4c9a"}


def text(h, obj):
    if isinstance(obj, h5py.Reference):
        obj = h[obj]
    return "".join(chr(int(v)) for v in obj[()].ravel(order="F") if v)


def mean(x, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(x, axis=axis)


def divide(a, b):
    a, b = np.broadcast_arrays(a, b)
    return np.divide(a, b, out=np.full(a.shape, np.nan), where=b != 0)


def rho(a, b):
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 3 or np.ptp(a[keep]) == 0 or np.ptp(b[keep]) == 0:
        return np.nan
    return float(spearmanr(a[keep], b[keep]).statistic)


def top_indices(x, n=50):
    # ROI index is the deterministic tie breaker; NaNs are last.
    return np.lexsort((np.arange(len(x)), -np.nan_to_num(x, nan=-np.inf)))[:n]


def longest_run(mask):
    edges = np.diff(np.r_[False, mask, False].astype(int))
    lengths = np.flatnonzero(edges == -1) - np.flatnonzero(edges == 1)
    return int(lengths.max()) if len(lengths) else 0


def save_csv(rows, out, name):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(out / name, index=False, float_format="%.10g")
    return frame


def source_record(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
            "sha256": digest.hexdigest()}


def read_model(model, input_root, out):
    folder = input_root / model / JOBS[model]
    path = folder / "observed_wm_interaction.mat"
    print(f"Reading {model}: {path}", flush=True)
    manifest = [source_record(path)]
    audit_rows = []
    with h5py.File(path, "r") as h:
        n, r, t, b = [int(h[k][()].item()) for k in ("S", "Rg", "T", "nBoot")]
        labels = [text(h, ref) for ref in h["contrast_labels"][()].ravel()]
        c = len(labels)
        assert (n, r, t, c, b) == (410, 489, 41, 6, 100), (model, n, r, t, c, b)
        fs = h["FeatureSelection"]
        assert bool(fs["enabled"][()].item())
        baseline = h["TA_obs"][()].astype(float)
        accuracy = fs["accuracy_obs"][()].astype(float)
        freq = fs["selection_frequency"][()].astype(float)
        baseline_reliance = h["Reliance_mean"][()].astype(float)
        selected = np.zeros((b, c, t, r), dtype=bool)
        eligible = np.zeros_like(selected)
        reliance = np.full(selected.shape, np.nan, dtype=np.float32)
        count = np.full((b, c, t), np.nan)
        inner_accuracy = np.full_like(count, np.nan)
        valid = np.zeros((b, c, t), dtype=bool)
        splits = []
        for i, ref in enumerate(fs["bootstrap"][()].ravel()):
            q = h[ref]
            selected[i] = q["selected"][()].astype(bool)
            eligible[i] = q["eligible"][()].astype(bool)
            reliance[i] = q["reliance"][()]
            count[i] = q["count"][()]
            inner_accuracy[i] = q["innerAccuracy"][()]
            valid[i] = q["valid"][()].astype(bool)
            split = h[h["split_list"][()].ravel()[i]]
            train = split["train_boot"][()].ravel().astype(int) - 1
            test = split["test"][()].ravel().astype(int) - 1
            unique_ids = q["uniqueSubjectIDs"][()].ravel().astype(int) - 1
            fold = q["subjectFold"][()].ravel().astype(int)
            assert len(np.intersect1d(train, test)) == 0
            assert np.array_equal(np.unique(train), unique_ids)
            assert len(unique_ids) == len(fold) and set(fold) == set(range(1, 6))
            assert np.all(selected[i] <= eligible[i])
            assert np.array_equal(selected[i].sum(axis=-1)[valid[i]], count[i][valid[i]])
            assert np.all(np.isnan(reliance[i][~selected[i]]))
            assert np.all(np.isfinite(reliance[i][selected[i]]))
            splits.append((train, test))
            audit_rows.append(dict(model=model, bootstrap=i + 1, training_rows=len(train),
                unique_training_subjects=len(unique_ids), oob_subjects=len(test),
                subject_overlap=0, inner_folds=len(np.unique(fold)), valid_cells=int(valid[i].sum())))
        assert np.allclose(freq, divide(selected.sum(axis=0), valid.sum(axis=0)[..., None]), equal_nan=True)
        assert np.allclose(accuracy, np.stack([h[z]["accuracy"][()] for z in fs["bootstrap"][()].ravel()]), equal_nan=True)
        assert np.allclose(fs["accuracy_difference_obs"][()], accuracy - baseline, equal_nan=True)
        config = {k: fs["config"][k][()].ravel().tolist() for k in fs["config"]}
        config["importance_repeats"] = int(h["meta/importance/importanceRepeats"][()].item())
        weights = h["W_obs"][()]
        pairdiff = h["PairDiff_obs"][()]
        oob_hits = h["TA_subj_hits"][()]
        # CSVs are independently joined by their explicit contrast/lag/ROI keys.
        for filename in ["selection_accuracy_wm_interaction.csv", "regional_selection_wm_interaction.csv",
                         "regional_reliance_wm_interaction.csv"]:
            manifest.append(source_record(folder / filename))
        acc_csv = pd.read_csv(folder / "selection_accuracy_wm_interaction.csv")
        sel_csv = pd.read_csv(folder / "regional_selection_wm_interaction.csv")
        rel_csv = pd.read_csv(folder / "regional_reliance_wm_interaction.csv")
        assert len(acc_csv) == c * t and len(sel_csv) == len(rel_csv) == c * t * r
        assert not acc_csv.duplicated(["contrast", "time_index"]).any()
        assert not sel_csv.duplicated(["contrast", "time_index", "roi_index"]).any()
        assert not rel_csv.duplicated(["contrast", "time_index", "roi_index"]).any()
        for ci, label in enumerate(labels):
            ac = acc_csv[acc_csv.contrast == label].sort_values("time_index")
            sc = sel_csv[sel_csv.contrast == label].sort_values(["time_index", "roi_index"])
            rc = rel_csv[rel_csv.contrast == label].sort_values(["time_index", "roi_index"])
            assert np.allclose(ac.baseline_accuracy, mean(baseline[:, ci], axis=0))
            assert np.allclose(ac.selected_accuracy, mean(accuracy[:, ci], axis=0))
            assert np.allclose(ac.paired_accuracy_difference, mean(accuracy[:, ci] - baseline[:, ci], axis=0))
            assert np.allclose(sc.selection_frequency.to_numpy().reshape(t, r), freq[ci])
            assert np.allclose(rc.accuracy_loss.to_numpy().reshape(t, r), baseline_reliance[ci], equal_nan=True)
            assert np.array_equal(ac.valid_boot, valid[:, ci].sum(axis=0))

    # Verify raw-data alignment numerically, rather than infer participant
    # identity from file-local row IDs. These checks do not establish HCP IDs.
    raw_path = ROOT / "data/task_residual" / f"WM{model}.mat"
    raw = loadmat(raw_path, variable_names=["Results"])["Results"]
    assert raw.shape == (r, t, n, 8)
    load_diffs = np.stack([(raw[:, :, :, j + 4] - raw[:, :, :, j]).transpose(2, 1, 0)
                           for j in range(4)])
    del raw
    score_errors = []
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for bi in [0, 49, 99]:
        tr, te = splits[bi]
        for ci, (a, z) in enumerate(pairs):
            for ti in [0, 10, 20, 30, 40]:
                A = load_diffs[a, :, ti].astype(float)
                B = load_diffs[z, :, ti].astype(float)
                X = np.concatenate([A[tr], B[tr]])
                mu = mean(X, axis=0)
                X = np.where(np.isnan(X), mu, X)
                sd = np.std(X, axis=0, ddof=1)
                keep = np.isfinite(mu) & (sd > 0)
                A = np.where(np.isnan(A), mu, A)
                B = np.where(np.isnan(B), mu, B)
                score = mean(((A[te][:, keep] - B[te][:, keep]) / sd[keep]) @ weights[bi, ci, ti, keep])
                score_errors.append(abs(score - pairdiff[bi, ci, ti]))
    assert max(score_errors) < 2e-5, (model, max(score_errors))
    manifest.append(source_record(raw_path))
    print(f"  {model}: score audit: {len(score_errors)} checks, max error={max(score_errors):.3g}", flush=True)
    del weights, pairdiff
    times = np.arange(t) * 0.72
    pop_folder = ROOT / "data/task_residual/jhpce_outputs_wm_population_omnibus_full" / labels[0] / model / POP_JOBS[model]
    pop_path = pop_folder / "temporal_omnibus_results.mat"
    pop_csv = pd.read_csv(pop_folder / "temporal_omnibus_roi_summary.csv")
    with h5py.File(pop_path, "r") as p:
        pop = p["results"]
        pop_mean = pop["mean_D"][()].astype(float)
        pop_sd = pop["sd_D"][()].astype(float)
        masks = pop["subjectMasks"][()].astype(bool)
        # roiIDs are MATLAB string objects (MCOS); the exported CSV carries
        # their explicit IDs, checked here without guessing MCOS internals.
        assert np.array_equal(pop_csv.roi_id.to_numpy(), np.arange(1, r + 1))
        assert np.allclose(pop["timeSec"][()].ravel(), times)
        assert masks.shape == (r, n)
        assert np.array_equal(masks.sum(axis=1), pop_csv.effective_n.to_numpy())
        # Omnibus uses whole-curve complete cases separately for each ROI;
        # the decoder instead imputes missing features using training means.
        raw_delta = load_diffs[0] - load_diffs[1]
        expected_mean = np.sum(np.where(masks.T[:, None, :], raw_delta, 0),
                               axis=0, dtype=np.float64) / masks.sum(axis=1)[None, :]
        pop_error = np.max(abs(expected_mean - pop_mean))
        assert pop_error < 1e-7, (model, pop_error)
        assert np.array_equal(pop["significant"][()].ravel().astype(bool), pop_csv.significant.astype(bool))
        assert np.allclose(pop["p_holm"][()].ravel(), pop_csv.p_holm)
        pop_has_scb = "simultaneousBands" in pop
        pop_nperm = int(pop["config/nPerm"][()].item())
    del load_diffs
    manifest.extend([source_record(pop_path), source_record(pop_folder / "temporal_omnibus_roi_summary.csv")])
    eliminated = valid & (count < eligible.sum(axis=-1))
    partial_count = eliminated.sum(axis=0)
    freq_partial = divide((selected & eliminated[..., None]).sum(axis=0), partial_count[..., None])
    selected_rel = divide(np.nansum(reliance, axis=0), np.isfinite(reliance).sum(axis=0))
    result = dict(model=model, n=n, r=r, t=t, b=b, labels=labels, times=times, baseline=baseline,
        accuracy=accuracy, frequency=freq, partial_frequency=freq_partial, partial_count=partial_count,
        selected=selected, eligible=eligible, reliance=reliance, baseline_reliance=baseline_reliance,
        selected_reliance=selected_rel, count=count, inner_accuracy=inner_accuracy, valid=valid,
        eliminated=eliminated, config=config, pop_csv=pop_csv, pop_mean=pop_mean, pop_sd=pop_sd,
        pop_has_scb=pop_has_scb, pop_nperm=pop_nperm,
        population_excluded_roi_count=int(np.sum(masks.sum(axis=1) < n)),
        population_effective_n_range=[int(masks.sum(axis=1).min()), int(masks.sum(axis=1).max())],
        max_score_error=max(score_errors), max_population_mean_error=float(pop_error),
        oob_hits_min=int(oob_hits.min()), oob_hits_max=int(oob_hits.max()),
        oob_counts=[len(te) for _, te in splits])
    return result, audit_rows, manifest


def summarize(results, out):
    summary_rows, temporal_rows, size_rows, regional_rows, top_rows = [], [], [], [], []
    pop_rows, pop_summary, stability_rows = [], [], []
    for x in results:
        model = x["model"]
        windows = {"all_lags_0_to_28.8s": np.ones(x["t"], dtype=bool),
                   "reference_4.32_to_8.64s": (x["times"] >= 4.32 - 1e-9) & (x["times"] <= 8.64 + 1e-9)}
        for ci, label in enumerate(x["labels"]):
            for ri in range(x["r"]):
                pf = x["partial_frequency"][ci, :, ri]
                adequate = x["partial_count"][ci] >= 20
                stable = (pf >= .8) & adequate
                peak = int(np.nanargmax(pf))
                stability_rows.append(dict(model=model, contrast=label, roi_index=ri + 1,
                    strong_lag_count=int(stable.sum()), longest_strong_run_lags=longest_run(stable),
                    adequate_reduced_fit_lags=int(adequate.sum()),
                    maximum_partial_frequency=float(pf[peak]), maximum_partial_frequency_sec=x["times"][peak],
                    reduced_fits_at_maximum=int(x["partial_count"][ci, peak]),
                    display_threshold=.8, minimum_reduced_fits=20))
            for ti, sec in enumerate(x["times"]):
                counts = x["count"][:, ci, ti]
                delta = x["accuracy"][:, ci, ti] - x["baseline"][:, ci, ti]
                temporal_rows.append(dict(model=model, contrast=label, time_index=ti + 1, time_sec=sec,
                    baseline_accuracy=mean(x["baseline"][:, ci, ti]), selected_accuracy=mean(x["accuracy"][:, ci, ti]),
                    difference_pp=100 * mean(delta), difference_resample_sd_pp=100 * np.nanstd(delta, ddof=1),
                    selected_count_median=np.nanmedian(counts), selected_count_q25=np.nanquantile(counts, .25),
                    selected_count_q75=np.nanquantile(counts, .75), selected_count_mean=mean(counts),
                    all_regions_fraction=1 - mean(x["eliminated"][:, ci, ti]),
                    valid_fits=int(x["valid"][:, ci, ti].sum()),
                    positive_baseline_reliance_regions=int(np.sum(x["baseline_reliance"][ci, ti] > 0)),
                    raw_selection_reliance_spearman=rho(x["frequency"][ci, ti], x["baseline_reliance"][ci, ti]),
                    partial_selection_reliance_spearman=rho(x["partial_frequency"][ci, ti], x["baseline_reliance"][ci, ti])))
            for window, mask in windows.items():
                delta = x["accuracy"][:, ci, mask] - x["baseline"][:, ci, mask]
                counts = x["count"][:, ci, mask]
                baseline_curve = mean(x["baseline"][:, ci], axis=0)
                selected_curve = mean(x["accuracy"][:, ci], axis=0)
                summary_rows.append(dict(model=model, contrast=label, window=window, lags=int(mask.sum()),
                    baseline_accuracy=mean(x["baseline"][:, ci, mask]), selected_accuracy=mean(x["accuracy"][:, ci, mask]),
                    difference_pp=100 * mean(delta), selected_count_median=np.nanmedian(counts),
                    selected_count_mean=mean(counts), all_regions_fraction=1 - mean(x["eliminated"][:, ci, mask]),
                    lag_cells_with_positive_difference=int(np.sum(mean(delta, axis=0) > 0)),
                    baseline_peak_accuracy=float(np.max(baseline_curve)), baseline_peak_sec=float(x["times"][np.argmax(baseline_curve)]),
                    selected_peak_accuracy=float(np.max(selected_curve)), selected_peak_sec=float(x["times"][np.argmax(selected_curve)])))
                for size in np.unique(counts[np.isfinite(counts)]):
                    size_rows.append(dict(model=model, contrast=label, window=window, regions_retained=int(size),
                        fits=int(np.sum(counts == size)), total_valid_fit_lag_cells=int(np.isfinite(counts).sum()),
                        fraction=float(np.mean(counts == size))))
                f = mean(x["frequency"][ci, mask], axis=0)
                pden = x["eliminated"][:, ci, mask].sum()
                pf = divide((x["selected"][:, ci, mask] & x["eliminated"][:, ci, mask, None]).sum(axis=(0, 1)), pden)
                br = mean(x["baseline_reliance"][ci, mask], axis=0)
                sr = divide(np.nansum(x["reliance"][:, ci, mask], axis=(0, 1)),
                            np.isfinite(x["reliance"][:, ci, mask]).sum(axis=(0, 1)))
                rank = np.empty(x["r"], int); rank[top_indices(pf, x["r"])] = np.arange(1, x["r"] + 1)
                frame = pd.DataFrame(dict(model=model, contrast=label, window=window,
                    roi_index=np.arange(1, x["r"] + 1), roi_name=[f"Parcel_{j:03}" for j in range(1, x["r"] + 1)],
                    mean_selection_frequency=f, partial_fit_selection_frequency=pf,
                    partial_fit_lag_cells=int(pden), mean_baseline_reliance_pp=100 * br,
                    selected_model_reliance_when_selected_pp=100 * sr, partial_frequency_rank=rank))
                regional_rows.extend(frame.to_dict("records"))
                top_rows.extend(frame[frame.partial_frequency_rank <= 20].to_dict("records"))
        # Matched mean-difference/omnibus table: only Body-versus-Face is available.
        pop = x["pop_csv"].sort_values("roi_id")
        assert list(pop.roi_id) == list(range(1, x["r"] + 1))
        sig = pop.significant.to_numpy(dtype=bool)
        for ti, sec in enumerate(x["times"]):
            for ri in range(x["r"]):
                pop_rows.append(dict(model=model, contrast=x["labels"][0], roi_index=ri + 1,
                    roi_name=pop.roi_name.iloc[ri], time_index=ti + 1, time_sec=sec,
                    mean_paired_difference=float(x["pop_mean"][ti, ri]),
                    standardized_mean_difference=float(divide(x["pop_mean"][ti, ri], x["pop_sd"][ti, ri])),
                    omnibus_S=float(pop.S_obs.iloc[ri]), omnibus_p_holm=float(pop.p_holm.iloc[ri]),
                    omnibus_significant=bool(sig[ri]), effective_n=int(pop.effective_n.iloc[ri]),
                    selection_frequency=float(x["frequency"][0, ti, ri]),
                    partial_fit_selection_frequency=float(x["partial_frequency"][0, ti, ri]),
                    baseline_reliance_pp=100 * float(x["baseline_reliance"][0, ti, ri]),
                    selected_model_reliance_when_selected_pp=100 * float(x["selected_reliance"][0, ti, ri]),
                    simultaneous_band_status="not_present_in_saved_mat"))
        f = mean(x["frequency"][0], axis=0)
        pf = divide((x["selected"][:, 0] & x["eliminated"][:, 0, :, None]).sum(axis=(0, 1)), x["eliminated"][:, 0].sum())
        br = mean(x["baseline_reliance"][0], axis=0)
        top50 = top_indices(pf)
        # Thresholds are exploratory display choices with sensitivity analysis.
        for threshold in [.5, .7, .8, .9]:
            stable = pf >= threshold
            pop_summary.append(dict(model=model, partial_frequency_threshold=threshold,
                planned_rois=x["r"], omnibus_significant_rois=int(sig.sum()),
                stable_rois=int(stable.sum()), significant_and_stable=int((sig & stable).sum()),
                significant_not_stable=int((sig & ~stable).sum()), not_significant_stable=int((~sig & stable).sum()),
                neither=int((~sig & ~stable).sum()), top50_partial_frequency_significant=int(sig[top50].sum()),
                partial_frequency_vs_omnibus_S_spearman=rho(pf, pop.S_obs.to_numpy()),
                raw_frequency_vs_omnibus_S_spearman=rho(f, pop.S_obs.to_numpy()),
                baseline_reliance_vs_omnibus_S_spearman=rho(br, pop.S_obs.to_numpy()),
                mean_partial_frequency_significant=mean(pf[sig]),
                mean_partial_frequency_not_significant=mean(pf[~sig])))
    summary = save_csv(summary_rows, out, "performance_summary.csv")
    temporal = save_csv(temporal_rows, out, "performance_by_lag.csv")
    save_csv(size_rows, out, "subset_size_distribution.csv")
    save_csv(stability_rows, out, "regional_temporal_stability.csv")
    regional = save_csv(regional_rows, out, "regional_summary.csv")
    save_csv(top_rows, out, "top20_regions_by_partial_selection.csv")
    matched = save_csv(pop_rows, out, "population_prediction_matched.csv")
    comparison = save_csv(pop_summary, out, "population_overlap_sensitivity.csv")
    cross_rows = []
    for ci, label in enumerate(results[0]["labels"]):
        for window in windows:
            vectors = {m: regional[(regional.model == m) & (regional.contrast == label) & (regional.window == window)]
                       .sort_values("roi_index") for m in MODELS}
            for i, a in enumerate(MODELS):
                for b in MODELS[i + 1:]:
                    A, B = vectors[a], vectors[b]
                    ai, bi = set(top_indices(A.partial_fit_selection_frequency.to_numpy())), set(top_indices(B.partial_fit_selection_frequency.to_numpy()))
                    cross_rows.append(dict(model_a=a, model_b=b, contrast=label, window=window,
                        raw_frequency_spearman=rho(A.mean_selection_frequency.to_numpy(), B.mean_selection_frequency.to_numpy()),
                        partial_frequency_spearman=rho(A.partial_fit_selection_frequency.to_numpy(), B.partial_fit_selection_frequency.to_numpy()),
                        baseline_reliance_spearman=rho(A.mean_baseline_reliance_pp.to_numpy(), B.mean_baseline_reliance_pp.to_numpy()),
                        top50_overlap=len(ai & bi), top50_jaccard=len(ai & bi) / len(ai | bi)))
    cross = save_csv(cross_rows, out, "cross_hrf_consistency.csv")
    return summary, temporal, regional, matched, comparison, cross


def plots(results, out):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 170, "axes.titleweight": "semibold"})
    labels = results[0]["labels"]
    titles = [z.replace("_LoadDiff", "").replace("_vs_", " vs ") for z in labels]
    fig, axes = plt.subplots(3, 6, figsize=(22, 10), sharex=True, sharey=True, constrained_layout=True)
    for row, x in enumerate(results):
        for ci in range(6):
            ax = axes[row, ci]
            ax.plot(x["times"], 100 * mean(x["baseline"][:, ci], axis=0), color="#777777", ls="--", label="All regions")
            ax.plot(x["times"], 100 * mean(x["accuracy"][:, ci], axis=0), color=COLORS[x["model"]], label="Nested RFE")
            ax.axhline(50, color="#bbbbbb", lw=.8)
            ax.set_ylim(45, 100); ax.set_xticks([0, 7.2, 14.4, 21.6, 28.8])
            if row == 0: ax.set_title(titles[ci])
            if ci == 0: ax.set_ylabel(x["model"] + "\nOOB accuracy (%)")
            if row == 2: ax.set_xlabel("Lag (seconds)")
    axes[0, 0].legend(fontsize=9, loc="lower right")
    fig.suptitle("Held-out condition decoding: all regions and nested feature selection\nMeans across 100 subject bootstraps; chance reference is descriptive, not a significance test")
    fig.savefig(out / "accuracy_comparison.png"); plt.close(fig)
    for kind in ["difference", "count"]:
        fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, constrained_layout=True)
        for ci, ax in enumerate(axes.ravel()):
            for x in results:
                if kind == "difference":
                    curve = 100 * mean(x["accuracy"][:, ci] - x["baseline"][:, ci], axis=0)
                    ax.plot(x["times"], curve, color=COLORS[x["model"]], label=x["model"])
                else:
                    curve = np.nanmedian(x["count"][:, ci], axis=0)
                    q = np.nanquantile(x["count"][:, ci], [.25, .75], axis=0)
                    ax.plot(x["times"], curve, color=COLORS[x["model"]], label=x["model"])
                    ax.fill_between(x["times"], q[0], q[1], color=COLORS[x["model"]], alpha=.08)
            if kind == "difference": ax.axhline(0, color="#777777", lw=.8)
            else: ax.set_ylim(0, 510)
            ax.set_title(titles[ci]); ax.set_xlabel("Lag (seconds)")
            ax.set_ylabel("Selected minus baseline (pp)" if kind == "difference" else "Regions retained (median and IQR)")
        axes[0, 0].legend(fontsize=9)
        fig.suptitle("Paired held-out accuracy difference" if kind == "difference" else "Selected subset size across bootstrap fits")
        fig.savefig(out / ("accuracy_difference.png" if kind == "difference" else "subset_size.png")); plt.close(fig)
    for x in results:
        for partial in [False, True]:
            fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
            for ci, ax in enumerate(axes.ravel()):
                values = x["partial_frequency"] if partial else x["frequency"]
                im = ax.imshow(values[ci].T, aspect="auto", origin="lower", vmin=0, vmax=1,
                    extent=[-.36, 29.16, .5, 489.5], cmap="viridis", interpolation="nearest")
                ax.set_title(titles[ci]); ax.set_xlabel("Lag (seconds)"); ax.set_ylabel("ROI index (atlas order)")
            fig.colorbar(im, ax=axes.ravel().tolist(), shrink=.8, label="Selection frequency")
            title = "Selection conditional on dropping regions" if partial else "Selection across all successful fits"
            fig.suptitle(f"{x['model']}: {title}")
            fig.savefig(out / f"selection_heatmap_{x['model']}{'_partial' if partial else ''}.png"); plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for x, ax in zip(results, axes):
        pf = divide((x["selected"][:, 0] & x["eliminated"][:, 0, :, None]).sum(axis=(0, 1)), x["eliminated"][:, 0].sum())
        p = x["pop_csv"].sort_values("roi_id")
        sig = p.significant.astype(bool).to_numpy()
        for status, color, name in [(False, "#a0a0a0", "Omnibus not significant"), (True, COLORS[x["model"]], "Omnibus significant")]:
            mask = sig == status
            ax.scatter(np.log10(p.S_obs.to_numpy()[mask]), pf[mask], c=color, s=14, alpha=.7, label=name)
        ax.set_title(x["model"]); ax.set_xlabel("log10 whole-curve omnibus statistic")
        ax.set_ylabel("Selection frequency in reduced-subset fits")
    axes[0].legend(fontsize=8)
    fig.suptitle("Body versus Face load interaction: population evidence and predictive selection\nEach point is one ROI; frequencies pooled across all lags and reduced-subset fits")
    fig.savefig(out / "population_vs_selection.png"); plt.close(fig)


def markdown_table(frame, columns, digits=3):
    # Avoid an optional dependency on tabulate.
    rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in frame.iterrows():
        vals = []
        for col in columns:
            v = row[col]
            vals.append(f"{v:.{digits}f}" if isinstance(v, (float, np.floating)) else str(v))
        rows.append("| " + " | ".join(vals) + " |")
    return "\n".join(rows)


def report(results, summary, temporal, regional, comparison, cross, out):
    all_lags = summary[summary.window == "all_lags_0_to_28.8s"]
    overview = []
    for x in results:
        overview.append(dict(model=x["model"], baseline_pct=100 * mean(x["baseline"]),
            selected_pct=100 * mean(x["accuracy"]), difference_pp=100 * mean(x["accuracy"] - x["baseline"]),
            median_regions=float(np.nanmedian(x["count"])), mean_regions=mean(x["count"]),
            all_regions_pct=100 * (1 - mean(x["eliminated"])),
            mean_regions_when_reduced=mean(x["count"][x["eliminated"]])))
    ov = save_csv(overview, out, "overall_model_summary.csv")
    stability = pd.read_csv(out / "regional_temporal_stability.csv")
    body_face = "Body_LoadDiff_vs_Face_LoadDiff"
    ranking = regional[(regional.contrast == body_face) & (regional.window == "all_lags_0_to_28.8s")]
    leader_rows = ranking.sort_values(["model", "partial_frequency_rank"]).groupby("model").head(5)
    narrative = [
        "Nested RFE does not improve average held-out accuracy: across all contrasts and lags, the mean change is −0.166 pp for cHRF, −0.172 pp for cHRFderiv and −0.341 pp for sHRF. The selected and all-region accuracy curves largely track each other. This is evidence of modest compression at a small average predictive cost, not evidence that the selected procedure is equivalent or superior. The only positive contrast-wide mean change is cHRF Body versus Tool (+0.181 pp); 17 of 18 model/contrast averages are negative.",
        "The median subset size is 200 in all three models. All 489 regions are retained in 42.5%, 42.5% and 39.5% of fit-lag-contrast cells, respectively. When elimination actually occurs, the mean subset is much smaller (see table). Thus the evidence does not support one tiny fixed regional set that replaces the full decoder at every lag and contrast.",
        "cHRF and cHRFderiv have higher descriptive average decoding accuracy than sHRF, and their regional selections agree more closely. Across six contrasts, their partial-frequency rank correlations range from 0.800 to 0.848 and their top-50 overlap is 33–38 regions. Comparisons involving sHRF have correlations about 0.396–0.487. This is sensitivity to HRF specification, not independent replication or a formal test of which HRF is best.",
        "For the matched Body-versus-Face interaction, population evidence and prediction align substantially but do not coincide. Of the top 50 regions by partial-fit selection frequency, 48/50 are omnibus significant for cHRF, 50/50 for cHRFderiv and 36/50 for sHRF. For context, the corresponding significant proportions over all 489 regions are 332/489, 304/489 and 114/489. Descriptive rank correlations between omnibus S and partial selection frequency are 0.642, 0.663 and 0.613. These observations support overlap between population differences and useful predictive regions; they do not establish unique or causal regional contributions.",
        "Selection is time-dependent. In the all-lag pooled Body-versus-Face ranking, no ROI reaches 80% selection among reduced-subset fits. That does not mean no ROI is stable at a particular lag: many ROI-lag cells exceed 80%. regional_temporal_stability.csv distinguishes isolated high-frequency cells from consecutive-lag persistence, requiring at least 20 reduced fits as an exploratory display guard.",
        "ROI 14 is the leading pooled Body-versus-Face candidate for cHRF and cHRFderiv. Its unconditional selection frequencies are 81.0% and 76.0%, but reduced-fit frequencies are 67.6% and 61.3%. Perturbing it in the baseline model reduces accuracy by an average 0.887 and 0.791 pp, respectively. It is a candidate for follow-up, not a uniquely identified anatomical mechanism. The sHRF ranking is different, led by ROI 382.",
        "Larger isolated accuracy gains exist (for example +5.01 pp for cHRF Body versus Tool at 21.60 s). These were identified after inspecting many lags and contrasts, so they are exploratory observations and should not replace the full-grid accuracy summaries or be used to claim an independently confirmed improvement."
    ]
    lines = ["# SVM feature-selection downstream analysis", "", "Analysis date: 2026-10-26.", "",
        "## Main findings", "", markdown_table(ov, list(ov.columns)), "",
        "Accuracy averages give each of the six contrasts and 41 sampled lags equal weight. They are descriptive summaries, not 73,800 independent observations. Differences are percentage points (pp).",
        "", "## Interpretation", "", "\n\n".join(narrative),
        "", "## Accuracy and subset size by contrast", "",
        markdown_table(all_lags, ["model", "contrast", "baseline_accuracy", "selected_accuracy", "difference_pp", "selected_count_median", "all_regions_fraction"]), "",
        "Feature selection is useful for compression only to the extent that fewer regions retain acceptable held-out accuracy. The current results estimate the performance of nested selection, not that of a single fixed regional subset. A small mean loss is not an equivalence or noninferiority result.", "",
        "Peak values in performance_summary.csv are selected after viewing these data and are exploratory. The primary summaries use every lag. A secondary 4.32–8.64 second reference window is carried forward from the repository README; it is not presented as a preregistered confirmatory analysis.", "",
        "## Selection stability", "",
        "High unconditional selection frequency can arise when the inner procedure keeps every region. We therefore also report selection conditional on a fit retaining fewer than its eligible regions. This is a descriptive diagnostic of selection under elimination, not a second validated feature-selection procedure and not a p-value. The conditioning changes the target population of fits.", "",
        "regional_summary.csv contains both frequencies and baseline-model reliance. Selected-model reliance is averaged only over fits where that region was selected. It must not be interpreted as reliance in every fit or added across regions. The perturbation preserves subject pairs but disrupts cross-region dependence.", "",
        "Leading Body-versus-Face candidates, pooled across all lags:", "",
        markdown_table(leader_rows, ["model", "roi_index", "mean_selection_frequency", "partial_fit_selection_frequency", "mean_baseline_reliance_pp"]), "",
        "Regions with at least one high-frequency Body-versus-Face lag (exploratory ≥80% partial frequency and at least 20 reduced fits):", "",
        markdown_table(stability[stability.contrast == body_face].groupby("model", sort=False).agg(
            regions_with_strong_lag=("strong_lag_count", lambda v: int((v > 0).sum())),
            regions_with_3_consecutive_lags=("longest_strong_run_lags", lambda v: int((v >= 3).sum())),
            maximum_consecutive_lags=("longest_strong_run_lags", "max")).reset_index(),
            ["model", "regions_with_strong_lag", "regions_with_3_consecutive_lags", "maximum_consecutive_lags"]), "",
        "## Cross-HRF consistency", "",
        markdown_table(cross[cross.window == "all_lags_0_to_28.8s"], ["model_a", "model_b", "contrast", "partial_frequency_spearman", "top50_overlap", "top50_jaccard"]), "",
        "Rank correlations and top-50 overlaps are descriptive. Top 50 is a reporting convention, not an inferential threshold. Correspondence uses the common 489-ROI order. Numeric raw-data checks support correspondence to each model's local raw file; real HCP participant IDs are unavailable, so shared participant identity across model files cannot be established from row numbers alone.", "",
        "## Matched population inference: Body versus Face load interaction", "",
        markdown_table(comparison[comparison.partial_frequency_threshold == .8], ["model", "omnibus_significant_rois", "stable_rois", "significant_and_stable", "not_significant_stable", "partial_frequency_vs_omnibus_S_spearman"]), "",
        "The 0.80 frequency cut is exploratory; the full sensitivity table also uses 0.50, 0.70 and 0.90. Population significance is the saved Holm-corrected whole-curve sign-flip result across 489 planned ROIs per model. Its independence and null-symmetry assumptions are inherited from the original run. It is not lag-specific significance. Population mean contrasts, standardized differences, ROI omnibus results, selection frequencies and model reliance are joined in population_prediction_matched.csv.", "",
        "No matched simultaneous confidence bands are present in these saved population result files. Category load-difference bands, and synthetic coverage-validation bands, cannot be substituted for this interaction. Band-based overlap remains unavailable without matched saved bands or an additional inferential run.", "",
        "## Verification and limits", "",
        "- All three full runs have 410 subjects, 489 ROIs, 41 lags, six contrasts and 100 outer bootstrap fits. The two-bootstrap cHRF pilot is excluded.",
        "- Source accuracy, selection and reliance CSVs were cross-checked against MAT arrays by explicit keys. Selected-mask counts, eligibility, saved accuracy differences and reliance missingness were checked.",
        "- Every saved outer split has disjoint training and held-out subject rows. Inner-fold subject IDs match unique outer-training rows. This checks subject separation, not HCP family independence.",
        "- Raw WM condition slots were reconstructed with single-precision arithmetic. Across each HRF model, 90 saved paired-score summaries were reproduced from raw features, training-only scales and saved weights (within floating-point tolerance). The raw Body-versus-Face population mean was also matched across all ROI-lag cells using the saved ROI-specific whole-curve complete-case masks.",
        "- Missing-data handling differs: the SVM imputes missing features from training means, whereas the population test excludes subjects with incomplete curves for each ROI. The matched table retains effective N, and source_manifest.json records affected ROI counts and effective-N ranges. These differences matter when interpreting individual regions.",
        "- Time is (index−1) × 0.72 seconds, numerically matched to the saved population time vector.",
        "- Generic parcel labels are used. No verified local ROI-to-network/anatomical label map was available, so network or anatomical interpretation is not invented.",
        "- No new p-values, confidence intervals for selection, above-chance tests, or tests treating bootstrap fits as independent were computed. Comparing methods on these same data is not independent replication.", "",
        "## Reproduce", "", "```powershell", "python code/residual_contrast/analysis/svm_downstream/analyze_svm_results.py", "```", "",
        "Input paths and SHA-256 hashes are recorded in source_manifest.json. Original inputs are unchanged. Figures show all contrasts and lags; no top-region or peak choice was used to refit a decoder."]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return ov


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=ROOT / "data/feature_selection/feature_selection")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/svm_downstream_20261026")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results, audits, manifest = [], [], []
    for model in MODELS:
        x, a, sources = read_model(model, args.input_root, args.out)
        results.append(x); audits.extend(a); manifest.extend(sources)
    save_csv(audits, args.out, "split_audit.csv")
    summary, temporal, regional, matched, comparison, cross = summarize(results, args.out)
    plots(results, args.out)
    ov = report(results, summary, temporal, regional, comparison, cross, args.out)
    metadata = {"analysis_date": "2026-10-26", "sources": manifest,
        "checks": [{"model": x["model"], "config": x["config"],
            "score_checks": 90, "max_paired_score_error": x["max_score_error"],
            "max_raw_population_mean_error": x["max_population_mean_error"],
            "oob_subject_hits_range": [x["oob_hits_min"], x["oob_hits_max"]],
            "oob_size_range": [min(x["oob_counts"]), max(x["oob_counts"])],
            "population_permutations": x["pop_nperm"], "matched_bands_present": x["pop_has_scb"],
            "population_excluded_roi_count": x["population_excluded_roi_count"],
            "population_effective_n_range": x["population_effective_n_range"]}
            for x in results]}
    # Encode the all-features candidate explicitly rather than nonstandard JSON Infinity.
    for item in metadata["checks"]:
        item["config"]["featureCounts"] = ["all" if not np.isfinite(v) else int(v)
                                                for v in item["config"]["featureCounts"]]
    (args.out / "source_manifest.json").write_text(json.dumps(metadata, indent=2, allow_nan=False), encoding="utf-8")
    print(ov.to_string(index=False), flush=True)
    print("Population:", comparison[comparison.partial_frequency_threshold == .8].to_string(index=False), flush=True)
    print(f"Saved outputs to {args.out}", flush=True)


if __name__ == "__main__":
    main()

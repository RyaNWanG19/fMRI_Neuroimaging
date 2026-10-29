"""Audit and summarize saved WM simultaneous bands; never rerun primary tests.

Python 3: numpy, pandas, scipy, h5py, matplotlib, scikit-learn, nibabel.
The atlas exporter is a separate MATLAB step. All output is derived analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from threadpoolctl import threadpool_limits

MODELS = ["cHRF", "cHRFderiv", "sHRF"]
CATEGORIES = ["body", "faces", "places", "tools"]
JOBS = {"cHRF": "36689156", "cHRFderiv": "36689158", "sHRF": "36689155"}
COLORS = ["#1675a9", "#d47726", "#6a4aa2"]
NROI, NT = 489, 41
TIMES = np.arange(NT) * .72
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 170})


def decode(ds):
    return "".join(chr(int(x)) for x in ds[()].ravel())


def scalar(group, key):
    return group[key][()].item()


def holm(p):
    ix = np.argsort(p, kind="stable")
    out = np.empty_like(p)
    out[ix] = np.minimum(1, np.maximum.accumulate(p[ix] * np.arange(len(p), 0, -1)))
    return out


def save_csv(frame, path):
    frame.to_csv(path, index=False, float_format="%.12g")


def finish(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def blocks(mask):
    edges = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1) - 1))


def audit_model(root, output, atlas, model):
    job = root / "data/AllFourCategory_LoadDiff/AllFourCategory_LoadDiff" / model / ("job_" + JOBS[model])
    band_csv = job / "simultaneous_confidence_bands.csv"
    mat_file = job / "temporal_omnibus_results.mat"
    bands = pd.read_csv(band_csv).sort_values(["test_index", "lag_index"]).reset_index(drop=True)
    summary = pd.read_csv(job / "temporal_omnibus_roi_summary.csv")
    assert len(bands) == NROI * 4 * NT and len(summary) == NROI * 4
    assert np.array_equal(bands.test_index, np.repeat(np.arange(1, 1957), NT))
    assert np.array_equal(bands.lag_index, np.tile(np.arange(1, 42), 1956))
    assert np.allclose(bands.time_sec, np.tile(TIMES, 1956))
    assert np.array_equal(summary.atlas_roi_index, np.tile(np.arange(1, 490), 4))
    assert np.array_equal(summary.category, np.repeat(CATEGORIES, NROI))
    assert np.array_equal(bands.category, np.repeat(summary.category, NT))
    assert np.array_equal(bands.atlas_roi_index, np.repeat(summary.atlas_roi_index, NT))
    assert not bands.duplicated(["test_index", "lag_index"]).any()
    matrices = {col: bands[col].to_numpy().reshape(1956, NT) for col in [
        "mean_difference", "standard_error", "lower_simultaneous", "upper_simultaneous",
        "band_eligible", "supported_positive", "supported_negative"]}
    with h5py.File(mat_file, "r") as f:
        r, cfg = f["results"], f["results/config"]
        b = r["simultaneousBands"]
        q = scalar(b, "criticalValue")
        audit = {"model": model, "job": JOBS[model], "source_mat": str(mat_file.relative_to(root)),
                 "mat_sha256": hashlib.sha256(mat_file.read_bytes()).hexdigest(),
                 "source_band_csv": str(band_csv.relative_to(root)), "n_permutations": int(scalar(cfg, "nPerm")),
                 "run_mode": decode(cfg["runMode"]), "alpha": scalar(cfg, "alpha"),
                 "bootstrap_method": decode(b["method"]),
                 "bootstrap_replicates": int(scalar(b["config"], "nBootstrap")), "critical_value": q,
                 "degenerate_bootstrap_replicates": int(scalar(b, "degenerateBootstrapReplicates")),
                 "planned_cells": int(scalar(b, "plannedCellCount")),
                 "calibrated_cells": int(scalar(b, "calibratedCellCount")),
                 "inference_valid": bool(scalar(r, "inferenceValid")),
                 "band_inference_valid": bool(scalar(b, "inferenceValid")),
                 "band_inference_status": " | ".join(sorted(bands.inference_status.unique())),
                 "independence_assumed": bool(scalar(cfg["resampling"], "independentSubjectsAssumed")),
                 "independence_verified_flag": bool(scalar(cfg["resampling"], "independentSubjectsVerified")),
                 "null_symmetry_assumed": bool(scalar(cfg["resampling"], "nullSymmetryAssumed")),
                 "subject_count": r["subjectMasks"].shape[1],
                 "effective_n_min": int(r["effectiveN"][()].min()),
                 "effective_n_max": int(r["effectiveN"][()].max())}
        assert audit["n_permutations"] == 99999 and audit["bootstrap_replicates"] == 9999
        assert audit["run_mode"] == "inference" and audit["inference_valid"] and audit["band_inference_valid"]
        assert audit["bootstrap_method"] == "subject-bootstrap-t" and np.isfinite(q)
        assert audit["planned_cells"] == 80196 and audit["alpha"] == .05
        maxima = b["bootstrapMaxima"][()].ravel()
        rank = int(np.ceil((len(maxima) + 1) * .95)) - 1
        assert np.sort(maxima)[rank] == q and np.isfinite(maxima).all()
        # Order-statistic sensitivity, not a confidence interval for actual coverage.
        tail_se = np.sqrt(len(maxima) * .95 * .05)
        audit["critical_order_sensitivity_low"] = float(np.sort(maxima)[max(0, int(rank - 2 * tail_se))])
        audit["critical_order_sensitivity_high"] = float(np.sort(maxima)[min(len(maxima)-1, int(rank + 2 * tail_se))])
        mu = r["mean_D"][()].T
        sd = r["sd_D"][()].T
        tobs = r["T_obs"][()].T
        mask = r["subjectMasks"][()].T.astype(bool)
        n = r["effectiveN"][()].ravel()
        assert np.array_equal(n, mask.sum(axis=0))
        assert np.allclose(summary.p_holm, holm(summary.p_raw.to_numpy()), atol=1e-11)
        assert np.array_equal(summary.significant.to_numpy(bool), summary.p_holm <= .05)
        for key, field in [("mean_difference", "mean"), ("standard_error", "standardError"),
                           ("lower_simultaneous", "lower"), ("upper_simultaneous", "upper"),
                           ("band_eligible", "eligible"), ("supported_positive", "positive"),
                           ("supported_negative", "negative")]:
            assert np.allclose(matrices[key], b[field][()].T, equal_nan=True)
        assert np.allclose(matrices["standard_error"], sd / np.sqrt(n[:, None]))
        eligible = matrices["band_eligible"].astype(bool)
        assert np.array_equal(matrices["supported_positive"].astype(bool), eligible & (matrices["lower_simultaneous"] > 0))
        assert np.array_equal(matrices["supported_negative"].astype(bool), eligible & (matrices["upper_simultaneous"] < 0))
        assert np.allclose(matrices["lower_simultaneous"][eligible], (mu-q*sd/np.sqrt(n[:,None]))[eligible])
        assert np.allclose(matrices["upper_simultaneous"][eligible], (mu+q*sd/np.sqrt(n[:,None]))[eligible])
        assert np.array_equal(bands.effective_n, np.repeat(n, NT))

    raw_path = root / "data/task_residual" / ("WM" + model + ".mat")
    raw = loadmat(raw_path, variable_names=["Results"])["Results"]
    assert raw.shape == (NROI, NT, audit["subject_count"], 8)
    ns = raw.shape[2]
    D = np.empty((ns, 4*NROI, NT), dtype=float)
    zero_conditions = []
    for c in range(4):
        D[:, c*NROI:(c+1)*NROI, :] = (raw[:, :, :, c+4].astype(float)-raw[:, :, :, c].astype(float)).transpose(2,0,1)
    for c in range(8):
        zero_conditions.append(int(np.all(raw[:,:,:,c] == 0, axis=(0,1)).sum()))
    del raw
    observed_mask = np.isfinite(D).all(axis=2)
    assert np.array_equal(mask, observed_mask), "Subject masks differ from current local source."
    D[~observed_mask] = 0
    direct_mean = D.sum(axis=0) / n[:,None]
    centered = np.where(observed_mask[:,:,None], D-direct_mean, 0)
    direct_sd = np.sqrt((centered**2).sum(axis=0)/(n[:,None]-1))
    direct_t = direct_mean/(direct_sd/np.sqrt(n[:,None]))
    audit.update({"source_raw": str(raw_path.relative_to(root)), "source_raw_bytes": raw_path.stat().st_size,
                  "max_source_mean_error": float(np.abs(mu-direct_mean).max()),
                  "max_source_sd_error": float(np.abs(sd-direct_sd).max()),
                  "max_source_t_error": float(np.abs(tobs-direct_t).max()),
                  "zero_filled_whole_condition_subject_counts": zero_conditions,
                  "source_has_extraction_flags_or_real_subject_ids": False,
                  "excluded_subject_test_pairs": int((~mask).sum())})
    assert audit["max_source_mean_error"] < 1e-10
    assert audit["max_source_sd_error"] < 1e-10
    assert audit["max_source_t_error"] < 1e-8
    assert np.allclose((direct_t**2).sum(axis=1), summary.S_obs, rtol=1e-10)
    del centered

    pos = matrices["supported_positive"].astype(bool)
    neg = matrices["supported_negative"].astype(bool)
    signed = pos.astype(np.int8) - neg.astype(np.int8)
    any_support = pos | neg
    evidence = summary.merge(atlas, on="atlas_roi_index", validate="many_to_one")
    evidence["model"] = model
    evidence["band_positive_lags"] = pos.sum(axis=1)
    evidence["band_negative_lags"] = neg.sum(axis=1)
    evidence["band_supported_lags"] = any_support.sum(axis=1)
    evidence["band_any"] = any_support.any(axis=1)
    evidence["band_biphasic"] = pos.any(axis=1) & neg.any(axis=1)
    evidence["mean_curve_rms"] = np.sqrt(np.mean(mu**2, axis=1))
    peak = np.argmax(np.abs(mu), axis=1)
    evidence["peak_abs_mean"] = np.abs(mu[np.arange(1956),peak])
    evidence["peak_signed_mean"] = mu[np.arange(1956),peak]
    evidence["peak_lag_sec"] = TIMES[peak]
    evidence["first_supported_lag_sec"] = [TIMES[np.flatnonzero(x)[0]] if x.any() else np.nan for x in any_support]
    evidence["last_supported_lag_sec"] = [TIMES[np.flatnonzero(x)[-1]] if x.any() else np.nan for x in any_support]
    evidence["evidence_class"] = np.select([
        evidence.significant.astype(bool) & evidence.band_any,
        evidence.significant.astype(bool), evidence.band_any], ["omnibus_and_bands", "omnibus_only", "bands_only"], default="neither")
    # Guaranteed signed magnitude at sampled lags; descriptive ranking of
    # simultaneous bounds, not an independent significance test.
    lower, upper = matrices["lower_simultaneous"], matrices["upper_simultaneous"]
    evidence["max_band_distance_from_zero"] = np.maximum(np.maximum(lower, -upper),0).max(axis=1)
    audit["band_cells_unstable_under_order_sensitivity"] = int(np.count_nonzero(
        (np.abs(tobs)>audit["critical_order_sensitivity_low"]) !=
        (np.abs(tobs)>audit["critical_order_sensitivity_high"])))
    audit["all_source_checks_passed"] = True
    segments = []
    for i in range(1956):
        for direction, mask_row in [("positive", pos[i]), ("negative", neg[i])]:
            for first, last in blocks(mask_row):
                segments.append({"model": model, "category": evidence.category.iloc[i],
                                 "atlas_roi_index": int(evidence.atlas_roi_index.iloc[i]),
                                 "atlas_label": evidence.atlas_label.iloc[i], "direction": direction,
                                 "first_sample_sec": TIMES[first], "last_sample_sec": TIMES[last],
                                 "supported_sample_count": last-first+1,
                                 "mean_effect_over_samples": float(mu[i,first:last+1].mean()),
                                 "boundary_interpretation": "sampled support run; not exact onset/duration"})
    print(f"AUDIT PASS {model}: q={q:.6f}, N={n.min()}-{n.max()}, mean error={audit['max_source_mean_error']:.2g}", flush=True)
    return {"model":model, "job":job, "audit":audit, "bands":bands, "summary":summary,
            "matrices":matrices, "mean":mu, "sd":sd, "n":n, "mask":mask,
            "D":D, "signed":signed, "evidence":evidence, "segments":segments}


def temporal_plots(result, figdir):
    model = result["model"]
    for c, category in enumerate(CATEGORIES):
        ii = slice(c*NROI,(c+1)*NROI)
        mu = result["mean"][ii]
        signed = result["signed"][ii]
        vmax = float(np.quantile(np.abs(mu), .99))
        fig, axes = plt.subplots(1,2,figsize=(12,7),constrained_layout=True,sharey=True)
        im = axes[0].imshow(mu,aspect="auto",origin="upper",interpolation="nearest",cmap="RdBu_r",vmin=-vmax,vmax=vmax,
                            extent=[-.36,29.16,489.5,.5])
        plt.colorbar(im,ax=axes[0],label="Mean 2-back minus 0-back residual (99th-percentile color limit)",shrink=.8)
        axes[0].set_title("Estimated mean difference")
        im=axes[1].imshow(signed,aspect="auto",origin="upper",interpolation="nearest",cmap=ListedColormap(["#2465a5","#f1f1f1","#c83832"]),
                          vmin=-1,vmax=1,extent=[-.36,29.16,489.5,.5])
        cb=plt.colorbar(im,ax=axes[1],ticks=[-1,0,1],shrink=.8)
        cb.ax.set_yticklabels(["Negative support","Band includes zero","Positive support"])
        axes[1].set_title("Joint simultaneous-band decisions")
        for ax in axes:
            ax.set_xlabel("Lag since onset (s)"); ax.set_xticks([0,7.2,14.4,21.6,28.8])
        axes[0].set_ylabel("Atlas ROI index (original order)")
        fig.suptitle(f"{category.capitalize()} | {model} | correction over 80,196 sampled cells within this model")
        finish(fig,figdir/f"heatmap_{model}_{category}.png")


def functional_clusters(result, outdir, figdir, bootstraps):
    """Exploratory ROI clustering; all categories concatenated, not selected ROIs.

    Weighted PCA is the grid approximation to FPCA. All four categories are
    retained in a common subject resample. PCA and clustering are refitted in
    every replicate at fixed selected K; K selection is separately reported.
    """
    model, D, mask = result["model"], result["D"], result["mask"]
    ns = D.shape[0]
    rng = np.random.default_rng(20261026)
    weights = rng.multinomial(ns, np.full(ns,1/ns), size=bootstraps)
    counts = weights @ mask.astype(float)
    boot_means = (weights @ D.reshape(ns,-1)).reshape(bootstraps,1956,NT)/counts[:,:,None]
    boot_means = boot_means.reshape(bootstraps,4,NROI,NT).transpose(0,2,1,3).reshape(bootstraps,NROI,4*NT)
    mean = result["mean"].reshape(4,NROI,NT).transpose(1,0,2).reshape(NROI,4*NT)
    time_weights = np.full(NT,.72); time_weights[[0,-1]] *= .5
    basis_weight = np.sqrt(np.tile(time_weights,4))
    selected_summaries, assignments, candidates = [], [], []
    result["clusters"]={}
    for mode in ["amplitude", "shape"]:
        def transform(x):
            y=x*basis_weight
            if mode=="shape":
                y=y/np.maximum(np.linalg.norm(y,axis=1,keepdims=True),1e-12)
            return y
        X=transform(mean)
        pca=PCA(n_components=.95,svd_solver="full")
        scores=pca.fit_transform(X)
        models={}
        for k in range(2,7):
            km=KMeans(n_clusters=k,n_init=30,random_state=8301).fit(scores)
            sizes=np.bincount(km.labels_,minlength=k)
            score=float(silhouette_score(scores,km.labels_))
            candidates.append({"model":model,"mode":mode,"k":k,"silhouette":score,
                               "smallest_cluster":int(sizes.min()),"fpca_modes":scores.shape[1]})
            models[k]=(score,km)
        k=max(models,key=lambda kk:models[kk][0])
        silhouette, km=models[k]
        labels=km.labels_
        stability=np.zeros(NROI); coassignment=np.zeros((NROI,NROI)); aris=[]
        for b in range(bootstraps):
            xb=transform(boot_means[b])
            sb=PCA(n_components=.95,svd_solver="full").fit_transform(xb)
            lb=KMeans(n_clusters=k,n_init=10,random_state=8301+b).fit_predict(sb)
            contingency=np.zeros((k,k),int)
            np.add.at(contingency,(labels,lb),1)
            rr,cc=linear_sum_assignment(-contingency)
            mapping=np.zeros(k,int); mapping[cc]=rr
            aligned=mapping[lb]
            stability += aligned==labels
            coassignment += lb[:,None]==lb[None,:]
            aris.append(adjusted_rand_score(labels,lb))
        stability /= bootstraps; coassignment /= bootstraps
        stats={"model":model,"mode":mode,"selected_k":k,"fpca_modes":scores.shape[1],
               "variance_explained":float(pca.explained_variance_ratio_.sum()),
               "silhouette":silhouette,"bootstrap_replicates":bootstraps,
               "median_bootstrap_ARI":float(np.median(aris)),
               "bootstrap_ARI_025":float(np.quantile(aris,.025)),
               "bootstrap_ARI_975":float(np.quantile(aris,.975)),
               "median_ROI_assignment_stability":float(np.median(stability)),
               "ROIs_with_assignment_stability_ge_080":int((stability>=.8).sum())}
        selected_summaries.append(stats)
        mapping=result["evidence"].drop_duplicates("atlas_roi_index").sort_values("atlas_roi_index")
        for roi in range(NROI):
            assignments.append({"model":model,"mode":mode,"atlas_roi_index":roi+1,
                                "atlas_label":mapping.atlas_label.iloc[roi],
                                "atlas_network":mapping.atlas_network.iloc[roi],
                                "cluster":int(labels[roi]+1),"assignment_stability":stability[roi]})
        np.savez_compressed(outdir/f"cluster_{model}_{mode}.npz",labels=labels+1,
                            fpca_scores=scores,fpca_loadings=pca.components_,
                            explained_variance_ratio=pca.explained_variance_ratio_,
                            subject_bootstrap_ari=aris,coassignment=coassignment)
        fig,axs=plt.subplots(k,4,figsize=(14,2.3*k),squeeze=False,sharex=True)
        limits=float(np.quantile(np.abs(mean),.995))
        for g in range(k):
            ids=labels==g
            for c,cat in enumerate(CATEGORIES):
                ax=axs[g,c]
                curves=mean[ids,c*NT:(c+1)*NT]
                if mode=="shape":
                    curves=mean[ids,c*NT:(c+1)*NT]/np.maximum(np.linalg.norm(mean[ids]*basis_weight,axis=1,keepdims=True),1e-12)
                lo,hi=np.quantile(curves,[.1,.9],axis=0)
                ax.fill_between(TIMES,lo,hi,color="#aac7d7",alpha=.5)
                ax.plot(TIMES,np.median(curves,axis=0),color="#18465d",lw=1.6)
                ax.axhline(0,color="grey",ls=":",lw=.7)
                if g==0: ax.set_title(cat)
                if c==0: ax.set_ylabel(f"Group {g+1}, n={ids.sum()}\nMedian stability={np.median(stability[ids]):.2f}")
                if g==k-1: ax.set_xlabel("Lag (s)")
                if mode=="amplitude": ax.set_ylim(-limits,limits)
        fig.suptitle(f"Exploratory {mode} groups | {model} | K={k}, silhouette={silhouette:.2f}\n"
                     "Medians and 10–90% between-region ranges; shading is NOT a confidence band",fontsize=12)
        fig.tight_layout(rect=(0,0,1,.94))
        finish(fig,figdir/f"clusters_{model}_{mode}.png")
        result["clusters"][mode]={"labels":labels,"stability":stability,"stats":stats}
        print(f"CLUSTERS {model}/{mode}: K={k}, silhouette={silhouette:.3f}, median ARI={np.median(aris):.3f}",flush=True)
    return selected_summaries, assignments, candidates


def model_comparisons(results, outdir, figdir):
    rows, overlaps, consensus = [], [], []
    means=np.stack([r["mean"] for r in results])
    signs=np.stack([r["signed"] for r in results])
    evidence=[r["evidence"] for r in results]
    for c,cat in enumerate(CATEGORIES):
        ix=slice(c*NROI,(c+1)*NROI)
        support=signs[:,ix]!=0
        samepos=(signs[:,ix]==1).all(axis=0)
        sameneg=(signs[:,ix]==-1).all(axis=0)
        opposite=(signs[:,ix].max(axis=0)==1)&(signs[:,ix].min(axis=0)==-1)
        for roi in range(NROI):
            consensus.append({"category":cat,"atlas_roi_index":roi+1,
                "atlas_label":evidence[0].atlas_label.iloc[c*NROI+roi],
                "atlas_network":evidence[0].atlas_network.iloc[c*NROI+roi],
                "models_with_any_supported_lag":int(support[:,:,].any(axis=2)[:,roi].sum()),
                "same_positive_lags_all_models":int(samepos[roi].sum()),
                "same_negative_lags_all_models":int(sameneg[roi].sum()),
                "opposite_supported_direction_lags":int(opposite[roi].sum())})
        for a,b in [(0,1),(0,2),(1,2)]:
            aa,bb=support[a],support[b]
            union=(aa|bb).sum(); intersection=(aa&bb).sum()
            both=aa&bb
            agree=both&(signs[a,ix]==signs[b,ix]); conflict=both&~agree
            x,y=means[a,ix],means[b,ix]
            cors=np.array([np.corrcoef(xx,yy)[0,1] for xx,yy in zip(x,y)])
            rows.append({"category":cat,"model_a":MODELS[a],"model_b":MODELS[b],
                         "support_cell_jaccard":intersection/union if union else np.nan,
                         "jointly_supported_cells":int(both.sum()),"same_direction_cells":int(agree.sum()),
                         "opposite_direction_cells":int(conflict.sum()),
                         "opposite_fraction_among_jointly_supported":conflict.sum()/both.sum() if both.any() else np.nan,
                         "median_ROI_curve_correlation":float(np.nanmedian(cors)),
                         "mean_curve_rms_ratio_b_over_a":float(np.median(np.sqrt((y*y).mean(1))/np.sqrt((x*x).mean(1))))})
        flags=np.stack([ev.band_any.to_numpy()[ix] for ev in evidence],axis=1)
        for code in range(8):
            pattern=np.array([(code>>bit)&1 for bit in range(3)],bool)
            overlaps.append({"category":cat,"cHRF":int(pattern[0]),"cHRFderiv":int(pattern[1]),
                             "sHRF":int(pattern[2]),"regions":int(np.all(flags==pattern,axis=1).sum())})
    save_csv(pd.DataFrame(rows),outdir/"hrf_pairwise_comparison.csv")
    save_csv(pd.DataFrame(overlaps),outdir/"hrf_region_overlap_patterns.csv")
    consensus=pd.DataFrame(consensus)
    save_csv(consensus,outdir/"hrf_direction_consensus.csv")

    # Deliberate post-result examples: first unused ROI by atlas index in
    # each evidence pattern. Patterns can overlap across different lags.
    examples=[]
    fig,axs=plt.subplots(4,3,figsize=(14,12),sharex=True)
    for c,cat in enumerate(CATEGORIES):
        cc=consensus[consensus.category==cat]
        used=set()
        patterns=[("Opposite supported directions",cc[cc.opposite_supported_direction_lags>0]),
                  ("Same-direction support in all models",cc[(cc.same_positive_lags_all_models+cc.same_negative_lags_all_models)>0]),
                  ("No supported lag in any model",cc[cc.models_with_any_supported_lag==0])]
        for col,(name,subset) in enumerate(patterns):
            ax=axs[c,col]
            subset=subset[~subset.atlas_roi_index.isin(used)]
            if subset.empty:
                ax.text(.5,.5,"No region in this pattern",ha="center",transform=ax.transAxes)
                ax.set_title(f"{cat}: {name}",fontsize=9); continue
            roi=int(subset.atlas_roi_index.iloc[0]); ii=c*NROI+roi-1
            used.add(roi)
            label=evidence[0].atlas_label.iloc[ii]
            for j,r in enumerate(results):
                mat=r["matrices"]
                ax.plot(TIMES,r["mean"][ii],color=COLORS[j],label=MODELS[j],lw=1.2)
                ax.fill_between(TIMES,mat["lower_simultaneous"][ii],mat["upper_simultaneous"][ii],color=COLORS[j],alpha=.12)
                sig=r["signed"][ii]!=0
                ax.scatter(TIMES[sig],r["mean"][ii,sig],color=COLORS[j],s=7)
            examples.append({"category":cat,"pattern":name,"atlas_roi_index":roi,"atlas_label":label})
            ax.axhline(0,color="grey",ls=":",lw=.8)
            ax.set_title(f"{cat}: {label} (ROI {roi})\n{name}",fontsize=9)
            if c==3: ax.set_xlabel("Lag (s)")
            if col==0: ax.set_ylabel("Mean residual difference")
            if c==0 and col==0: ax.legend(fontsize=8)
    fig.suptitle("Selected after analysis: first unused atlas ROI in each available pattern\n"
                 "Model-specific simultaneous bands; no inferential comparison between models",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.95))
    finish(fig,figdir/"selected_region_curves.png")
    save_csv(pd.DataFrame(examples),outdir/"selected_examples.csv")
    return pd.DataFrame(rows),consensus


def anatomy(results, atlas, outdir, figdir):
    img=nib.load(outdir/"atlas/canlab2018_labels.nii")
    labels=np.rint(img.get_fdata()).astype(int)
    assert set(range(1,490)).issubset(set(np.unique(labels)))
    mapdir=outdir/"anatomy"; mapdir.mkdir(exist_ok=True)
    for c,cat in enumerate(CATEGORIES):
        fig,axs=plt.subplots(3,4,figsize=(12,8),constrained_layout=True)
        for k,r in enumerate(results):
            counts=np.stack([(r["signed"][c*NROI:(c+1)*NROI]>0).sum(1),
                             (r["signed"][c*NROI:(c+1)*NROI]<0).sum(1)],axis=1)
            lookup=np.vstack([np.zeros((1,2)),counts])
            volume=lookup[labels].astype(np.float32)
            nifti=nib.Nifti1Image(volume,img.affine,img.header)
            nifti.set_data_dtype(np.float32)
            nib.save(nifti,mapdir/f"{r['model']}_{cat}_supported_lag_counts_pos_neg.nii.gz")
            for j,zmm in enumerate([-40,-10,20,50]):
                # Atlas affine is axis aligned; find axial voxel at chosen MNI z.
                z=int(round(nib.affines.apply_affine(np.linalg.inv(img.affine),[0,0,zmm])[2]))
                z=np.clip(z,0,labels.shape[2]-1)
                ax=axs[k,j]
                background=np.where(labels[:,:,z]>0,.85,np.nan)
                ax.imshow(background.T,origin="lower",interpolation="nearest",cmap="gray",vmin=0,vmax=1)
                vals=volume[:,:,z,:].sum(axis=-1).T
                im=ax.imshow(np.ma.masked_where(vals==0,vals),origin="lower",interpolation="nearest",cmap="viridis",vmin=0,vmax=41)
                ax.axis("off"); ax.set_title(f"{r['model']} | z={zmm} mm",fontsize=9)
        fig.colorbar(im,ax=axs.ravel().tolist(),label="Number of supported sampled lags (either direction)",shrink=.6)
        fig.suptitle(f"{cat}: descriptive parcel maps, not voxelwise tests\nRepresentative axial slices; atlas-index mapping follows project convention",fontsize=12)
        finish(fig,figdir/f"anatomy_{cat}.png")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[4])
    parser.add_argument("--output",type=Path,default=Path("outputs/wm_scb_downstream_20261026"))
    parser.add_argument("--cluster-bootstrap",type=int,default=100)
    args=parser.parse_args(); root=args.root.resolve()
    outdir=args.output.resolve(); outdir.mkdir(parents=True,exist_ok=True)
    figdir=outdir/"figures"; figdir.mkdir(exist_ok=True)
    atlas=pd.read_csv(outdir/"atlas/atlas_regions.csv")
    assert len(atlas)==489 and atlas.atlas_roi_index.is_unique
    results=[]; cluster_stats=[]; cluster_assign=[]; cluster_candidates=[]
    with threadpool_limits(limits=2):
        for model in MODELS:
            result=audit_model(root,outdir,atlas,model)
            temporal_plots(result,figdir)
            cs,ca,cc=functional_clusters(result,outdir,figdir,args.cluster_bootstrap)
            cluster_stats.extend(cs); cluster_assign.extend(ca); cluster_candidates.extend(cc)
            del result["D"]
            results.append(result)
    evidence=pd.concat([r["evidence"] for r in results],ignore_index=True)
    save_csv(evidence,outdir/"regional_evidence.csv")
    save_csv(pd.DataFrame([s for r in results for s in r["segments"]]),outdir/"supported_sample_runs.csv")
    save_csv(pd.DataFrame([r["audit"] for r in results]),outdir/"run_audit.csv")
    save_csv(pd.DataFrame(cluster_stats),outdir/"clustering_summary.csv")
    save_csv(pd.DataFrame(cluster_assign),outdir/"region_cluster_assignments.csv")
    save_csv(pd.DataFrame(cluster_candidates),outdir/"cluster_k_candidates.csv")
    summary=[]; time_summary=[]; network_summary=[]
    endpoint_summary=[]; window_summary=[]
    for r in results:
        ev=r["evidence"]
        for c,cat in enumerate(CATEGORIES):
            ee=ev[ev.category==cat]
            signed=r["signed"][c*NROI:(c+1)*NROI]
            support=signed!=0
            counts=support.sum(axis=0)
            endpoint_summary.append({"model":r["model"],"category":cat,
                "regions_any":int(support.any(axis=1).sum()),
                "regions_any_excluding_final_sample":int(support[:,:-1].any(axis=1).sum()),
                "regions_only_final_sample":int((support[:,-1]&~support[:,:-1].any(axis=1)).sum()),
                "regions_only_boundary_samples":int((support.any(axis=1)&~support[:,1:-1].any(axis=1)).sum()),
                "peak_support_count":int(counts.max()),"peak_support_lag_sec":TIMES[counts.argmax()]})
            for label,ii in [("0_to_under_7.2",np.arange(10)),("7.2_to_under_14.4",np.arange(10,20)),
                             ("14.4_to_under_21.6",np.arange(20,30)),("21.6_through_28.8",np.arange(30,41))]:
                window_summary.append({"model":r["model"],"category":cat,"descriptive_window":label,
                    "positive_supported_cells":int((signed[:,ii]>0).sum()),
                    "negative_supported_cells":int((signed[:,ii]<0).sum()),
                    "regions_any":int(support[:,ii].any(axis=1).sum())})
            summary.append({"model":r["model"],"category":cat,"regions":489,
                "omnibus_regions":int(ee.significant.sum()),"band_regions":int(ee.band_any.sum()),
                "positive_regions":int((ee.band_positive_lags>0).sum()),
                "negative_regions":int((ee.band_negative_lags>0).sum()),
                "biphasic_regions":int(ee.band_biphasic.sum()),
                "omnibus_only_regions":int((ee.evidence_class=="omnibus_only").sum()),
                "bands_only_regions":int((ee.evidence_class=="bands_only").sum()),
                "neither_regions":int((ee.evidence_class=="neither").sum()),
                "positive_supported_cells":int((signed==1).sum()),
                "negative_supported_cells":int((signed==-1).sum()),
                "median_supported_lags_among_supported_regions":float(ee.loc[ee.band_any,"band_supported_lags"].median()),
                "median_mean_curve_rms_all_regions":float(ee.mean_curve_rms.median()),
                "regions_supported_at_lag_zero":int((signed[:,0]!=0).sum()),
                "regions_supported_at_last_lag":int((signed[:,-1]!=0).sum())})
            for t in range(NT):
                time_summary.append({"model":r["model"],"category":cat,"time_sec":TIMES[t],
                    "positive_regions":int((signed[:,t]==1).sum()),"negative_regions":int((signed[:,t]==-1).sum()),
                    "supported_regions":int((signed[:,t]!=0).sum())})
            for network,nn in ee.groupby("atlas_network",dropna=False):
                network_summary.append({"model":r["model"],"category":cat,"atlas_network":network,
                    "regions_in_network":len(nn),"supported_regions":int(nn.band_any.sum()),
                    "supported_region_fraction":float(nn.band_any.mean()),
                    "positive_regions":int((nn.band_positive_lags>0).sum()),
                    "negative_regions":int((nn.band_negative_lags>0).sum()),
                    "biphasic_regions":int(nn.band_biphasic.sum()),
                    "supported_cells":int(nn.band_supported_lags.sum()),
                    "supported_cell_fraction":float(nn.band_supported_lags.sum()/(len(nn)*41)),
                    "median_curve_rms":float(nn.mean_curve_rms.median())})
    summary=pd.DataFrame(summary); ts=pd.DataFrame(time_summary); nets=pd.DataFrame(network_summary)
    assert np.array_equal(summary.omnibus_regions+summary.bands_only_regions,
                          summary.band_regions+summary.omnibus_only_regions)
    assert np.array_equal(summary.positive_regions+summary.negative_regions-summary.biphasic_regions,
                          summary.band_regions)
    assert np.array_equal(summary.band_regions+summary.omnibus_only_regions+summary.neither_regions,
                          summary.regions)
    save_csv(summary,outdir/"category_summary.csv"); save_csv(ts,outdir/"temporal_support_counts.csv")
    save_csv(nets,outdir/"network_summary.csv")
    save_csv(pd.DataFrame(endpoint_summary),outdir/"endpoint_support_summary.csv")
    save_csv(pd.DataFrame(window_summary),outdir/"descriptive_window_summary.csv")
    pairwise,consensus=model_comparisons(results,outdir,figdir)
    anatomy(results,atlas,outdir,figdir)
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,constrained_layout=True)
    for c,cat in enumerate(CATEGORIES):
        ax=axes.flat[c]
        for k,m in enumerate(MODELS):
            ss=ts[(ts.model==m)&(ts.category==cat)]
            ax.plot(ss.time_sec,ss.positive_regions,color=COLORS[k],label=m+" positive")
            ax.plot(ss.time_sec,-ss.negative_regions,color=COLORS[k],ls="--",label=m+" negative")
        ax.axhline(0,color="grey",lw=.8); ax.set_title(cat); ax.set_xlabel("Lag (s)")
        ax.set_ylabel("Supported regions (+positive / −negative)"); ax.set_ylim(-489,489)
    axes.flat[0].legend(fontsize=7,ncol=2)
    fig.suptitle("Number of regions with simultaneous bands excluding zero at each sampled lag")
    finish(fig,figdir/"temporal_support_counts.png")
    fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    x=np.arange(4)
    for k,m in enumerate(MODELS):
        ss=summary[summary.model==m].set_index("category").loc[CATEGORIES]
        ax.bar(x+(k-1)*.24,ss.band_regions,width=.22,color=COLORS[k],label=m)
        for i,n in enumerate(ss.band_regions): ax.text(x[i]+(k-1)*.24,n+5,str(n),ha="center",fontsize=9)
    ax.set_xticks(x,CATEGORIES); ax.set_ylim(0,520); ax.set_ylabel("Regions with at least one supported sampled lag")
    ax.set_title("Time-localized evidence by category and HRF model"); ax.legend()
    finish(fig,figdir/"supported_region_counts.png")
    # Cross-model cluster agreement is descriptive; labels are arbitrary.
    cluster_agreement=[]
    for mode in ["amplitude","shape"]:
        for a,b in [(0,1),(0,2),(1,2)]:
            cluster_agreement.append({"mode":mode,"model_a":MODELS[a],"model_b":MODELS[b],
                "adjusted_rand_index":adjusted_rand_score(results[a]["clusters"][mode]["labels"],results[b]["clusters"][mode]["labels"])})
    save_csv(pd.DataFrame(cluster_agreement),outdir/"cluster_hrf_agreement.csv")
    assignments=pd.DataFrame(cluster_assign)
    composition=assignments.groupby(["model","mode","cluster","atlas_network"]).size().rename("regions").reset_index()
    totals=assignments.groupby(["model","mode","cluster"]).size().rename("cluster_size").reset_index()
    composition=composition.merge(totals,on=["model","mode","cluster"])
    composition["within_cluster_fraction"]=composition.regions/composition.cluster_size
    save_csv(composition,outdir/"cluster_network_composition.csv")
    # Network plot uses full denominators and reports descriptive proportions.
    network_order=nets.groupby("atlas_network").regions_in_network.first().sort_values(ascending=False).head(20).index
    mat=nets.groupby(["atlas_network","model"]).supported_region_fraction.mean().unstack().reindex(index=network_order,columns=MODELS)
    fig,ax=plt.subplots(figsize=(7,9),constrained_layout=True)
    im=ax.imshow(mat,aspect="auto",interpolation="nearest",vmin=0,vmax=1,cmap="viridis")
    ax.set_xticks(range(3),MODELS); ax.set_yticks(range(len(mat)),mat.index)
    ax.set_title("Regional support by atlas network\nMean fraction across four categories; 20 largest groups",fontsize=11)
    fig.colorbar(im,ax=ax,label="Fraction with at least one supported lag",shrink=.6)
    finish(fig,figdir/"network_support.png")
    payload={"audit":[r["audit"] for r in results],"category_summary":summary.to_dict("records"),
             "pairwise_hrf":pairwise.to_dict("records"),"clustering":cluster_stats,
             "cluster_hrf_agreement":cluster_agreement,
             "consensus_counts":consensus.groupby("category").agg(
                regions_with_same_direction_all_models=("same_positive_lags_all_models",lambda x:int(((x+consensus.loc[x.index,"same_negative_lags_all_models"])>0).sum())),
                regions_with_opposite_supported_directions=("opposite_supported_direction_lags",lambda x:int((x>0).sum()))).reset_index().to_dict("records")}
    (outdir/"analysis_summary.json").write_text(json.dumps(payload,indent=2),encoding="utf-8")
    print(summary.to_string(index=False),flush=True)
    print("COMPLETE:",outdir,flush=True)


if __name__=="__main__":
    main()

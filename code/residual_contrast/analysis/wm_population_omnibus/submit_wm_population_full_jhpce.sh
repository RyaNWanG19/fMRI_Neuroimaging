#!/bin/bash
# Three independent model jobs; Holm scope is 489 ROIs per model/contrast.
# Independence confirmed by the data owner in the analysis conversation.
# No manifest was supplied; row IDs are local to each model file.
set -euo pipefail
export REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
export DATA_DIR="${DATA_DIR:-/users/rwang}"
export OUTPUT_ROOT="${OUTPUT_ROOT:-/users/rwang/jhpce_outputs_wm_population_omnibus_full}"
export SMOKE_TEST=0 NPERM=99999 MAKE_PLOTS=0
export RUN_DIAGNOSTICS=1 NBOOTSTRAP=99999
export INDEPENDENT_SUBJECTS_VERIFIED=1
export NULL_SYMMETRY_VERIFIED=0 NULL_SYMMETRY_ASSUMED=1
export RESAMPLING_VERIFICATION_NOTE='Data owner confirms distinct unrelated participants (2026-09-20 conversation). No subject manifest available. Whole-curve null symmetry accepted as a modeling assumption; not empirically verified.'
export CONTRAST_NAME=Body_LoadDiff_vs_Face_LoadDiff
mkdir -p /users/rwang/logs
for HRF_MODEL in cHRF cHRFderiv sHRF; do
    export HRF_MODEL
    sbatch --export=ALL "${REPO_DIR}/code/residual_contrast/analysis/wm_population_omnibus/run_wm_population_omnibus_jhpce.sbatch"
done

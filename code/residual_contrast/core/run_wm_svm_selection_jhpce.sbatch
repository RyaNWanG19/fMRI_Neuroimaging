#!/bin/bash
#SBATCH --job-name=wm_svm_rfe
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/users/rwang/logs/wm_svm_rfe_%A_%a.log
#SBATCH --error=/users/rwang/logs/wm_svm_rfe_%A_%a.err
set -euo pipefail
export REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
export DATA_DIR="${DATA_DIR:-/users/rwang}"
export OUTPUT_ROOT="${OUTPUT_ROOT:-${REPO_DIR}/outputs/wm_svm_selection}"
export NBOOT="${NBOOT:-100}"
export INNER_FOLDS="${INNER_FOLDS:-5}"
export FEATURE_COUNTS="${FEATURE_COUNTS:-10 25 50 100 200 Inf}"
export RUN_IMPORTANCE="${RUN_IMPORTANCE:-1}"
export IMPORTANCE_REPEATS="${IMPORTANCE_REPEATS:-20}"
models=(cHRF cHRFderiv sHRF)
index="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$index" =~ ^[0-2]$ ]] || { echo 'Expected array index 0-2' >&2; exit 2; }
export HRF_MODEL="${models[$index]}"
export DATA_FILE="${DATA_DIR}/WM${HRF_MODEL}.mat"
MATLAB_BIN="${MATLAB_BIN:-/jhpce/shared/jhpce/core/matlab/R2025a/bin/matlab}"
[[ -x "$MATLAB_BIN" ]] || { echo 'Set MATLAB_BIN to an executable MATLAB path' >&2; exit 127; }
[[ -f "$DATA_FILE" ]] || { echo "Missing data: $DATA_FILE" >&2; exit 2; }
cd "$REPO_DIR"
"$MATLAB_BIN" -batch "addpath(fullfile(getenv('REPO_DIR'),'code','residual_contrast','core')); run_wm_svm_selection_jhpce_env;"

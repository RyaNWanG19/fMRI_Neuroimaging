#!/bin/bash
#SBATCH --job-name=wm_load_omnibus
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=/users/rwang/logs/wm_load_%A_%a.log
#SBATCH --error=/users/rwang/logs/wm_load_%A_%a.err
set -euo pipefail
export REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
export DATA_DIR="${DATA_DIR:-/users/rwang}"
export OUTPUT_ROOT="${OUTPUT_ROOT:-${REPO_DIR}/data/task_residual/wm_load_omnibus}"
export WM_ANALYSIS=load-differences
models=(cHRF cHRFderiv sHRF)
index="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$index" =~ ^[0-2]$ ]] || { echo 'Expected array index 0-2' >&2; exit 2; }
export HRF_MODEL="${models[$index]}"
# Reuse executable checks, assumption gates, environment setup and batch quoting.
bash "${REPO_DIR}/code/residual_contrast/analysis/wm_population_omnibus/run_wm_population_omnibus_jhpce.sbatch"

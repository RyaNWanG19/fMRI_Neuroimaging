#!/bin/bash
#SBATCH --job-name=wm_omnibus
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rwang/logs/wm_omnibus_%j.log
#SBATCH --error=/users/rwang/logs/wm_omnibus_%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
DATA_DIR="${DATA_DIR:-${REPO_DIR}/data/task_residual}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${DATA_DIR}/wm_population_omnibus}"
HRF_MODEL="${HRF_MODEL:-cHRF}"
CONTRAST_NAME="${CONTRAST_NAME:-Body_LoadDiff_vs_Face_LoadDiff}"
SMOKE_TEST="${SMOKE_TEST:-0}"
MAKE_PLOTS="${MAKE_PLOTS:-0}"
NPERM="${NPERM:-}"
SUBJECT_IDS_FILE="${SUBJECT_IDS_FILE:-}"
INDEPENDENT_SUBJECTS_VERIFIED="${INDEPENDENT_SUBJECTS_VERIFIED:-0}"
NULL_SYMMETRY_VERIFIED="${NULL_SYMMETRY_VERIFIED:-0}"
NULL_SYMMETRY_ASSUMED="${NULL_SYMMETRY_ASSUMED:-0}"
RUN_DIAGNOSTICS="${RUN_DIAGNOSTICS:-0}"
NBOOTSTRAP="${NBOOTSTRAP:-99999}"
RESAMPLING_VERIFICATION_NOTE="${RESAMPLING_VERIFICATION_NOTE:-}"
MATLAB_BIN="${MATLAB_BIN:-/jhpce/shared/jhpce/core/matlab/R2025a/bin/matlab}"

# Slurm batch shells may not initialize Environment Modules even when the
# login shell does. Use the same MATLAB module convention as the project's
# earlier JHPCE jobs, while permitting an explicit executable override.
if [[ -z "${MATLAB_BIN}" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        if [[ -r /etc/profile.d/modules.sh ]]; then
            # shellcheck disable=SC1091
            source /etc/profile.d/modules.sh
        elif [[ -r /usr/share/Modules/init/bash ]]; then
            # shellcheck disable=SC1091
            source /usr/share/Modules/init/bash
        fi
    fi

    if ! command -v matlab >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
        if module load matlab 2>/dev/null; then
            :
        elif module load matlab/R2025a 2>/dev/null; then
            :
        fi
    fi

    if command -v matlab >/dev/null 2>&1; then
        MATLAB_BIN="$(command -v matlab)"
    elif [[ -f /users/rwang/MATLAB && -x /users/rwang/MATLAB ]]; then
        MATLAB_BIN=/users/rwang/MATLAB
    elif [[ -x /users/rwang/MATLAB/bin/matlab ]]; then
        MATLAB_BIN=/users/rwang/MATLAB/bin/matlab
    fi
fi

if [[ -z "${MATLAB_BIN}" || ! -x "${MATLAB_BIN}" ]]; then
    echo "MATLAB executable unavailable on compute node." >&2
    echo "Run 'module avail matlab' or submit with MATLAB_BIN=/full/path/to/matlab." >&2
    exit 127
fi

case "${HRF_MODEL,,}" in
    chrf)
        DATA_FILE="${DATA_DIR}/WMcHRF.mat"
        ;;
    chrfderiv)
        DATA_FILE="${DATA_DIR}/WMcHRFderiv.mat"
        ;;
    shrf)
        DATA_FILE="${DATA_DIR}/WMsHRF.mat"
        ;;
    *)
        echo "HRF_MODEL must be cHRF, cHRFderiv, or sHRF; got: ${HRF_MODEL}" >&2
        exit 2
        ;;
esac

if [[ ! -f "${DATA_FILE}" ]]; then
    echo "WM residual data file not found: ${DATA_FILE}" >&2
    exit 2
fi

if [[ "${SMOKE_TEST}" != "1" ]]; then
    if [[ -n "${SUBJECT_IDS_FILE}" && ! -f "${SUBJECT_IDS_FILE}" ]]; then
        echo "SUBJECT_IDS_FILE was specified but does not exist." >&2
        exit 2
    fi
    if [[ "${INDEPENDENT_SUBJECTS_VERIFIED}" != "1" ]]; then
        echo "Set INDEPENDENT_SUBJECTS_VERIFIED=1 only after verifying independent subject rows." >&2
        exit 2
    fi
    if [[ "${NULL_SYMMETRY_VERIFIED}" != "1" && "${NULL_SYMMETRY_ASSUMED}" != "1" ]]; then
        echo "Set NULL_SYMMETRY_ASSUMED=1 to document the accepted modeling assumption." >&2
        exit 2
    fi
    if [[ -z "${RESAMPLING_VERIFICATION_NOTE}" ]]; then
        echo "RESAMPLING_VERIFICATION_NOTE is required for a full inferential run." >&2
        exit 2
    fi
fi

mkdir -p "${OUTPUT_ROOT}"
cd "${REPO_DIR}"

export REPO_DIR DATA_FILE OUTPUT_ROOT HRF_MODEL CONTRAST_NAME SMOKE_TEST
export MAKE_PLOTS NPERM SUBJECT_IDS_FILE INDEPENDENT_SUBJECTS_VERIFIED
export NULL_SYMMETRY_VERIFIED RESAMPLING_VERIFICATION_NOTE
export NULL_SYMMETRY_ASSUMED RUN_DIAGNOSTICS NBOOTSTRAP

echo "MATLAB_BIN: ${MATLAB_BIN}"

MATLAB_RUNNER="addpath(fullfile(getenv('REPO_DIR'),'code','residual_contrast','analysis','wm_population_omnibus')); run_wm_population_omnibus_jhpce_env;"
"${MATLAB_BIN}" -batch "${MATLAB_RUNNER}"

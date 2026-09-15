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
RESAMPLING_VERIFICATION_NOTE="${RESAMPLING_VERIFICATION_NOTE:-}"

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
    if [[ -z "${SUBJECT_IDS_FILE}" || ! -f "${SUBJECT_IDS_FILE}" ]]; then
        echo "A valid SUBJECT_IDS_FILE is required for a full inferential run." >&2
        exit 2
    fi
    if [[ "${INDEPENDENT_SUBJECTS_VERIFIED}" != "1" ]]; then
        echo "Set INDEPENDENT_SUBJECTS_VERIFIED=1 only after verifying independent subject rows." >&2
        exit 2
    fi
    if [[ "${NULL_SYMMETRY_VERIFIED}" != "1" ]]; then
        echo "Set NULL_SYMMETRY_VERIFIED=1 only after justifying sign-flip null symmetry." >&2
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

matlab -batch "
moduleDir = fullfile(getenv('REPO_DIR'), 'code', 'residual_contrast', ...
    'analysis', 'wm_population_omnibus');
addpath(moduleDir);

cfg = struct();
cfg.repoRoot = getenv('REPO_DIR');
cfg.dataFile = getenv('DATA_FILE');
cfg.hrfModelName = getenv('HRF_MODEL');
cfg.contrastName = getenv('CONTRAST_NAME');
cfg.smokeTest = strcmp(getenv('SMOKE_TEST'), '1');
cfg.makePlots = strcmp(getenv('MAKE_PLOTS'), '1');

nPermText = getenv('NPERM');
if ~isempty(nPermText)
    cfg.nPerm = str2double(nPermText);
    assert(isfinite(cfg.nPerm) && cfg.nPerm >= 1 && cfg.nPerm == floor(cfg.nPerm), ...
        'NPERM must be a positive integer.');
end

jobID = getenv('SLURM_JOB_ID');
if isempty(jobID)
    jobID = datestr(now, 'yyyymmddTHHMMSSFFF');
end
cfg.outputDir = fullfile(getenv('OUTPUT_ROOT'), cfg.contrastName, ...
    cfg.hrfModelName, ['job_' jobID]);

if ~cfg.smokeTest
    idPayload = load(getenv('SUBJECT_IDS_FILE'));
    candidateNames = {'subjectIDs', 'subject_ids', 'subjIDs'};
    subjectIDs = [];
    for k = 1:numel(candidateNames)
        if isfield(idPayload, candidateNames{k})
            subjectIDs = idPayload.(candidateNames{k});
            break
        end
    end
    assert(~isempty(subjectIDs), ...
        'SUBJECT_IDS_FILE must contain subjectIDs, subject_ids, or subjIDs.');
    cfg.subjectIDs = subjectIDs;
end

cfg.resampling.independentSubjectsVerified = ...
    strcmp(getenv('INDEPENDENT_SUBJECTS_VERIFIED'), '1');
cfg.resampling.nullSymmetryVerified = ...
    strcmp(getenv('NULL_SYMMETRY_VERIFIED'), '1');
cfg.resampling.verificationNote = getenv('RESAMPLING_VERIFICATION_NOTE');

results = run_wm_population_omnibus(cfg);
fprintf('WM temporal-omnibus results saved to: %s\\n', results.config.outputDir);
"

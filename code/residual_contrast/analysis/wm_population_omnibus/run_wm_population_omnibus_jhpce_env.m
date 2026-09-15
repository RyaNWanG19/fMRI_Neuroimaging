function results = run_wm_population_omnibus_jhpce_env()
%RUN_WM_POPULATION_OMNIBUS_JHPCE_ENV Build configuration from Slurm env.

repoDir = required_env('REPO_DIR');
dataFile = required_env('DATA_FILE');
outputRoot = required_env('OUTPUT_ROOT');

cfg = struct();
cfg.repoRoot = repoDir;
cfg.dataFile = dataFile;
cfg.hrfModelName = required_env('HRF_MODEL');
cfg.contrastName = required_env('CONTRAST_NAME');
cfg.smokeTest = env_flag('SMOKE_TEST');
cfg.makePlots = env_flag('MAKE_PLOTS');

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
cfg.outputDir = fullfile(outputRoot, cfg.contrastName, ...
    cfg.hrfModelName, ['job_' jobID]);

if ~cfg.smokeTest
    idPayload = load(required_env('SUBJECT_IDS_FILE'));
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

cfg.resampling.independentSubjectsVerified = env_flag('INDEPENDENT_SUBJECTS_VERIFIED');
cfg.resampling.nullSymmetryVerified = env_flag('NULL_SYMMETRY_VERIFIED');
cfg.resampling.verificationNote = getenv('RESAMPLING_VERIFICATION_NOTE');

results = run_wm_population_omnibus(cfg);
fprintf('WM temporal-omnibus results saved to: %s\n', results.config.outputDir);
end

function value = required_env(name)
value = getenv(name);
if isempty(value)
    error('run_wm_population_omnibus_jhpce_env:MissingEnvironmentVariable', ...
        'Required environment variable %s is empty.', name);
end
end

function tf = env_flag(name)
tf = strcmp(getenv(name), '1');
end

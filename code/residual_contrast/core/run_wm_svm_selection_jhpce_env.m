function run_wm_svm_selection_jhpce_env()
% Cluster runner: independent output per HRF model and Slurm job/array task.
repoDir=required_env('REPO_DIR');
dataFile=required_env('DATA_FILE');
outputRoot=required_env('OUTPUT_ROOT');
model=required_env('HRF_MODEL');
cfg=struct('runFeatureSelection',true);
cfg.featureSelection=struct('innerFolds',env_number('INNER_FOLDS',5));
text=getenv('FEATURE_COUNTS');
if ~isempty(text)
    % Parse plain numbers without evaluating arbitrary MATLAB expressions.
    tokens=regexp(strtrim(text),'[ ,]+','split');
    cfg.featureSelection.featureCounts=cellfun(@str2double,tokens);
end
cfg.featureSelection=svm_rfe_config(cfg.featureSelection);
cfg.importanceRepeats=env_number('IMPORTANCE_REPEATS',20);
flag=getenv('RUN_IMPORTANCE');
if isempty(flag), flag='1'; end
assert(any(strcmp(flag,{'0','1'})), 'RUN_IMPORTANCE must be 0 or 1.');
cfg.runImportance=strcmp(flag,'1');
nBoot=env_number('NBOOT',100);
validateattributes(nBoot,{'numeric'},{'scalar','integer','positive'});
setenv('NBOOT',num2str(nBoot));
workers=env_number('SLURM_CPUS_PER_TASK',2);
validateattributes(workers,{'numeric'},{'scalar','integer','positive'});
maxNumCompThreads(1);
jobID=getenv('SLURM_JOB_ID');
if isempty(jobID), jobID=datestr(now,'yyyymmddTHHMMSSFFF'); end
taskID=getenv('SLURM_ARRAY_TASK_ID');
if isempty(taskID), taskID='0'; end
tag=['job_' jobID '_task_' taskID];
outputDir=fullfile(outputRoot,model,tag);
assert(~exist(outputDir,'dir'), 'Output directory already exists: %s',outputDir);
mkdir(outputDir);
addpath(fullfile(repoDir,'code','residual_contrast','core'));
pool=gcp('nocreate');
if isempty(pool)
    cluster=parcluster('local');
    storage=fullfile(tempdir,['svm_rfe_' tag]);
    if ~exist(storage,'dir'), mkdir(storage); end
    cluster.JobStorageLocation=storage;
    parpool(cluster,min(workers,cluster.NumWorkers));
end
fprintf('Model=%s; outer bootstraps=%d; inner folds=%d; outputs=%s\n', ...
    model,nBoot,cfg.featureSelection.innerFolds,outputDir);
PredictContrast_WM(dataFile,outputDir,cfg);
end

function text=required_env(name)
text=getenv(name);
assert(~isempty(text),'Required environment variable missing: %s',name);
end

function value=env_number(name,default)
text=getenv(name);
if isempty(text), value=default; else, value=str2double(text); end
validateattributes(value,{'numeric'},{'scalar','finite'},mfilename,name);
end

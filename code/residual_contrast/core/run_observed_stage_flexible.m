function run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, TASKNAME, cfg)
% run_observed_stage_flexible
% Stage A (Observed): subject-bootstrap (OOB) pairwise decoding per timepoint,
% linear SVM weights, and Haufe patterns.
%
% UPDATED OUTPUTS
%   - TA_*             : temporal accuracy
%   - PMZ_*            : mean paired score difference divided by its SD
%   - PairDiff_*       : raw paired score difference summaries
%   - TA_subj_oob      : subject-level OOB accuracy trajectories
%   - PairDiff_subj_oob: subject-level OOB paired score-difference trajectories
%
% INPUTS
%   DATA_FILE   : .mat file containing variable 'Results'
%   OUTDIR      : output directory
%   COMP2       : Kx2 matrix of pairwise condition indices (e.g., [1 2; 1 3; 2 3])
%                 If empty and C==2, it will default to [1 2].
%   cond_labels : cellstr of length C (optional; used for labels only)
%   TASKNAME    : char/string, for output file naming and metadata
%   cfg         : optional struct; runImportance=false, importanceRepeats=20,
%                 importanceSeed=2718. Importance is OOB paired half-swap
%                 accuracy loss, not population inference. Set
%                 runFeatureSelection=true for nested subject-grouped SVM-RFE;
%                 see svm_rfe_config for cfg.featureSelection defaults.

if nargin < 6, cfg = struct(); end
if ~isfield(cfg,'runImportance'), cfg.runImportance = false; end
if ~isfield(cfg,'importanceRepeats'), cfg.importanceRepeats = 20; end
if ~isfield(cfg,'importanceSeed'), cfg.importanceSeed = 2718; end
if ~isfield(cfg,'runFeatureSelection'), cfg.runFeatureSelection = false; end
if ~isfield(cfg,'featureSelection'), cfg.featureSelection = struct(); end
cfg.featureSelection = svm_rfe_config(cfg.featureSelection);
validateattributes(cfg.runImportance, {'logical','numeric'}, {'scalar','binary'});
validateattributes(cfg.runFeatureSelection, {'logical','numeric'}, {'scalar','binary'});
validateattributes(cfg.importanceRepeats, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.importanceSeed, {'numeric'}, {'scalar','integer','nonnegative','<=',2^32-1});

if nargin < 2 || isempty(OUTDIR), OUTDIR = fullfile(pwd,'outputs'); end
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end
if nargin < 5 || isempty(TASKNAME), TASKNAME = 'task'; end

fprintf('[%s] Stage A (Observed OOB) starting...\n', char(TASKNAME));

%% ---------- 1) Load & canonicalize data ----------
tmp = load(DATA_FILE, 'Results', 'AlreadyCanonical');
assert(isfield(tmp,'Results'), 'Results missing in MAT file: %s', DATA_FILE);
Rraw = tmp.Results;
sz = size(Rraw); sz(end+1:4) = 1;

if isfield(tmp,'AlreadyCanonical') && tmp.AlreadyCanonical
    % Wrapper guarantees internal canonical [S x R x T x C]
    R = Rraw;
elseif sz(1)==489 && sz(2)==41
    % File is [R x T x S x C] -> permute to [S x R x T x C]
    R = permute(Rraw, [3 1 2 4]);
elseif sz(2)==489 && sz(3)==41
    % Already [S x R x T x C]
    R = Rraw;
else
    % Fallback robust inference
    R = toR(Rraw);
end
clear Rraw tmp

[S, Rg, T, C] = size(R);
fprintf('Data dims (internal): S=%d, R=%d, T=%d, C=%d\n', S, Rg, T, C);

% Hard guards (adapt if needed per task)
assert(Rg >= 200 && Rg <= 600, 'Regions (Rg) out of range: %d', Rg);
assert(T  >= 20  && T  <= 1000, 'Timepoints (T) out of range: %d', T);
assert(C  >= 2   && C  <= 12,   'Conditions (C) out of range: %d', C);

% COMP2 default for C==2
if (~exist('COMP2','var') || isempty(COMP2)) && C==2
    COMP2 = [1 2];
end
assert(all(COMP2(:) >= 1 & COMP2(:) <= C), 'COMP2 contains invalid condition indices.');

nCon = size(COMP2, 1);

%% ---------- 2) Bootstrap design (OOB) ----------
% Bootstrap iterations from env, default 100
nBootEnv = str2double(getenv('NBOOT'));
nBoot = 100;
if ~isnan(nBootEnv) && nBootEnv > 0, nBoot = nBootEnv; end
fprintf('Bootstrap iterations: %d\n', nBoot);

minOOB = max(20, ceil(0.10 * S));

rng(42, 'twister');  % reproducible splits
split_list = cell(nBoot,1);
for b = 1:nBoot
    while true
        boot_idx = randi(S, S, 1);           % WITH replacement
        oob_mask = true(S,1);
        oob_mask(boot_idx) = false;
        out_of_bag = find(oob_mask);
        if numel(out_of_bag) >= minOOB, break; end
    end
    split_list{b}.train_boot = boot_idx;     % duplicates kept
    split_list{b}.test       = out_of_bag;   % unique OOB subjects
end

%% ---------- 3) Preallocation ----------
TA_obs       = zeros(T, nCon, nBoot);                    % [T x Con x B]
PMZ_obs      = nan(T, nCon, nBoot);                     % [T x Con x B]
PairDiff_obs = zeros(T, nCon, nBoot);                   % [T x Con x B]

W_obs     = zeros(Rg, T, nCon, nBoot, 'single');        % [R x T x Con x B]
Haufe_obs = zeros(Rg, T, nCon, nBoot, 'single');
Beta_obs  = zeros(T, nBoot, nCon, Rg, 'single');        % [T x B x Con x R]
% NaN marks features dropped by training preprocessing or failed fits.
Reliance_obs = nan(Rg, T, nCon, nBoot, 'single');
OOB_count = zeros(nBoot,1);
FS_boot = cell(nBoot,1);

% NEW: subject-level OOB trajectories
TA_subj_obs       = nan(S, T, nCon, nBoot, 'single');   % [S x T x Con x B]
PairDiff_subj_obs = nan(S, T, nCon, nBoot, 'single');   % [S x T x Con x B]

%% ---------- 4) Main loop ----------
% Rconst = parallel.pool.Constant(@() R);

parfor b = 1:nBoot
    % Rloc = Rconst.Value;
    Rloc = R;

    train_boot = split_list{b}.train_boot;   % [S x 1], with duplicates
    testSubs   = split_list{b}.test;         % OOB unique
    OOB_count(b) = numel(testSubs);

    nTr = numel(train_boot);
    nTe = numel(testSubs);
    Ytr = [ones(nTr,1); -ones(nTr,1)];
    Yte = [ones(nTe,1); -ones(nTe,1)]; %#ok<NASGU>

    TA_b       = zeros(T, nCon);
    PMZ_b      = nan(T, nCon);
    PairDiff_b = zeros(T, nCon);

    W_b = zeros(Rg, T, nCon, 'single');
    H_b = zeros(Rg, T, nCon, 'single');
    Reliance_b = nan(Rg, T, nCon, 'single');
    FS_b = [];
    if cfg.runFeatureSelection
        FS_b = empty_selection_boot(Rg,T,nCon);
    end

    % NEW: subject-level OOB trajectories for this bootstrap
    TA_subj_b       = nan(S, T, nCon, 'single');
    PairDiff_subj_b = nan(S, T, nCon, 'single');

    for k = 1:nCon
        A  = COMP2(k,1);
        Bc = COMP2(k,2);

        % Slice contrast once
        X1_AllT = Rloc(:,:,:,A);             % [S x R x T]
        X2_AllT = Rloc(:,:,:,Bc);

        for t = 1:T
            % Subjects x Regions at time t
            X1t = X1_AllT(:,:,t);
            X2t = X2_AllT(:,:,t);

            % Build train/test (keep bootstrap multiplicity in train)
            Xtr = [X1t(train_boot,:); X2t(train_boot,:)];
            Xte = [X1t(testSubs,:);   X2t(testSubs,:)];

            % Train-only standardization + NaN impute + constant-drop
            [XtrZ, mu, sd, keep, ok_tr] = standardize_impute_train(double(Xtr));
            if ~ok_tr
                TA_b(t,k)       = NaN;
                PMZ_b(t,k)      = NaN;
                PairDiff_b(t,k) = NaN;
                % subject-level arrays stay NaN
                % maps remain zeros
                continue
            end

            XteZ = standardize_impute_apply(double(Xte), mu, sd, keep);

            mdl = fitcsvm(XtrZ, Ytr, ...
                'KernelFunction', 'linear', ...
                'BoxConstraint', 1, ...
                'ClassNames', [-1 1]);

            [yhat, score] = predict(mdl, XteZ);

            if size(score,2) == 1
                s = score;
            else
                poscol = find(mdl.ClassNames == 1, 1);
                s = score(:, poscol);
            end

            % Split held-out predictions into positive and negative observations
            yhat_pos = yhat(1:nTe);
            yhat_neg = yhat(nTe+1:end);

            s_pos = s(1:nTe);
            s_neg = s(nTe+1:end);

            % Subject-level accuracy for each held-out pair:
            % 1   = both correct
            % 0.5 = one correct
            % 0   = neither correct
            acc_subj = 0.5 * ((yhat_pos == 1) + (yhat_neg == -1));

            % Subject-level paired score difference
            d_subj = s_pos - s_neg;

            % Bootstrap-level summaries
            TA_b(t,k)       = mean(acc_subj, 'omitnan');
            PairDiff_b(t,k) = mean(d_subj, 'omitnan');
            PMZ_b(t,k)      = compute_pmz(d_subj);

            % NEW: save subject-level OOB values for this bootstrap
            TA_subj_b(testSubs, t, k)       = single(acc_subj);
            PairDiff_subj_b(testSubs, t, k) = single(d_subj);

            % Raw SVM weights
            w_keep = mdl.Beta;
            w_full = zeros(Rg,1,'single');
            w_full(keep) = single(w_keep);

            % Haufe (fast): a = (X' * (X*w)) / (n-1) in z-space, then .* SD_train
            v = XtrZ * w_keep;                           % [2*nTr x 1]
            a_keep = (XtrZ' * v) / max(1, size(XtrZ,1)-1);
            a_full = zeros(Rg,1,'single');
            a_full(keep) = single(a_keep .* sd(keep)'); % back to raw units

            W_b(:,t,k) = w_full;
            H_b(:,t,k) = a_full;

            % Independent per-fit stream; also used for selected-model reliance.
            swapSeed = mod(double(cfg.importanceSeed) + ...
                (b-1)*T*nCon + (k-1)*T + t-1, 2^32);
            if cfg.runImportance
                loss = svm_paired_region_reliance(mdl, XteZ, ...
                    cfg.importanceRepeats, swapSeed);
                Reliance_b(keep,t,k) = single(loss);
            end

            if cfg.runFeatureSelection
                fsCfg = cfg.featureSelection;
                % Same inner subject folds across lags/contrasts in a given
                % outer bootstrap, independent of all OOB observations.
                fsCfg.seed = mod(double(fsCfg.seed)+b-1,2^32);
                [selectedMdl, prep, info] = svm_nested_rfe( ...
                    double(X1t(train_boot,:)), double(X2t(train_boot,:)), ...
                    train_boot, fsCfg);
                FS_b.uniqueSubjectIDs = info.uniqueSubjectIDs;
                FS_b.subjectFold = info.subjectFold;
                FS_b.audit{t,k} = rmfield(info, ...
                    {'selected','eligible','uniqueSubjectIDs','subjectFold','rowFold'});
                if info.valid
                    selected = info.selected;
                    % Apply only the outer-training preprocessing and mask.
                    testX = double(Xte(:,selected));
                    muSelected = prep.mu(selected);
                    sdSelected = prep.sd(selected);
                    for j=1:size(testX,2)
                        testX(isnan(testX(:,j)),j) = muSelected(j);
                    end
                    testZ = (testX-muSelected)./sdSelected;
                    [yp,sp] = predict(selectedMdl,testZ);
                    FS_b.accuracy(t,k) = mean(0.5*((yp(1:nTe)==1) + ...
                        (yp(nTe+1:end)==-1)));
                    FS_b.scoreDifference(t,k) = mean(sp(1:nTe,2)-sp(nTe+1:end,2));
                    FS_b.selected(:,t,k) = selected(:);
                    FS_b.eligible(:,t,k) = info.eligible(:);
                    FS_b.valid(t,k) = true;
                    FS_b.count(t,k) = info.chosenCount;
                    FS_b.requestedCount(t,k) = info.chosenRequestedCount;
                    FS_b.innerAccuracy(t,k) = info.innerAccuracy;
                    FS_b.weights(selected,t,k) = single(selectedMdl.Beta);
                    if cfg.runImportance
                        FS_b.reliance(selected,t,k) = single( ...
                            svm_paired_region_reliance(selectedMdl,testZ, ...
                            cfg.importanceRepeats,swapSeed));
                    end
                end
            end
        end
    end

    TA_obs(:,:,b)         = TA_b;
    PMZ_obs(:,:,b)        = PMZ_b;
    PairDiff_obs(:,:,b)   = PairDiff_b;

    TA_subj_obs(:,:,:,b)       = TA_subj_b;
    PairDiff_subj_obs(:,:,:,b) = PairDiff_subj_b;

    W_obs(:,:,:,b)     = W_b;
    Haufe_obs(:,:,:,b) = H_b;
    Reliance_obs(:,:,:,b) = Reliance_b;
    FS_boot{b} = FS_b;
    Beta_obs(:, b, :, :) = permute(W_b, [2 3 1]);  % [T x Con x R] -> [T x 1 x Con x R]
end

%% ---------- 5) Aggregate ----------
TA_mean = mean(TA_obs, 3, 'omitnan');                 % [T x Con]
TA_se   = std(TA_obs,  0, 3, 'omitnan');              % bootstrap SD (= SE of statistic)

PMZ_mean = mean(PMZ_obs, 3, 'omitnan');               % [T x Con]
PMZ_se   = std(PMZ_obs,  0, 3, 'omitnan');

PairDiff_mean = mean(PairDiff_obs, 3, 'omitnan');     % [T x Con]
PairDiff_se   = std(PairDiff_obs,  0, 3, 'omitnan');

% Subject-level OOB trajectories averaged across the bootstraps in which
% each subject was held out
TA_subj_oob       = mean(TA_subj_obs,       4, 'omitnan');   % [S x T x Con]
PairDiff_subj_oob = mean(PairDiff_subj_obs, 4, 'omitnan');   % [S x T x Con]

TA_subj_hits       = sum(~isnan(TA_subj_obs),       4);      % [S x T x Con]
PairDiff_subj_hits = sum(~isnan(PairDiff_subj_obs), 4);      % [S x T x Con]

Haufe_mean = mean(Haufe_obs, 4, 'omitnan');          % [R x T x Con]
Haufe_std  = std(Haufe_obs,  0, 4, 'omitnan');       % bootstrap SD (= SE of statistic)

Reliance_mean = mean(Reliance_obs, 4, 'omitnan');
Reliance_resample_sd = std(Reliance_obs, 0, 4, 'omitnan');
Reliance_valid_boot = sum(isfinite(Reliance_obs), 4);

% Labels
if ~exist('cond_labels','var') || isempty(cond_labels)
    cond_labels = arrayfun(@(i) sprintf('cond%d', i), 1:C, 'uni', 0);
end
contrast_labels = arrayfun(@(i) sprintf('%s_vs_%s', ...
    safe_label(cond_labels, COMP2(i,1)), safe_label(cond_labels, COMP2(i,2))), ...
    (1:nCon)', 'uni', 0);

FeatureSelection = struct('enabled',logical(cfg.runFeatureSelection), ...
    'config',cfg.featureSelection);
if cfg.runFeatureSelection
    FeatureSelection.bootstrap = FS_boot;
    FeatureSelection.method = 'nested subject-grouped batched linear SVM-RFE; C fixed; max inner mean accuracy, smallest count on ties';
    FeatureSelection.interpretation = 'predictive selection stability, not regional significance';
    accuracy = cellfun(@(x) x.accuracy,FS_boot,'UniformOutput',false);
    accuracy = cat(3,accuracy{:});
    valid = cellfun(@(x) x.valid,FS_boot,'UniformOutput',false);
    valid = cat(3,valid{:});
    selected = cellfun(@(x) x.selected,FS_boot,'UniformOutput',false);
    selected = cat(4,selected{:});
    eligible = cellfun(@(x) x.eligible,FS_boot,'UniformOutput',false);
    eligible = cat(4,eligible{:});
    FeatureSelection.accuracy_obs = accuracy;
    FeatureSelection.accuracy_mean = mean(accuracy,3,'omitnan');
    FeatureSelection.accuracy_difference_obs = accuracy-TA_obs;
    FeatureSelection.accuracy_difference_mean = mean(accuracy-TA_obs,3,'omitnan');
    FeatureSelection.valid_boot = sum(valid,3);
    FeatureSelection.selected_boot = sum(selected,4);
    FeatureSelection.eligible_boot = sum(eligible,4);
    denominator = reshape(FeatureSelection.valid_boot,[1 T nCon]);
    FeatureSelection.selection_frequency = FeatureSelection.selected_boot./denominator;
    FeatureSelection.selection_frequency_when_eligible = ...
        FeatureSelection.selected_boot./FeatureSelection.eligible_boot;
end

%% ---------- 6) Save ----------
outf = fullfile(OUTDIR, ['observed_' lower(char(TASKNAME)) '.mat']);

meta.timestamp    = datestr(now,30);
meta.task         = char(TASKNAME);
meta.size_R       = [S Rg T C];
meta.nBoot        = nBoot;
meta.oob_min      = minOOB;
meta.oob_counts   = OOB_count;
meta.comp2        = COMP2;
meta.haufe_method = 'X''*(X*w)/(n-1) in z-space, then .* SD_train';
meta.standardize  = 'train-only mean/std; NaN->train-mean; drop constant features';
meta.importance = cfg;
meta.importance.method = 'OOB accuracy loss under independent 50% within-subject condition swaps, one region at a time';
meta.importance.interpretation = 'fitted decoder reliance; no p-values, causal interpretation, or feature selection';
meta.importance.uncertainty = 'resample SD is descriptive; overlapping OOB sets are not independent replicates';
meta.featureSelection = cfg.featureSelection;
meta.featureSelection.enabled = logical(cfg.runFeatureSelection);
meta.metrics      = struct( ...
    'TA', 'mean subject-level OOB paired accuracy (0, 0.5, 1)', ...
    'PMZ', 'mean paired score difference divided by its sample SD within each bootstrap OOB set', ...
    'PairDiff', 'mean paired score difference within each bootstrap OOB set', ...
    'PairDiff_subj_oob', 'subject-level OOB paired score difference, averaged across bootstraps where subject was OOB');
meta.software     = struct('matlab', version, 'parallel', tryget(@() ver('parallel')));
try
    p = gcp('nocreate');
    meta.nWorkers = p.NumWorkers;
catch
    meta.nWorkers = NaN;
end

save(outf, ...
    'TA_obs','TA_se','TA_mean', ...
    'PMZ_obs','PMZ_se','PMZ_mean', ...
    'PairDiff_obs','PairDiff_se','PairDiff_mean', ...
    'TA_subj_oob','TA_subj_hits', ...
    'PairDiff_subj_oob','PairDiff_subj_hits', ...
    'W_obs','Haufe_obs','Haufe_mean','Haufe_std', ...
    'Reliance_obs','Reliance_mean','Reliance_resample_sd','Reliance_valid_boot', ...
    'FeatureSelection','split_list', ...
    'Beta_obs','COMP2','contrast_labels','cond_labels', ...
    'S','Rg','T','C','nBoot','TASKNAME','meta', '-v7.3');

if cfg.runImportance
    [roi_index, time_index, contrast_index] = ndgrid(1:Rg, 1:T, 1:nCon);
    contrast = reshape(contrast_labels(contrast_index(:)), [], 1);
    accuracy_loss = double(Reliance_mean(:));
    resample_sd = double(Reliance_resample_sd(:));
    valid_boot = Reliance_valid_boot(:);
    summary = table(roi_index(:), time_index(:), contrast, ...
        accuracy_loss, resample_sd, valid_boot, 'VariableNames', ...
        {'roi_index','time_index','contrast','accuracy_loss','resample_sd','valid_boot'});
    writetable(summary, fullfile(OUTDIR, ['regional_reliance_' lower(char(TASKNAME)) '.csv']));
end

if cfg.runFeatureSelection
    [roi_index,time_index,contrast_index] = ndgrid(1:Rg,1:T,1:nCon);
    contrast = reshape(contrast_labels(contrast_index(:)),[],1);
    selection_frequency = FeatureSelection.selection_frequency(:);
    selected_boot = FeatureSelection.selected_boot(:);
    eligible_boot = FeatureSelection.eligible_boot(:);
    valid_boot = repmat(reshape(FeatureSelection.valid_boot,[1 T nCon]),Rg,1,1);
    summary = table(roi_index(:),time_index(:),contrast,selection_frequency, ...
        selected_boot,eligible_boot,valid_boot(:),'VariableNames', ...
        {'roi_index','time_index','contrast','selection_frequency', ...
        'selected_boot','eligible_boot','valid_boot'});
    writetable(summary,fullfile(OUTDIR,['regional_selection_' lower(char(TASKNAME)) '.csv']));
    [time_index,contrast_index] = ndgrid(1:T,1:nCon);
    contrast = reshape(contrast_labels(contrast_index(:)),[],1);
    baseline_accuracy = TA_mean(:);
    selected_accuracy = FeatureSelection.accuracy_mean(:);
    paired_accuracy_difference = FeatureSelection.accuracy_difference_mean(:);
    valid_boot = FeatureSelection.valid_boot(:);
    comparison = table(time_index(:),contrast,baseline_accuracy,selected_accuracy, ...
        paired_accuracy_difference,valid_boot,'VariableNames', ...
        {'time_index','contrast','baseline_accuracy','selected_accuracy', ...
        'paired_accuracy_difference','valid_boot'});
    writetable(comparison,fullfile(OUTDIR,['selection_accuracy_' lower(char(TASKNAME)) '.csv']));
end

fprintf('[%s] Stage A saved: %s\n', char(TASKNAME), outf);
end

%% ================= helpers =================
function fs = empty_selection_boot(Rg,T,nCon)
fs = struct('selected',false(Rg,T,nCon),'eligible',false(Rg,T,nCon), ...
    'valid',false(T,nCon),'count',nan(T,nCon),'requestedCount',nan(T,nCon), ...
    'innerAccuracy',nan(T,nCon),'accuracy',nan(T,nCon), ...
    'scoreDifference',nan(T,nCon),'weights',zeros(Rg,T,nCon,'single'), ...
    'reliance',nan(Rg,T,nCon,'single'),'uniqueSubjectIDs',[], 'subjectFold',[]);
fs.audit = cell(T,nCon);
end

function R = toR(Results)
% Robustly infer [S x R x T x C] from arbitrary 4D layout.
sz = size(Results);
sz(end+1:4) = 1;

% If it already looks canonical, accept as-is
if sz(2) >= 200 && sz(4) <= 12 && sz(3) >= 20
    R = Results;
    return;
end

cands = perms(1:4);
best = [];
bestscore = -inf;

for i = 1:size(cands,1)
    p = cands(i,:);
    Xsz = sz(p);

    Ssz = Xsz(1);
    Rsz = Xsz(2);
    Tsz = Xsz(3);
    Csz = Xsz(4);

    ok = (Ssz >= 20) && (Rsz >= 50) && (Tsz >= 20) && (Csz <= 12);
    if ~ok, continue; end

    sc = 0;
    if Rsz == 489, sc = sc + 100; end           % strong prior for atlas size
    sc = sc + min(Rsz,600)/10;                  % prefer larger R
    sc = sc + (12 - Csz);                       % prefer small C
    sc = sc + max(0, 400 - abs(Tsz-300))/100;   % plausible T
    sc = sc + max(0, 500 - abs(Ssz-300))/200;   % plausible S

    if sc > bestscore
        bestscore = sc;
        best = p;
    end
end

if ~isempty(best)
    R = permute(Results, best);
else
    warning('toR: fallback to raw order; could not confidently infer dims.');
    R = Results;
end
end

function [XtrZ, mu_full, sd_full, keepIdx, ok] = standardize_impute_train(Xtr)
% Train-only standardization with NaN impute to train mean; drop constants.
mu_full = mean(Xtr, 1, 'omitnan');
Xtr_f   = Xtr;

nanmask = isnan(Xtr_f);
if any(nanmask(:))
    for j = 1:size(Xtr_f,2)
        if any(nanmask(:,j))
            Xtr_f(nanmask(:,j),j) = mu_full(j);
        end
    end
end

sd_full = std(Xtr_f, 0, 1);
keepIdx = ~isnan(mu_full) & (sd_full > 0);

if ~any(keepIdx)
    ok   = false;
    XtrZ = [];
else
    ok = true;
    mu_full(~keepIdx) = 0;
    sd_full(~keepIdx) = 1;
    XtrZ = (Xtr_f(:,keepIdx) - mu_full(keepIdx)) ./ sd_full(keepIdx);
end
end

function XteZ = standardize_impute_apply(Xte, mu_full, sd_full, keepIdx)
Xte_f = Xte;

nanmask = isnan(Xte_f);
if any(nanmask(:))
    for j = 1:size(Xte_f,2)
        if any(nanmask(:,j))
            Xte_f(nanmask(:,j),j) = mu_full(j);
        end
    end
end

XteZ = (Xte_f(:,keepIdx) - mu_full(keepIdx)) ./ sd_full(keepIdx);
end

function z = compute_pmz(d)
% Mean paired score difference divided by its sample SD.
% Returns NaN if there are too few observations or zero variance.
d = d(:);
d = d(~isnan(d));

if numel(d) < 2
    z = NaN;
    return;
end

mu = mean(d);
sdv = std(d, 0);   % sample SD

if sdv <= eps
    if abs(mu) <= eps
        z = 0;
    else
        z = NaN;
    end
else
    z = mu / sdv;
end
end

function s = safe_label(lbls, i)
if ~isempty(lbls) && i <= numel(lbls)
    s = lbls{i};
else
    s = sprintf('cond%d', i);
end
end

function v = tryget(f)
try
    v = feval(f);
catch
    v = [];
end
end

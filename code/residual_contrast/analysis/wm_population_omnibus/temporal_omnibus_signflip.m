function results = temporal_omnibus_signflip(D, subjectIDs, roiIDs, roiNames, timeSec, cfg)
%TEMPORAL_OMNIBUS_SIGNFLIP ROI-wise temporal omnibus sign-flip test.
%
% Tests, for every ROI, whether the population mean paired-difference curve
% is zero at all analyzed times. The statistic is S = sum_t T(t)^2 and its
% null distribution is generated with one sign per independent subject,
% shared across every ROI and time point.
%
% Required inputs
%   D          paired differences in the order specified by
%              cfg.inputDimensionOrder (default/inferred output is S x R x T)
%   subjectIDs one ID per subject
%   roiIDs     one ID per ROI (may be [])
%   roiNames   one name per ROI (may be [])
%   timeSec    one time value in seconds per analyzed sample
%   cfg        configuration struct; important fields are:
%       nPerm                  default 99999
%       seed                   default 20260914
%       alpha                  default 0.05
%       outputDir              required when saveOutputs=true
%       contrastName, hrfModelName
%       inputDimensionOrder    e.g. 'roi-time-subject'
%       runMode                'inference', 'smoke-test', or 'verification'
%       resampling.scheme      must be 'independent-subject-signflip'
%       resampling.independentSubjectsVerified (required true for inference)
%
% Complete-case missingness is fixed separately for each ROI: a subject is
% retained only if all analyzed times are finite. Untestable ROIs remain in
% the Holm family with internal p=1. No classifier observations are formed.

if nargin < 6 || ~isstruct(cfg)
    error('temporal_omnibus_signflip:InvalidConfig', 'cfg must be a struct.');
end
cfg = fill_defaults(cfg);
validate_resampling(cfg);

[D, inputInfo] = canonicalize_input(D, subjectIDs, timeSec, cfg);
D = double(D);
[nSubject, nROI, nTime] = size(D);

subjectIDs = text_vector(subjectIDs, nSubject, 'subjectIDs');
if numel(unique(subjectIDs)) ~= nSubject
    error('temporal_omnibus_signflip:DuplicateSubjectIDs', ...
        'subjectIDs must uniquely identify the independent observational units.');
end
timeSec = double(timeSec(:));
if numel(timeSec) ~= nTime || any(~isfinite(timeSec))
    error('temporal_omnibus_signflip:InvalidTime', ...
        'timeSec must contain one finite value per time point.');
end
if isempty(roiIDs)
    roiIDs = string((1:nROI)');
else
    roiIDs = text_vector(roiIDs, nROI, 'roiIDs');
end
if isempty(roiNames)
    roiNames = string(compose('Parcel_%03d', (1:nROI)'));
else
    roiNames = text_vector(roiNames, nROI, 'roiNames');
end

nanMask = isnan(D);
infMask = isinf(D);
subjectMask = all(isfinite(D), 3); % [subject x ROI], fixed across permutations
effectiveN = sum(subjectMask, 1)';
testable = effectiveN >= 2;

meanD = nan(nROI, nTime);
sdD = nan(nROI, nTime);
tObs = nan(nROI, nTime);
sObs = nan(nROI, 1);
zeroVarianceObserved = false(nROI, nTime);
nonzeroMeanZeroVarianceObserved = false(nROI, nTime);
excludedSubjectIDs = cell(nROI, 1);

for r = 1:nROI
    keep = subjectMask(:, r);
    excludedSubjectIDs{r} = subjectIDs(~keep);
    if ~testable(r)
        continue
    end
    x = reshape(D(keep, r, :), effectiveN(r), nTime);
    meanD(r, :) = mean(x, 1);
    sdD(r, :) = std(x, 0, 1); % sample SD, N-1 normalization
    zeroVarianceObserved(r, :) = (sdD(r, :) == 0);
    finiteVariance = ~zeroVarianceObserved(r, :);
    tObs(r, finiteVariance) = meanD(r, finiteVariance) ./ ...
        (sdD(r, finiteVariance) ./ sqrt(effectiveN(r)));
    bothZero = zeroVarianceObserved(r, :) & (meanD(r, :) == 0);
    nonzeroDegenerate = zeroVarianceObserved(r, :) & (meanD(r, :) ~= 0);
    tObs(r, bothZero) = 0;
    tObs(r, nonzeroDegenerate) = sign(meanD(r, nonzeroDegenerate)) .* Inf;
    nonzeroMeanZeroVarianceObserved(r, :) = nonzeroDegenerate;
    sObs(r) = sum(tObs(r, :) .^ 2);
end

% sum(X.^2) is invariant to sign flipping. Excluded subject/ROI curves are
% set to zero after the observed summaries are computed. This permits one
% batched matrix multiplication across ROIs even when complete-case masks
% differ; multiplication by the corresponding global subject sign still
% selects exactly the retained signs for each ROI.
sumSq = nan(nROI, nTime);
for r = find(testable)'
    x = reshape(D(subjectMask(:, r), r, :), effectiveN(r), nTime);
    sumSq(r, :) = sum(x .^ 2, 1);
end
testableROI = find(testable);
Dfixed = D;
for r = 1:nROI
    Dfixed(~subjectMask(:, r), r, :) = 0;
end

exceedanceCount = zeros(nROI, 1);
permutationDegenerateCellCount = zeros(nROI, 1);
roundoffClampCellCount = zeros(nROI, 1);
if cfg.storePermutationStatistics
    sPerm = nan(cfg.nPerm, nROI);
else
    sPerm = [];
end
if cfg.storePermutationSigns
    permutationSigns = zeros(cfg.nPerm, nSubject, 'int8');
else
    permutationSigns = [];
end

stream = RandStream('mt19937ar', 'Seed', cfg.seed);
nBatches = ceil(cfg.nPerm / cfg.batchSize);
if cfg.verbose
    fprintf('Temporal omnibus: %d subjects, %d ROIs, %d times, %d permutations.\n', ...
        nSubject, nROI, nTime, cfg.nPerm);
end

for batch = 1:nBatches
    firstPerm = (batch - 1) * cfg.batchSize + 1;
    lastPerm = min(cfg.nPerm, batch * cfg.batchSize);
    nThis = lastPerm - firstPerm + 1;

    % rand is generated [subject x permutation] so a permutation occupies a
    % contiguous RNG segment and results do not depend on cfg.batchSize.
    signs = ones(nSubject, nThis);
    signs(rand(stream, nSubject, nThis) < 0.5) = -1;
    signs = signs'; % [permutation x subject]
    if cfg.storePermutationSigns
        permutationSigns(firstPerm:lastPerm, :) = int8(signs);
    end

    for rStart = 1:cfg.roiBatchSize:numel(testableROI)
        rStop = min(numel(testableROI), rStart + cfg.roiBatchSize - 1);
        ridx = testableROI(rStart:rStop);
        nR = numel(ridx);
        x = reshape(Dfixed(:, ridx, :), nSubject, nR * nTime);
        signedSum = signs * x;
        fixedSumSq = reshape(sumSq(ridx, :), 1, []);
        nByCell = reshape(repmat(effectiveN(ridx), 1, nTime), 1, []);

        varianceNumerator = nByCell .* fixedSumSq - signedSum .^ 2;
        roundoffTolerance = 64 .* eps(max(nByCell .* fixedSumSq, 1));
        substantiveNegative = varianceNumerator < -roundoffTolerance;
        if any(substantiveNegative(:))
            error('temporal_omnibus_signflip:NegativeVariance', ...
                ['Optimized variance calculation produced a negative value ', ...
                 'beyond floating-point tolerance.']);
        end
        roundoffMask = varianceNumerator < 0;
        varianceNumerator(roundoffMask) = 0;
        zeroVariance = (varianceNumerator == 0);

        tSquared = zeros(size(signedSum));
        positiveVariance = varianceNumerator > 0;
        tNumerator = (nByCell - 1) .* signedSum .^ 2;
        tSquared(positiveVariance) = tNumerator(positiveVariance) ./ ...
            varianceNumerator(positiveVariance);
        tSquared(zeroVariance & signedSum ~= 0) = Inf;
        sBatch = sum(reshape(tSquared, nThis, nR, nTime), 3);

        exceedanceCount(ridx) = exceedanceCount(ridx) + ...
            sum(sBatch >= reshape(sObs(ridx), 1, []), 1)';
        degenerateByROI = squeeze(sum(sum(reshape(zeroVariance, nThis, nR, nTime), 1), 3));
        clampedByROI = squeeze(sum(sum(reshape(roundoffMask, nThis, nR, nTime), 1), 3));
        permutationDegenerateCellCount(ridx) = permutationDegenerateCellCount(ridx) + degenerateByROI(:);
        roundoffClampCellCount(ridx) = roundoffClampCellCount(ridx) + clampedByROI(:);
        if cfg.storePermutationStatistics
            sPerm(firstPerm:lastPerm, ridx) = sBatch;
        end
    end

    if cfg.verbose && (batch == 1 || batch == nBatches || mod(batch, cfg.progressEveryBatches) == 0)
        fprintf('  completed %d/%d permutations\n', lastPerm, cfg.nPerm);
    end
end

% p=1 for untestable ROIs preserves the planned correction-family size.
pRaw = ones(nROI, 1);
pRaw(testable) = (1 + exceedanceCount(testable)) ./ (cfg.nPerm + 1);
pHolm = holm_adjust(pRaw);
significant = testable & (pHolm <= cfg.alpha);

mcSE = sqrt(pRaw .* (1 - pRaw) ./ (cfg.nPerm + 1));
pLow = max(1 / (cfg.nPerm + 1), pRaw - 2 .* mcSE);
pHigh = min(1, pRaw + 2 .* mcSE);
pLow(~testable) = 1;
pHigh(~testable) = 1;
holmLow = holm_adjust(pLow);
holmHigh = holm_adjust(pHigh);
monteCarloNearThreshold = testable & holmLow <= cfg.alpha & holmHigh >= cfg.alpha;

status = repmat("ok", nROI, 1);
status(effectiveN < nSubject & testable) = "finite_complete_case_exclusions";
status(~testable) = "untestable_fewer_than_2_complete_subjects";
hasObservedDegeneracy = any(zeroVarianceObserved, 2) & testable;
for r = find(hasObservedDegeneracy)'
    status(r) = status(r) + ";zero_variance_degeneracy";
end

inferenceValid = strcmpi(cfg.runMode, 'inference') && ...
    cfg.resampling.independentSubjectsVerified && cfg.resampling.nullSymmetryVerified;
if inferenceValid
    inferenceStatus = "inferential";
else
    inferenceStatus = "non-inferential_" + string(cfg.runMode);
end

correctionFamily = sprintf(['All %d planned ROIs for contrast %s and HRF model %s; ', ...
    'untestable ROIs retained with p=1.'], nROI, cfg.contrastName, cfg.hrfModelName);

summaryTable = table(roiIDs, roiNames, effectiveN, sObs, pRaw, pHolm, significant, ...
    repmat(string(cfg.contrastName), nROI, 1), repmat(string(cfg.hrfModelName), nROI, 1), ...
    status, repmat(inferenceStatus, nROI, 1), monteCarloNearThreshold, ...
    'VariableNames', {'roi_id','roi_name','effective_n','S_obs','p_raw','p_holm', ...
    'significant','contrast','hrf_model','data_quality_status','inference_status', ...
    'monte_carlo_near_threshold'});

results = struct();
results.mean_D = meanD;
results.sd_D = sdD;
results.T_obs = tObs;
results.S_obs = sObs;
results.p_raw = pRaw;
results.p_holm = pHolm;
results.significant = significant;
results.effectiveN = effectiveN;
results.subjectMasks = subjectMask;
results.excludedSubjectIDs = excludedSubjectIDs;
results.subjectIDs = subjectIDs;
results.roiIDs = roiIDs;
results.roiNames = roiNames;
results.timeSec = timeSec;
results.config = cfg;
results.inputMapping = inputInfo;
results.correctionFamily = correctionFamily;
results.minimumAttainableP = 1 / (cfg.nPerm + 1);
results.monteCarloSE_raw = mcSE;
results.monteCarloNearThreshold = monteCarloNearThreshold;
results.monteCarloNote = ['Raw-p Monte Carlo SE is reported. The near-threshold flag ', ...
    'is a two-SE sensitivity screen propagated through Holm; decisions close to ', ...
    'the corrected threshold can change with another finite permutation sample.'];
results.inferenceValid = inferenceValid;
results.inferenceStatus = inferenceStatus;
results.dataQuality = struct( ...
    'totalNaN', sum(nanMask(:)), ...
    'totalInf', sum(infMask(:)), ...
    'nanByROI', squeeze(sum(sum(nanMask, 1), 3)), ...
    'infByROI', squeeze(sum(sum(infMask, 1), 3)), ...
    'zeroVarianceObserved', zeroVarianceObserved, ...
    'nonzeroMeanZeroVarianceObserved', nonzeroMeanZeroVarianceObserved, ...
    'permutationDegenerateCellCount', permutationDegenerateCellCount, ...
    'roundoffClampCellCount', roundoffClampCellCount, ...
    'status', status);
results.summaryTable = summaryTable;
results.S_perm = sPerm;
results.permutationSigns = permutationSigns;
results.outputFiles = struct();

if cfg.saveOutputs
    results = save_outputs(results);
end
end

function cfg = fill_defaults(cfg)
cfg = set_default(cfg, 'nPerm', 99999);
cfg = set_default(cfg, 'seed', 20260914);
cfg = set_default(cfg, 'alpha', 0.05);
cfg = set_default(cfg, 'outputDir', '');
cfg = set_default(cfg, 'contrastName', 'unspecified_contrast');
cfg = set_default(cfg, 'hrfModelName', 'unspecified_HRF');
cfg = set_default(cfg, 'inputDimensionOrder', '');
cfg = set_default(cfg, 'runMode', 'inference');
cfg = set_default(cfg, 'batchSize', 200);
cfg = set_default(cfg, 'roiBatchSize', 64);
cfg = set_default(cfg, 'progressEveryBatches', 10);
cfg = set_default(cfg, 'verbose', true);
cfg = set_default(cfg, 'saveOutputs', true);
cfg = set_default(cfg, 'makePlots', true);
cfg = set_default(cfg, 'includePointwiseBands', true);
cfg = set_default(cfg, 'rankPlotTopN', 50);
cfg = set_default(cfg, 'curvesPerPage', 12);
cfg = set_default(cfg, 'storePermutationStatistics', false);
cfg = set_default(cfg, 'storePermutationSigns', false);
if ~isfield(cfg, 'resampling') || ~isstruct(cfg.resampling)
    cfg.resampling = struct();
end
cfg.resampling = set_default(cfg.resampling, 'scheme', 'independent-subject-signflip');
cfg.resampling = set_default(cfg.resampling, 'independentSubjectsVerified', false);
cfg.resampling = set_default(cfg.resampling, 'nullSymmetryVerified', false);
cfg.resampling = set_default(cfg.resampling, 'verificationNote', '');

validateattributes(cfg.nPerm, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.seed, {'numeric'}, {'scalar','integer','nonnegative'});
validateattributes(cfg.alpha, {'numeric'}, {'scalar','>',0,'<',1});
validateattributes(cfg.batchSize, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.roiBatchSize, {'numeric'}, {'scalar','integer','positive'});
end

function validate_resampling(cfg)
if ~strcmpi(cfg.resampling.scheme, 'independent-subject-signflip')
    error('temporal_omnibus_signflip:UnsupportedResampling', ...
        ['No dependence-aware resampling scheme is justified by this repository. ', ...
         'Only independent-subject-signflip is implemented.']);
end
if strcmpi(cfg.runMode, 'inference')
    if ~isequal(cfg.resampling.independentSubjectsVerified, true)
        error('temporal_omnibus_signflip:IndependenceUnresolved', ...
            ['Full inference requires cfg.resampling.independentSubjectsVerified=true ', ...
             'after verifying that rows are independent subjects (for example, an ', ...
             'unrelated HCP sample). Use runMode=''smoke-test'' only for non-inferential checks.']);
    end
    if ~isequal(cfg.resampling.nullSymmetryVerified, true)
        error('temporal_omnibus_signflip:NullSymmetryUnresolved', ...
            ['Exact subject-level sign flipping also requires an appropriate ', ...
             'symmetry assumption for the paired-difference curves. Set ', ...
             'cfg.resampling.nullSymmetryVerified=true only with a documented justification.']);
    end
    if isempty(strtrim(char(string(cfg.resampling.verificationNote))))
        error('temporal_omnibus_signflip:MissingVerificationNote', ...
            'Document the independence evidence in cfg.resampling.verificationNote.');
    end
elseif ~any(strcmpi(cfg.runMode, {'smoke-test','verification'}))
    error('temporal_omnibus_signflip:InvalidRunMode', ...
        'cfg.runMode must be inference, smoke-test, or verification.');
end
end

function [Dout, info] = canonicalize_input(Din, subjectIDs, timeSec, cfg)
if ~isnumeric(Din) || ~isreal(Din)
    error('D must be a real numeric array.');
end
sz = size(Din);
sz(end+1:3) = 1;
if numel(sz) > 3 && any(sz(4:end) ~= 1)
    error('D must have exactly three non-singleton dimensions.');
end
sz = sz(1:3);

if isempty(cfg.inputDimensionOrder)
    candidates = perms(1:3);
    ok = false(size(candidates, 1), 1);
    for i = 1:size(candidates, 1)
        proposed = sz(candidates(i, :));
        ok(i) = proposed(1) == numel(subjectIDs) && proposed(3) == numel(timeSec);
    end
    candidates = candidates(ok, :);
    if size(candidates, 1) ~= 1
        error('temporal_omnibus_signflip:AmbiguousDimensions', ...
            ['Could not uniquely infer D axes from subjectIDs and timeSec. ', ...
             'Set cfg.inputDimensionOrder explicitly.']);
    end
    permutation = candidates(1, :);
    orderLabel = 'inferred from subjectIDs and timeSec';
else
    token = regexprep(lower(char(cfg.inputDimensionOrder)), '[^a-z]', '');
    switch token
        case {'subjectroitime','subjectsroistime','srt'}
            permutation = [1 2 3];
        case {'subjecttimeroi','subjectstimerois','str'}
            permutation = [1 3 2];
        case {'roisubjecttime','roissubjectstime','rst'}
            permutation = [2 1 3];
        case {'roitimesubject','roistimesubjects','rts'}
            permutation = [3 1 2];
        case {'timesubjectroi','timesubjectsrois','tsr'}
            permutation = [2 3 1];
        case {'timeroisubject','timeroissubjects','trs'}
            permutation = [3 2 1];
        otherwise
            error('Unknown cfg.inputDimensionOrder: %s', cfg.inputDimensionOrder);
    end
    orderLabel = char(cfg.inputDimensionOrder);
end
Dout = permute(Din, permutation);
if size(Dout, 1) ~= numel(subjectIDs) || size(Dout, 3) ~= numel(timeSec)
    error('Dimension mapping does not agree with subjectIDs/timeSec lengths.');
end
info = struct('originalSize', sz, 'inputOrder', orderLabel, ...
    'permutationToSubjectsROIsTime', permutation, 'canonicalSize', size(Dout));
end

function out = text_vector(in, expectedN, label)
if ischar(in) && expectedN == 1
    out = string({in});
else
    try
        out = string(in(:));
    catch
        error('%s must be convertible to a string vector.', label);
    end
end
if numel(out) ~= expectedN || any(ismissing(out))
    error('%s must contain %d nonmissing values.', label, expectedN);
end
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

function results = save_outputs(results)
cfg = results.config;
if isempty(cfg.outputDir)
    error('cfg.outputDir is required when cfg.saveOutputs=true.');
end
if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

matFile = fullfile(cfg.outputDir, 'temporal_omnibus_results.mat');
csvFile = fullfile(cfg.outputDir, 'temporal_omnibus_roi_summary.csv');
rankFile = fullfile(cfg.outputDir, 'ranked_roi_omnibus_evidence.png');
protected = {matFile, csvFile};
if cfg.makePlots
    protected{end+1} = rankFile; %#ok<AGROW>
end
for i = 1:numel(protected)
    if isfile(protected{i})
        error('Refusing to overwrite existing result: %s', protected{i});
    end
end

writetable(results.summaryTable, csvFile);
results.outputFiles.csv = csvFile;
if cfg.makePlots
    make_rank_plot(results, rankFile);
    results.outputFiles.rankedPlot = rankFile;
    results.outputFiles.curvePlots = make_curve_plots(results, cfg.outputDir);
end
results.outputFiles.mat = matFile;
save(matFile, 'results', '-v7.3');
end

function make_rank_plot(results, outFile)
[~, order] = sort(results.p_holm, 'ascend');
nShow = min(results.config.rankPlotTopN, numel(order));
order = order(1:nShow);
y = -log10(results.p_holm(order));
isSig = results.significant(order);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 650]);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
hold on;
scatter(1:nShow, y, 42, [0.20 0.45 0.75], 'filled');
scatter(find(isSig), y(isSig), 54, [0.80 0.15 0.15], 'filled');
yline(-log10(results.config.alpha), '--k', sprintf('Holm alpha = %.3g', results.config.alpha));
xlabel(sprintf('ROI rank (top %d of %d)', nShow, numel(results.p_holm)));
ylabel('-log_{10}(Holm-adjusted temporal omnibus p)');
title(sprintf('%s | %s', results.config.contrastName, results.config.hrfModelName), ...
    'Interpreter', 'none');
grid on;
for j = 1:min(15, nShow)
    text(j, y(j), " " + results.roiNames(order(j)), 'FontSize', 7, ...
        'Interpreter', 'none', 'Rotation', 35);
end
exportgraphics(fig, outFile, 'Resolution', 180);
end

function files = make_curve_plots(results, outputDir)
sigIdx = find(results.significant);
files = strings(0, 1);
if isempty(sigIdx)
    return
end
nPer = results.config.curvesPerPage;
nPages = ceil(numel(sigIdx) / nPer);
for page = 1:nPages
    idx = sigIdx((page-1)*nPer + 1:min(page*nPer, numel(sigIdx)));
    nCol = min(3, numel(idx));
    nRow = ceil(numel(idx) / nCol);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 300*nRow]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    tl = tiledlayout(nRow, nCol, 'Padding', 'compact', 'TileSpacing', 'compact');
    for j = 1:numel(idx)
        r = idx(j);
        ax = nexttile(tl);
        hold(ax, 'on');
        if results.config.includePointwiseBands
            se = results.sd_D(r, :) ./ sqrt(results.effectiveN(r));
            lo = results.mean_D(r, :) - 1.96 .* se;
            hi = results.mean_D(r, :) + 1.96 .* se;
            fill(ax, [results.timeSec; flipud(results.timeSec)], ...
                [lo'; flipud(hi')], [0.75 0.85 0.95], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.55);
        end
        plot(ax, results.timeSec, results.mean_D(r, :), 'b-', 'LineWidth', 1.7);
        yline(ax, 0, 'k-');
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Mean residual contrast');
        title(ax, sprintf('%s | p_Holm=%.3g', char(results.roiNames(r)), results.p_holm(r)), ...
            'Interpreter', 'none');
        grid(ax, 'on');
    end
    if results.config.includePointwiseBands
        title(tl, ['Significant ROI trajectories; shaded 95% pointwise descriptive bands ' ...
            '(not simultaneous or selection-adjusted)']);
    else
        title(tl, 'Mean residual contrast trajectories for significant omnibus ROIs');
    end
    outFile = fullfile(outputDir, sprintf('significant_roi_curves_page_%02d.png', page));
    if isfile(outFile)
        error('Refusing to overwrite existing result: %s', outFile);
    end
    exportgraphics(fig, outFile, 'Resolution', 180);
    files(end+1, 1) = string(outFile); %#ok<AGROW>
    clear cleanup
end
end

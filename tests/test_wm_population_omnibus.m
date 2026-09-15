function test_wm_population_omnibus()
%TEST_WM_POPULATION_OMNIBUS Focused numerical/statistical implementation checks.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
moduleDir = fullfile(repoRoot, 'code', 'residual_contrast', 'analysis', 'wm_population_omnibus');
addpath(moduleDir);

baseCfg = struct('nPerm', 127, 'seed', 7719, 'alpha', 0.05, ...
    'runMode', 'verification', 'inputDimensionOrder', 'subjects-rois-time', ...
    'saveOutputs', false, 'makePlots', false, 'verbose', false, ...
    'batchSize', 17, 'roiBatchSize', 2, ...
    'storePermutationStatistics', true, 'storePermutationSigns', true);
baseCfg.resampling = struct('scheme', 'independent-subject-signflip', ...
    'independentSubjectsVerified', false, 'nullSymmetryVerified', false, ...
    'verificationNote', 'synthetic verification');

%% Inferential mode refuses unresolved exchangeability assumptions.
cfgBlocked = baseCfg;
cfgBlocked.runMode = 'inference';
try
    temporal_omnibus_signflip(randn(5, 2, 3), (1:5)', [], [], (0:2)', cfgBlocked);
    error('Expected unresolved-resampling check did not fire.');
catch ME
    assert(strcmp(ME.identifier, 'temporal_omnibus_signflip:IndependenceUnresolved'));
end

%% Optimized statistics equal direct loops.
rng(1001, 'twister');
D = randn(9, 4, 5);
D(2, 2, 3) = NaN;
ids = (101:109)';
r = temporal_omnibus_signflip(D, ids, [], [], (0:4)' .* 0.72, baseCfg);
[meanDirect, sdDirect, tDirect, sDirect] = direct_observed(D);
assert(equal_with_nan_inf(r.mean_D, meanDirect, 1e-12));
assert(equal_with_nan_inf(r.sd_D, sdDirect, 1e-12));
assert(equal_with_nan_inf(r.T_obs, tDirect, 1e-11));
assert(equal_with_nan_inf(r.S_obs, sDirect, 1e-10));
for b = 1:baseCfg.nPerm
    sb = direct_permutation(D, double(r.permutationSigns(b, :))');
    assert(equal_with_nan_inf(r.S_perm(b, :)', sb, 1e-9), ...
        'Optimized permutation statistic differs from direct loop.');
end

%% Holm hand calculation.
p = [0.01; 0.04; 0.03; 0.20];
assert(max(abs(holm_adjust(p) - [0.04; 0.09; 0.09; 0.20])) < 1e-15);

%% Reversing A-B flips mean/T, but not S or permutation p-values.
rNeg = temporal_omnibus_signflip(-D, ids, [], [], (0:4)' .* 0.72, baseCfg);
assert(equal_with_nan_inf(rNeg.mean_D, -r.mean_D, 1e-12));
assert(equal_with_nan_inf(rNeg.T_obs, -r.T_obs, 1e-11));
assert(equal_with_nan_inf(rNeg.S_obs, r.S_obs, 1e-10));
assert(isequal(rNeg.p_raw, r.p_raw) && isequal(rNeg.p_holm, r.p_holm));

%% Missingness, untestable ROIs, and zero-variance conventions.
Dq = randn(6, 4, 3);
Dq(1, 1, 2) = NaN;
Dq(:, 2, :) = 0;
Dq(:, 3, :) = 2;
Dq(2:end, 4, :) = NaN;
q = temporal_omnibus_signflip(Dq, (1:6)', [], [], (0:2)', baseCfg);
assert(q.effectiveN(1) == 5 && q.effectiveN(4) == 1);
assert(all(q.T_obs(2, :) == 0));
assert(all(isinf(q.T_obs(3, :))) && all(q.T_obs(3, :) > 0));
assert(q.p_raw(4) == 1 && ~q.significant(4));
assert(contains(q.dataQuality.status(4), 'untestable'));
assert(numel(q.excludedSubjectIDs{1}) == 1 && q.excludedSubjectIDs{1} == "1");

%% Same seed is reproducible, including across a different batch size.
cfg2 = baseCfg;
cfg2.batchSize = 31;
rAgain = temporal_omnibus_signflip(D, ids, [], [], (0:4)' .* 0.72, cfg2);
assert(isequal(r.permutationSigns, rAgain.permutationSigns));
assert(equal_with_nan_inf(r.S_perm, rAgain.S_perm, 1e-10));
assert(isequal(r.p_raw, rAgain.p_raw));

%% Sanity checks: null ROI and injected monophasic/biphasic effects.
% This is a power sanity check, not evidence or proof of type-I error control.
rng(404, 'twister');
n = 80; t = 8;
noise = 0.35 .* randn(n, 3, t);
Ds = noise;
Ds(:, 2, 3:5) = Ds(:, 2, 3:5) + 1.25;
antiNoise = 0.15 .* randn(n, 1);
Ds(:, 3, 2) = antiNoise + 1.5;
Ds(:, 3, 7) = -antiNoise - 1.5;
cfgSanity = baseCfg;
cfgSanity.nPerm = 999;
cfgSanity.batchSize = 100;
cfgSanity.storePermutationSigns = false;
cfgSanity.storePermutationStatistics = false;
sanity = temporal_omnibus_signflip(Ds, (1:n)', [], [], (0:t-1)', cfgSanity);
assert(sanity.p_holm(2) <= 0.05, 'Injected monophasic effect was not detected.');
assert(sanity.p_holm(3) <= 0.05, 'Injected sign-changing effect was not detected.');
assert(abs(mean(sanity.mean_D(3, [2 7]))) < 1e-12, ...
    'Biphasic synthetic effect should cancel under signed temporal averaging.');
assert(sanity.S_obs(3) > sanity.S_obs(1), ...
    'Omnibus statistic should retain evidence from opposite-signed phases.');

fprintf('All WM population temporal-omnibus verification checks passed.\n');
end

function [meanD, sdD, tObs, sObs] = direct_observed(D)
[~, nROI, nTime] = size(D);
meanD = nan(nROI, nTime);
sdD = nan(nROI, nTime);
tObs = nan(nROI, nTime);
sObs = nan(nROI, 1);
for r = 1:nROI
    keep = all(isfinite(reshape(D(:, r, :), size(D, 1), nTime)), 2);
    x = reshape(D(keep, r, :), sum(keep), nTime);
    if size(x, 1) < 2, continue; end
    meanD(r, :) = mean(x, 1);
    sdD(r, :) = std(x, 0, 1);
    tObs(r, :) = t_convention(meanD(r, :), sdD(r, :), size(x, 1));
    sObs(r) = sum(tObs(r, :) .^ 2);
end
end

function s = direct_permutation(D, signs)
[~, nROI, nTime] = size(D);
s = nan(nROI, 1);
for r = 1:nROI
    keep = all(isfinite(reshape(D(:, r, :), size(D, 1), nTime)), 2);
    x = reshape(D(keep, r, :), sum(keep), nTime) .* signs(keep);
    if size(x, 1) < 2, continue; end
    m = mean(x, 1);
    sd = std(x, 0, 1);
    tt = t_convention(m, sd, size(x, 1));
    s(r) = sum(tt .^ 2);
end
end

function t = t_convention(m, sd, n)
t = nan(size(m));
ok = sd > 0;
t(ok) = m(ok) ./ (sd(ok) ./ sqrt(n));
t(sd == 0 & m == 0) = 0;
t(sd == 0 & m ~= 0) = sign(m(sd == 0 & m ~= 0)) .* Inf;
end

function tf = equal_with_nan_inf(a, b, tolerance)
sameNaN = isnan(a) == isnan(b);
sameInf = isinf(a) == isinf(b) & (~isinf(a) | sign(a) == sign(b));
finite = isfinite(a) & isfinite(b);
delta = zeros(size(a));
delta(finite) = abs(a(finite) - b(finite));
tf = all(sameNaN(:)) && all(sameInf(:)) && all(delta(finite) <= tolerance);
end

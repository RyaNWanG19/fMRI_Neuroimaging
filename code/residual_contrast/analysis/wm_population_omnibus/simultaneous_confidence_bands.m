function out = simultaneous_confidence_bands(D, primary)
%SIMULTANEOUS_CONFIDENCE_BANDS Joint bands on the observed time grid.
% D: canonical [subject x test x time]; primary: temporal_omnibus_signflip
% result with matching complete-case masks and summaries. Tests may be
% flattened ROI-category pairs. No omnibus-based selection is performed.
%
% Default: centered whole-subject bootstrap-t, recomputing sample sizes/SE
% in each bootstrap sample. Subjects are sampled jointly across all cells.
% Optional 'gaussian-multiplier' uses fixed observed standard errors:
% Z_b(r,t) = sum_i g_bi (D_i(r,t)-mean_D(r,t)) / sqrt(sum_i residual_i^2).
% One g_bi ~ N(0,1) is shared over ALL tests/times for subject i. Conditional
% covariance equals the sample correlation of the mean estimators (including
% partially overlapping complete-case masks). Maxima cover both directions.
% Bands are mean_D +/- q*SE, q an upper empirical quantile of max(abs(Z_b)).
% Approximate inference for independent subject vectors and suitable
% moments/nondegenerate variances; NOT an exact sign-flip confidence set.
% Does not require joint null symmetry. Missingness must permit inference
% for the intended population. No smoothing or between-grid coverage claim.
%
% Optional primary.config.scb fields: method ('subject-bootstrap-t'),
% nBootstrap (9999), seed (primary seed
% + 32452843), batchSize (primary batchSize), testBatchSize (64),
% minSubjects (10; practical guard, not a guarantee of valid coverage).
% Ineligible cells get [-Inf,Inf], never a spurious zero-width interval.

cfg = primary.config;
if ~isfield(cfg,'scb'), cfg.scb = struct(); end
scb = cfg.scb;
scb = default(scb,'method','subject-bootstrap-t');
scb = default(scb,'nBootstrap',9999);
scb = default(scb,'seed',mod(double(cfg.seed)+32452843,2^32));
scb = default(scb,'batchSize',cfg.batchSize);
scb = default(scb,'testBatchSize',64);
scb = default(scb,'minSubjects',10);
validateattributes(scb.nBootstrap,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(scb.seed,{'numeric'},{'scalar','integer','nonnegative','<',2^32});
validateattributes(scb.batchSize,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(scb.testBatchSize,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(scb.minSubjects,{'numeric'},{'scalar','integer','>=',2,'finite'});
validateattributes(cfg.alpha,{'numeric'},{'scalar','>',0,'<',1,'finite'});
assert(any(strcmp(scb.method,{'subject-bootstrap-t','gaussian-multiplier'})), ...
    'scb.method must be subject-bootstrap-t or gaussian-multiplier.');
assert(isnumeric(D) && isreal(D),'D must be a real numeric array.');
[ns,m,nt] = size(D);
assert(numel(D)==ns*m*nt && size(primary.mean_D,1)==m && ...
    size(primary.mean_D,2)==nt && numel(primary.subjectIDs)==ns, ...
    'D must match the primary subject/test/time axes.');
time = primary.timeSec(:);
assert(numel(time)==nt && all(isfinite(time)) && all(diff(time)>0), ...
    'SCB time grid must be finite and strictly increasing.');
mask = all(isfinite(D),3);
assert(isequal(mask,primary.subjectMasks),'Complete-case masks do not match.');
if strcmpi(cfg.runMode,'inference')
    if ~isequal(cfg.resampling.independentSubjectsVerified,true) || ...
            isempty(strtrim(char(string(cfg.resampling.verificationNote))))
        error('simultaneous_confidence_bands:IndependenceUnresolved', ...
            'Inferential bands require documented independent subject rows.');
    end
end

% Each test uses its fixed complete-case subjects. Zero-padding centers
% preserves subject identity and overlap across every test in the family.
normalized = zeros(ns,m,nt);
se = nan(m,nt);
eligible = false(m,nt);
status = repmat("insufficient_complete_subjects",m,nt);
for r = 1:m
    n = sum(mask(:,r));
    if n < 2, continue; end
    x = double(reshape(D(mask(:,r),r,:),n,nt));
    y = x-mean(x,1);
    normY = vecnorm(y,2,1);
    se(r,:) = normY/sqrt(n*(n-1));
    if n < scb.minSubjects, continue; end
    active = normY>0 & isfinite(normY) & se(r,:)>0 & isfinite(se(r,:)) & ...
        isfinite(primary.mean_D(r,:));
    eligible(r,:) = active;
    status(r,:) = "zero_or_nonfinite_variance_unbounded";
    status(r,active) = "approximate_simultaneous_band";
    normalized(mask(:,r),r,active) = reshape(y(:,active)./normY(active),n,1,[]);
end

maxima = zeros(scb.nBootstrap,1);
q = Inf;
if any(eligible,'all')
    stream = RandStream('mt19937ar','Seed',scb.seed);
    stream.NormalTransform = 'Inversion'; % stable RNG across batch sizes
    for first = 1:scb.batchSize:scb.nBootstrap
        nb = min(scb.batchSize,scb.nBootstrap-first+1);
        if strcmp(scb.method,'gaussian-multiplier')
            weights = randn(stream,ns,nb)';
        else
            draws = randi(stream,ns,ns,nb);
            weights = full(sparse(draws(:),repelem((1:nb)',ns),1,ns,nb))';
        end
        batchMax = zeros(nb,1);
        for start = 1:scb.testBatchSize:m
            idx = start:min(m,start+scb.testBatchSize-1);
            e = reshape(normalized(:,idx,:),ns,[]);
            z = weights*e;
            if strcmp(scb.method,'subject-bootstrap-t')
                counts = repmat(weights*double(mask(:,idx)),1,nt);
                squares = weights*(e.^2);
                varianceNumerator = counts.*squares-z.^2;
                tolerance = 64*eps(max(counts.*squares,1));
                if any(varianceNumerator < -tolerance,'all')
                    error('simultaneous_confidence_bands:NegativeVariance', ...
                        'Substantively negative bootstrap variance.');
                end
                % Degenerate resampled eligible cells count as Inf, rather
                % than silently disappearing from the joint maximum.
                valid = counts>=2 & varianceNumerator>tolerance;
                eligibleBatch = reshape(eligible(idx,:),1,[]);
                invalid = ~valid & eligibleBatch;
                z(valid) = z(valid).*sqrt(counts(valid)-1)./sqrt(varianceNumerator(valid));
                z(invalid) = Inf;
                z(:,~eligibleBatch) = 0;
            end
            batchMax = max(batchMax,max(abs(z),[],2));
        end
        maxima(first:first+nb-1) = batchMax;
        if cfg.verbose && (first==1 || first+nb-1==scb.nBootstrap)
            fprintf('  simultaneous bands: %d/%d bootstrap replicates\n', ...
                first+nb-1,scb.nBootstrap);
        end
    end
    % Upper order statistic, no interpolation. If B is too small for
    % requested alpha, keep q=Inf rather than overstate MC resolution.
    rank = ceil((scb.nBootstrap+1)*(1-cfg.alpha));
    if rank<=scb.nBootstrap
        ordered = sort(maxima);
        q = ordered(rank);
    end
end
lower = -inf(m,nt); upper = inf(m,nt);
if isfinite(q)
    lower(eligible) = primary.mean_D(eligible)-q*se(eligible);
    upper(eligible) = primary.mean_D(eligible)+q*se(eligible);
elseif any(eligible,'all')
    status(eligible) = "bootstrap_critical_value_unavailable_unbounded";
end
positive = eligible & lower>0;
negative = eligible & upper<0;

out = struct();
out.method = scb.method;
out.config = scb;
out.alpha = cfg.alpha;
out.criticalValue = q;
out.bootstrapMaxima = maxima;
out.degenerateBootstrapReplicates = nnz(~isfinite(maxima));
out.mean = primary.mean_D;
out.standardError = se;
out.lower = lower; out.upper = upper;
out.eligible = eligible; out.cellStatus = status;
out.positive = positive; out.negative = negative;
out.timeSec = time;
out.plannedTestCount = m;
out.plannedCellCount = m*nt;
out.calibratedCellCount = nnz(eligible);
out.correctionFamily = sprintf(['Joint two-sided bands across all %d planned tests x %d sampled lags; ' ...
    'ineligible cells unbounded. HRF model %s; no correction across models.'], ...
    m,nt,cfg.hrfModelName);
out.inferenceValid = strcmpi(cfg.runMode,'inference') && isfinite(q) && any(eligible,'all');
out.inferenceStatus = "non-inferential_"+string(cfg.runMode);
if out.inferenceValid
    out.inferenceStatus = "approximate_conditional_on_independent_subjects_and_regular_mean_estimation";
end
out.assumptions = ['Independent subject vectors; suitable finite moments and nondegenerate variances; ' ...
    'adequate sample size for joint bootstrap approximation; fixed source missingness masks and valid population ' ...
    'interpretation after exclusions. No proof of independence/coverage from this flag.'];
out.coverageScope = 'Observed lags only; no smoothing, interpolation or exact onset/duration inference.';
out.monteCarloNote = 'Finite-bootstrap upper empirical quantile; repeat/increase B for borderline band crossings.';

% Test-major long table preserves all planned tests and ineligible cells.
testIndex = repelem((1:m)',nt);
lagIndex = repmat((1:nt)',m,1);
testIndex = testIndex(:); lagIndex = lagIndex(:);
out.table = table(testIndex,reshape(primary.roiIDs(testIndex),[],1), ...
    reshape(primary.roiNames(testIndex),[],1),reshape(primary.effectiveN(testIndex),[],1), ...
    lagIndex,reshape(time(lagIndex),[],1),reshape(out.mean',[],1), ...
    reshape(se',[],1),reshape(lower',[],1),reshape(upper',[],1), ...
    reshape(eligible',[],1),reshape(positive',[],1),reshape(negative',[],1), ...
    reshape(status',[],1),repmat(out.inferenceStatus,m*nt,1), ...
    'VariableNames',{'test_index','roi_id','roi_name','effective_n','lag_index', ...
    'time_sec','mean_difference','standard_error','lower_simultaneous','upper_simultaneous', ...
    'band_eligible','supported_positive','supported_negative','cell_status','inference_status'});
end

function s = default(s,name,value)
if ~isfield(s,name) || isempty(s.(name)), s.(name)=value; end
end

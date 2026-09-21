function out = omnibus_sensitivity(D, primary)
%OMNIBUS_SENSITIVITY Descriptive diagnostics and centered curve bootstrap.
% D is double [subjects x ROIs x time]. Fixed primary complete-case masks
% are reused. Each ROI resamples n centered whole curves with replacement.
% This approximates the zero-mean null without imposing sign symmetry.
% Requires iid retained curves within ROI, finite suitable moments and a
% nondegenerate sampling distribution; NOT an exact randomization test.
% ROI streams are independent: no joint/cross-ROI bootstrap is claimed.
% Holm still includes every planned ROI. No data-driven subject exclusions.

cfg = primary.config;
[~, m, nt] = size(D);
B = cfg.nBootstrap;
p = ones(m,1);
skew = nan(m,nt);
maxZ = nan(m,1);
influence = nan(m,1);
influentialID = strings(m,1);
status = repmat("untestable",m,1);
bootStats = [];
if cfg.storePermutationStatistics, bootStats = nan(B,m); end
for r = 1:m
    keep = primary.subjectMasks(:,r);
    x = reshape(D(keep,r,:),sum(keep),nt);
    n = size(x,1);
    if n < 2, continue; end
    mu = mean(x,1);
    sd = std(x,0,1);
    active = sd > 0;
    y = x - mu;
    if any(active)
        z = y(:,active)./sd(active);
        v = mean(y(:,active).^2,1);
        skew(r,active) = mean(y(:,active).^3,1)./(v.^1.5);
        maxZ(r) = max(abs(z),[],'all');
        % Euclidean displacement of leave-one-out mean in full-sample SE
        % units. Descriptive influence, not a significance/exclusion rule.
        shifts = sqrt(sum((z.*sqrt(n)./(n-1)).^2,2));
        [influence(r),idx] = max(shifts);
        ids = primary.subjectIDs(keep);
        influentialID(r) = ids(idx);
    end
    if any(sd == 0 & mu ~= 0)
        status(r) = "degenerate_nonzero_constant_curve_cell_bootstrap_unreliable";
        continue
    end
    status(r) = "approximate_centered_bootstrap";
    y = y(:,active);
    stream = RandStream('mt19937ar','Seed',mod(double(cfg.seed)+104729+r,2^32));
    exceed = 0;
    for first = 1:cfg.batchSize:B
        nb = min(cfg.batchSize,B-first+1);
        draws = randi(stream,n,n,nb); % each column is one whole-curve sample
        weights = full(sparse(draws(:),repelem((1:nb)',n),1,n,nb))';
        sums = weights*y;
        squares = weights*(y.^2);
        numer = n.*squares-sums.^2;
        tol = 64.*eps(max(n.*squares,1));
        if any(numer < -tol,'all')
            error('omnibus_sensitivity:NegativeVariance','Substantive negative bootstrap variance.');
        end
        tsq = (n-1).*sums.^2./numer;
        % Recompute near cancellation with direct centered sample SD.
        suspect = abs(numer) <= tol;
        [bb,tt] = find(suspect);
        for j = 1:numel(bb)
            values = y(draws(:,bb(j)),tt(j));
            ss = std(values,0);
            mm = mean(values);
            if ss == 0
                tsq(bb(j),tt(j)) = 0;
                if mm ~= 0, tsq(bb(j),tt(j)) = Inf; end
            else
                tsq(bb(j),tt(j)) = (mm/(ss/sqrt(n)))^2;
            end
        end
        sb = sum(tsq,2);
        exceed = exceed + sum(sb >= primary.S_obs(r));
        if cfg.storePermutationStatistics, bootStats(first:first+nb-1,r) = sb; end
    end
    p(r) = (1+exceed)/(B+1);
    if cfg.verbose && (mod(r,25)==0 || r==m)
        fprintf('  centered bootstrap completed ROI %d/%d\n',r,m);
    end
end
adj = holm_adjust(p);
sig = adj <= cfg.alpha & status == "approximate_centered_bootstrap";
maxSkew = max(abs(skew),[],2,'omitnan');
out.table = table(primary.roiIDs,primary.roiNames,primary.effectiveN, ...
    maxSkew,maxZ,influence,influentialID,primary.p_holm,p,adj,sig, ...
    sig~=primary.significant,status, ...
    'VariableNames',{'roi_id','roi_name','effective_n','max_abs_centered_skewness', ...
    'max_abs_standardized_residual','max_LOO_mean_shift_SE_norm','influential_row_id', ...
    'signflip_p_holm','bootstrap_p_raw','bootstrap_p_holm','bootstrap_significant', ...
    'decision_disagreement','bootstrap_status'});
out.centeredSkewnessByTime = skew;
out.S_boot = bootStats;
out.nBootstrap = B;
out.minimumAttainableP = 1/(B+1);
out.note = ['Approximate centered whole-curve iid bootstrap; independently sampled per ROI. ', ...
    'Finite moments and nondegeneracy required. Diagnostics neither prove symmetry nor ', ...
    'trigger exclusions. Disagreement can reflect assumptions or Monte Carlo noise. ', ...
    'Sign-flip remains the prespecified primary analysis; do not select whichever p is smaller.'];
out.rng = 'mt19937ar; ROI r seed mod(primary seed + 104729 + r,2^32)';
end

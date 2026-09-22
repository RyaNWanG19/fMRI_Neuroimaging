function fit = fit_signed_cp(X, rankK, cfg)
%FIT_SIGNED_CP Unconstrained three-way CP-ALS, using base MATLAB only.
% X is [region x time x subject]. A and B have unit-length columns;
% C carries signed subject amplitudes: X(r,t,s) ~ sum_k A(r,k)B(t,k)C(s,k).
% No centering, standardization, nonnegativity, or independence is imposed.
% Multiple random starts mitigate local optima; convergence is not uniqueness.
if nargin < 3, cfg = struct(); end
cfg = defaults(cfg, 'nStarts', 10, 'maxIter', 300, 'tol', 1e-7, 'seed', 20260921);
validateattributes(X, {'numeric'}, {'real','finite','nonempty'});
assert(ndims(X) <= 3, 'X must have region-time-subject order.');
validateattributes(rankK, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.nStarts, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.maxIter, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.tol, {'numeric'}, {'scalar','positive','finite'});
validateattributes(cfg.seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
oldRng = rng; cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(cfg.seed, 'twister');
X = double(X); [nr, nt, ns] = size(X);
energy = sum(X(:).^2);
assert(energy > 0 && isfinite(energy), 'X must have positive finite energy.');
% Global numerical rescaling preserves all relative amplitudes.
scale = sqrt(energy); X = X / scale;
X1 = reshape(X, nr, []);
X2 = reshape(permute(X, [2 1 3]), nt, []);
X3 = reshape(permute(X, [3 1 2]), ns, []);
best = Inf;
err = nan(cfg.nStarts,1); iterations = zeros(cfg.nStarts,1);
converged = false(cfg.nStarts,1); histories = cell(cfg.nStarts,1);
for start = 1:cfg.nStarts
    A = unitcols(randn(nr,rankK)); B = unitcols(randn(nt,rankK));
    C = randn(ns,rankK); previous = Inf;
    history = nan(cfg.maxIter,1);
    for it = 1:cfg.maxIter
        A = (X1 * khatri_rao(C,B)) * pinv((C'*C).*(B'*B));
        A = unitcols(A);
        B = (X2 * khatri_rao(C,A)) * pinv((C'*C).*(A'*A));
        B = unitcols(B);
        Z = khatri_rao(B,A);
        C = (X3 * Z) * pinv((B'*B).*(A'*A));
        % Direct residual avoids cancellation near an exact reconstruction.
        residual = X3 - C*Z';
        current = sum(residual(:).^2);
        history(it) = current;
        if ~isfinite(current), break; end
        if abs(previous-current) <= cfg.tol * max(previous,eps) && isfinite(previous)
            converged(start) = true; break;
        end
        previous = current;
    end
    err(start) = current; iterations(start) = it; histories{start} = history(1:it);
    if isfinite(current) && current < best
        best = current; bestA = A; bestB = B; bestC = C*scale; bestStart = start;
    end
end
assert(isfinite(best), 'All CP starts failed.');
A = bestA; B = bestB; C = bestC;
% Deterministic signs and order facilitate plotting; matching is still needed
% across fits. Components are NOT orthogonal, so their energies do not add.
for k = 1:rankK
    [~,i] = max(abs(A(:,k))); sg = sign(A(i,k)); if sg == 0, sg = 1; end
    A(:,k) = A(:,k)*sg; C(:,k) = C(:,k)*sg;
    [~,i] = max(abs(B(:,k))); sg = sign(B(i,k)); if sg == 0, sg = 1; end
    B(:,k) = B(:,k)*sg; C(:,k) = C(:,k)*sg;
end
[amplitudes, order] = sort(sqrt(sum(C.^2,1)), 'descend');
fit = struct('A',A(:,order),'B',B(:,order),'C',C(:,order), ...
    'rank',rankK,'relativeError',sqrt(best),'energyFit',1-best, ...
    'componentNorm',amplitudes,'bestStart',bestStart, ...
    'converged',converged(bestStart),'iterations',iterations(bestStart), ...
    'startRelativeError',sqrt(err),'startConverged',converged, ...
    'startIterations',iterations,'objectiveHistory',{histories},'config',cfg);
fit.cancellationRatio = sum(amplitudes) / max(sqrt(sum((C*khatri_rao(B,A)').^2,'all')),eps);
fit.possibleDegeneracy = fit.cancellationRatio > 10;
end

function Z = khatri_rao(U,V)
Z = reshape(reshape(V,size(V,1),1,[]) .* reshape(U,1,size(U,1),[]), [],size(U,2));
end
function A = unitcols(A)
A = A ./ max(sqrt(sum(A.^2,1)),eps);
end
function cfg = defaults(cfg,varargin)
for i=1:2:numel(varargin)
    if ~isfield(cfg,varargin{i}), cfg.(varargin{i})=varargin{i+1}; end
end
end

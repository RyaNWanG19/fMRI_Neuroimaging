function [loss, baselineAccuracy] = svm_paired_region_reliance(mdl, X, nRepeats, seed)
% Held-out model reliance for a linear binary SVM with paired rows.
% X = [condition-positive subjects; condition-negative SAME subjects],
% preprocessed using training statistics only. For each region and repeat,
% independently swap the two values with probability 0.5 within each subject.
% Return baseline minus mean perturbed accuracy (fraction, not percentage).
% This retains subject pairing and regional marginal values, but disrupts
% cross-region dependence. It is model reliance, not a regional hypothesis test.

validateattributes(X, {'numeric'}, {'2d','nonempty','finite'});
validateattributes(nRepeats, {'numeric'}, {'scalar','integer','positive'});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative','<=',2^32-1});
assert(mod(size(X,1),2)==0, 'X must contain two equal paired condition blocks.');
assert(isequal(mdl.ClassNames(:), [-1; 1]), 'Expected ClassNames [-1; 1].');
assert(strcmp(mdl.KernelParameters.Function,'linear'), 'Requires a linear SVM.');
assert(ischar(mdl.ScoreTransform) && strcmp(mdl.ScoreTransform,'none'), ...
    'Requires untransformed linear scores.');
w = mdl.Beta(:)';
assert(numel(w)==size(X,2), 'Feature count does not match model.');
n = size(X,1)/2;
[labels, scores] = predict(mdl, X);
baselineAccuracy = mean(0.5*((labels(1:n)==1)+(labels(n+1:end)==-1)));
% Positive-class score changes exactly by beta_j * change_in_feature_j.
% MATLAB resolves zero-score ties to the first class (-1).
delta = (X(n+1:end,:) - X(1:n,:)) .* w;
stream = RandStream('mt19937ar','Seed',seed);
perturbedAccuracy = zeros(1,size(X,2));
for q = 1:nRepeats
    changes = delta .* (rand(stream,n,size(X,2)) < 0.5);
    positive = scores(1:n,2) + changes;
    negative = scores(n+1:end,2) - changes;
    perturbedAccuracy = perturbedAccuracy + ...
        mean(0.5*((positive>0)+(negative<=0)),1);
end
loss = (baselineAccuracy - perturbedAccuracy/nRepeats)';
end

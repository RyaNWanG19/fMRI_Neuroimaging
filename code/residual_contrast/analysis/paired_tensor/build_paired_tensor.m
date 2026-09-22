function [D, audit] = build_paired_tensor(Results, weights, Flags, cfg)
%BUILD_PAIRED_TENSOR Form a linear contrast of [ROI x time x subject x condition].
% Complete subjects only: any nonfinite source entry in a used condition
% excludes the subject. Flags, when present, must be [condition x subject].
% Entirely zero source conditions are ambiguous without Flags; default excludes
% these subjects, with explicit audit. Zero CONTRASTS are valid and retained.
if nargin < 3, Flags = []; end
if nargin < 4, cfg = struct(); end
if ~isfield(cfg,'excludeAmbiguousZero'), cfg.excludeAmbiguousZero = true; end
validateattributes(Results,{'numeric'},{'real','nonempty'});
assert(ndims(Results)<=4, 'Expected ROI-time-subject-condition order.');
weights = double(weights(:));
assert(numel(weights)==size(Results,4) && all(isfinite(weights)) && any(weights~=0), ...
    'One finite weight per condition is required, with a nonzero contrast.');
assert(abs(sum(weights)) < 1e-10*sum(abs(weights)), 'Paired contrast weights must sum to zero.');
ns = size(Results,3); used = find(weights~=0);
badFinite = false(ns,1); badFlag = false(ns,1); zeroSource = false(ns,1);
if ~isempty(Flags)
    assert(isequal(size(Flags),[size(Results,4), ns]), 'Flags must be condition x subject.');
    badFlag = any(~isfinite(Flags(used,:)) | Flags(used,:)~=1,1)';
end
D = zeros(size(Results,1),size(Results,2),ns);
for c = used'
    X = double(Results(:,:,:,c));
    flat = reshape(X,[],ns);
    badFinite = badFinite | any(~isfinite(flat),1)';
    zeroSource = zeroSource | all(flat==0,1)';
    D = D + weights(c)*X;
end
ambiguousZero = zeroSource & isempty(Flags);
keep = ~(badFinite | badFlag | (ambiguousZero & cfg.excludeAmbiguousZero));
audit = table((1:ns)',keep,badFinite,badFlag,zeroSource,ambiguousZero, ...
    'VariableNames',{'sourceRow','retained','nonfiniteSource','invalidFlag', ...
    'allZeroSourceCondition','ambiguousZeroWithoutFlags'});
D = D(:,:,keep);
assert(size(D,3)>=2, 'Fewer than two complete subjects remain.');
assert(any(D(:)~=0), 'The retained contrast is entirely zero.');
end

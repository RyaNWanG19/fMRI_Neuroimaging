function cfg = svm_rfe_config(cfg)
% Defaults for nested, subject-grouped linear SVM-RFE.
if nargin==0 || isempty(cfg), cfg=struct(); end
assert(isstruct(cfg) && isscalar(cfg), 'featureSelection must be a scalar struct.');
defaults=struct('innerFolds',5,'featureCounts',[10 25 50 100 200 Inf], ...
    'boxConstraint',1,'seed',31415);
names=fieldnames(defaults);
for j=1:numel(names)
    if ~isfield(cfg,names{j}), cfg.(names{j})=defaults.(names{j}); end
end
validateattributes(cfg.innerFolds,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cfg.boxConstraint,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.seed,{'numeric'},{'scalar','integer','nonnegative','<=',2^32-1});
v=cfg.featureCounts;
assert(isnumeric(v) && isvector(v) && ~isempty(v) && ...
    all((isfinite(v) & v>=1 & v==floor(v)) | v==Inf), ...
    'featureCounts must contain positive integers or Inf (all eligible regions).');
cfg.featureCounts=unique(double(v(:)'));
end

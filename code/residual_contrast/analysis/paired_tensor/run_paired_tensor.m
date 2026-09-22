function results = run_paired_tensor(cfg)
%RUN_PAIRED_TENSOR Descriptive Option B pilot with CP rank sweep and SVD baseline.
% Default: WM (2back body-0back body)-(2back faces-0back faces), cHRF.
% Custom data: cfg.dataFile, cfg.contrastWeights, cfg.contrastName; input must
% be Results [ROI x time x subject x condition]. No dimension inference.
% No inferential p-values or automatic selection of a biological rank.
if nargin < 1, cfg = struct(); end
root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
cfg = defaults(cfg,'repoRoot',root,'hrfModelName','cHRF', ...
    'contrastName','Body_LoadDiff_vs_Face_LoadDiff','ranks',1:6, ...
    'nStarts',10,'maxIter',300,'tol',1e-7,'seed',20260921, ...
    'nSplitHalf',2,'TR',0.72,'saveOutputs',true,'makePlots',true, ...
    'excludeAmbiguousZero',true);
customFile = isfield(cfg,'dataFile');
if ~customFile
    names = {'cHRF','cHRFderiv','sHRF'};
    i = find(strcmpi(cfg.hrfModelName,names),1);
    assert(~isempty(i),'Unknown HRF model.'); cfg.hrfModelName = names{i};
    cfg.dataFile = fullfile(cfg.repoRoot,'data','task_residual',['WM' names{i} '.mat']);
end
if ~isfield(cfg,'contrastWeights')
    assert(~customFile, 'Custom dataFile requires explicit contrastWeights.');
    assert(strcmp(cfg.contrastName,'Body_LoadDiff_vs_Face_LoadDiff'), ...
        'Other contrasts require explicit contrastWeights.');
    cfg.contrastWeights = [-1 1 0 0 1 -1 0 0];
end
validateattributes(cfg.ranks,{'numeric'},{'vector','integer','positive','finite'});
assert(numel(unique(cfg.ranks))==numel(cfg.ranks),'Ranks must be unique.');
validateattributes(cfg.nSplitHalf,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.TR,{'numeric'},{'scalar','positive','finite'});
if cfg.nSplitHalf>0, assert(max(cfg.ranks)<=8,'Exact split-half matching supports ranks <= 8.'); end
variables = whos('-file',cfg.dataFile);
toLoad = intersect({'Results','Flags'},{variables.name});
payload = load(cfg.dataFile,toLoad{:});
assert(isfield(payload,'Results'),'MAT file must contain Results.');
if ~customFile
    assert(size(payload.Results,1)==489 && size(payload.Results,2)==41 && size(payload.Results,4)==8, ...
        'Default WM data must be 489 x 41 x subject x 8.');
end
flags = []; if isfield(payload,'Flags'), flags = payload.Flags; end
[D,audit] = build_paired_tensor(payload.Results,cfg.contrastWeights,flags,cfg);
sourceSize = size(payload.Results); clear payload
if ~isfield(cfg,'subjectIDs')
    ids = string(compose('row_%03d',audit.sourceRow));
    idSource = 'Ordinal source rows only; not verified cross-file subject identifiers';
else
    ids = string(cfg.subjectIDs(:)); idSource = 'cfg.subjectIDs';
end
assert(numel(ids)==height(audit) && numel(unique(ids))==numel(ids) && ...
    all(~ismissing(ids)) && all(strlength(ids)>0), 'Subject IDs must be nonempty and unique, one per source row.');
audit.subjectID = ids;
if ~isfield(cfg,'timeSec'), cfg.timeSec = (0:size(D,2)-1)'*cfg.TR; end
assert(numel(cfg.timeSec)==size(D,2) && all(isfinite(cfg.timeSec)) && all(diff(cfg.timeSec)>0), ...
    'timeSec must be finite, increasing, and match the time dimension.');
if ~isfield(cfg,'roiNames'), cfg.roiNames = compose('Parcel_%03d',(1:size(D,1))'); end
assert(numel(cfg.roiNames)==size(D,1),'roiNames must match the ROI dimension.');
assert(max(cfg.ranks)<=min(size(D,1)*size(D,2),size(D,3)), 'Ranks exceed matrix baseline dimensions.');
if cfg.nSplitHalf>0, assert(size(D,3)>=4,'Split-half requires at least four subjects.'); end
fprintf('Option B: %s / %s; %d of %d subjects retained; %d ambiguous zero-source subjects.\n', ...
    cfg.contrastName,cfg.hrfModelName,size(D,3),height(audit),sum(audit.ambiguousZeroWithoutFlags));
% Uncentered subject x (ROI*time) SVD: same target and units as CP.
% Eigenvalues of the smaller subject Gram matrix avoid a large full SVD.
M = reshape(D,[],size(D,3)); energy = sum(M(:).^2);
eigvals = sort(max(real(eig(M'*M)),0),'descend');
baselineFit = cumsum(eigvals)/energy;
oldRng = rng; cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(cfg.seed,'twister');
splits = cell(cfg.nSplitHalf,2);
for b=1:cfg.nSplitHalf
    p = randperm(size(D,3)); mid = floor(numel(p)/2);
    splits{b,1}=p(1:mid); splits{b,2}=p(mid+1:end);
end
n = numel(cfg.ranks); fits = cell(n,1); stability = cell(n,1);
cpFit = nan(n,1); svdFit = nan(n,1); converged = false(n,1);
degenerate = false(n,1); splitMedian = nan(n,1);
for j=1:n
    k = cfg.ranks(j); fitCfg = cfg; fitCfg.seed = cfg.seed+1000*k;
    fprintf('Fitting rank %d (%d starts)...\n',k,cfg.nStarts);
    fits{j} = fit_signed_cp(D,k,fitCfg);
    cpFit(j)=fits{j}.energyFit; svdFit(j)=baselineFit(k);
    converged(j)=fits{j}.converged; degenerate(j)=fits{j}.possibleDegeneracy;
    scores = nan(cfg.nSplitHalf,k); permutations = zeros(cfg.nSplitHalf,k);
    splitConverged = false(cfg.nSplitHalf,2);
    splitDegenerate = false(cfg.nSplitHalf,2);
    for b=1:cfg.nSplitHalf
        fitCfg.seed = cfg.seed+1000*k+2*b;
        left = fit_signed_cp(D(:,:,splits{b,1}),k,fitCfg);
        fitCfg.seed = fitCfg.seed+1;
        right = fit_signed_cp(D(:,:,splits{b,2}),k,fitCfg);
        similarity = abs(left.A'*right.A).*abs(left.B'*right.B);
        % Exact one-to-one matching of absolute spatial-temporal congruence.
        candidates = perms(1:k); total = zeros(size(candidates,1),1);
        for a=1:k, total=total+similarity(a,candidates(:,a))'; end
        [~,best]=max(total); match=candidates(best,:);
        for a=1:k, scores(b,a)=similarity(a,match(a)); end
        permutations(b,:)=match;
        splitConverged(b,:)=[left.converged,right.converged];
        splitDegenerate(b,:)=[left.possibleDegeneracy,right.possibleDegeneracy];
    end
    stability{j}=struct('matchedCongruence',scores,'permutations',permutations, ...
        'converged',splitConverged,'possibleDegeneracy',splitDegenerate);
    if ~isempty(scores), splitMedian(j)=median(scores(:)); end
    fprintf('  energy fit %.4f; SVD %.4f; converged %d; split congruence %.4f\n', ...
        cpFit(j),svdFit(j),converged(j),splitMedian(j));
end
summary = table(cfg.ranks(:),cpFit,svdFit,converged,degenerate,splitMedian, ...
    'VariableNames',{'rank','cpEnergyFit','svdEnergyFit','bestStartConverged', ...
    'possibleDegeneracy','medianSplitHalfCongruence'});
results=struct('config',cfg,'fits',{fits},'rankSummary',summary, ...
    'stability',{stability},'splitSubjectIndices',{splits},'subjectAudit',audit, ...
    'retainedSubjectIDs',ids(audit.retained),'sourceSize',sourceSize, ...
    'subjectIDSource',idSource,'flagsAvailable',~isempty(flags), ...
    'meanContrast',mean(D,3),'tensorSize',size(D), ...
    'interpretation','Descriptive, uncentered signed CP; no inference or automatic rank selection');
if cfg.saveOutputs
    if ~isfield(cfg,'outputDir')
        cfg.outputDir=fullfile(cfg.repoRoot,'data','task_residual','paired_tensor', ...
            regexprep(cfg.contrastName,'[^A-Za-z0-9_-]','_'),cfg.hrfModelName,datestr(now,'yyyymmddTHHMMSSFFF'));
    end
    if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
    results.config=cfg;
    save(fullfile(cfg.outputDir,'paired_tensor_results.mat'),'results','-v7.3');
    writetable(summary,fullfile(cfg.outputDir,'rank_summary.csv'));
    writetable(audit,fullfile(cfg.outputDir,'subject_audit.csv'));
    for j=1:n
        f=fits{j}; names=cellstr(compose('component_%02d',1:f.rank));
        spatial=array2table(f.A,'VariableNames',names); spatial.roiName=string(cfg.roiNames(:));
        temporal=array2table(f.B,'VariableNames',names); temporal.timeSec=cfg.timeSec(:);
        subject=array2table(f.C,'VariableNames',names); subject.subjectID=ids(audit.retained);
        writetable(spatial,fullfile(cfg.outputDir,sprintf('rank_%02d_spatial.csv',f.rank)));
        writetable(temporal,fullfile(cfg.outputDir,sprintf('rank_%02d_temporal.csv',f.rank)));
        writetable(subject,fullfile(cfg.outputDir,sprintf('rank_%02d_subject.csv',f.rank)));
    end
    if cfg.makePlots, plot_paired_tensor(results); end
    fprintf('Saved to %s\n',cfg.outputDir);
end
end
function cfg=defaults(cfg,varargin)
for i=1:2:numel(varargin)
    if ~isfield(cfg,varargin{i}), cfg.(varargin{i})=varargin{i+1}; end
end
end

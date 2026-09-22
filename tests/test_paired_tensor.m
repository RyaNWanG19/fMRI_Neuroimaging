function test_paired_tensor()
% Numerical recovery, contrast polarity/missingness, reproducibility and IO.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','analysis','paired_tensor'));
rng(71); A=randn(9,2); B=randn(7,2); C=randn(18,2);
X=zeros(9,7,18);
for s=1:18, X(:,:,s)=A*diag(C(s,:))*B'; end
cfg=struct('nStarts',4,'maxIter',500,'tol',1e-9,'seed',23);
before=rng; f=fit_signed_cp(X,2,cfg); after=rng;
assert(isequal(before,after),'Solver changed caller RNG.');
assert(f.relativeError<1e-5,'Failed to recover exact signed rank-two tensor.');
Y=zeros(size(X));
for s=1:18, Y(:,:,s)=f.A*diag(f.C(s,:))*f.B'; end
assert(sqrt(sum((X(:)-Y(:)).^2)/sum(X(:).^2))<1e-5);
assert(max(abs(sum(f.A.^2,1)-1))<1e-12 && max(abs(sum(f.B.^2,1)-1))<1e-12);
again=fit_signed_cp(X,2,cfg); assert(isequal(f.C,again.C));
negative=fit_signed_cp(-X,2,cfg);
assert(abs(f.energyFit-negative.energyFit)<1e-9,'Contrast reversal changed fit.');
single=fit_signed_cp(reshape((1:9)'*(1:7),9,7,1),1,cfg);
assert(single.relativeError<1e-10);
Results=randn(9,7,18,4);
Results(:,:,:,1)=X+Results(:,:,:,2);
[D,audit]=build_paired_tensor(Results,[1 -1 0 0]);
assert(max(abs(D(:)-X(:)))<1e-12 && all(audit.retained));
Results(:,:,2,1)=NaN; Results(:,:,3,1)=0;
% A true zero difference is retained; unused missing conditions are irrelevant.
Results(:,:,4,1)=Results(:,:,4,2); Results(:,:,:,4)=NaN;
[D,audit]=build_paired_tensor(Results,[1 -1 0 0]);
assert(~audit.retained(2) && ~audit.retained(3) && audit.retained(4));
assert(size(D,3)==16);
Flags=ones(4,18); Flags(1,5)=0;
[~,audit]=build_paired_tensor(Results,[1 -1 0 0],Flags);
assert(audit.retained(3) && ~audit.retained(5));
% End-to-end custom-file path, flags, ID tracking, split matching and exports.
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
save(fullfile(folder,'input.mat'),'Results','Flags');
runCfg=struct('dataFile',fullfile(folder,'input.mat'), ...
    'contrastWeights',[1 -1 0 0],'contrastName','synthetic', ...
    'ranks',[1 2],'nStarts',2,'maxIter',80,'nSplitHalf',1, ...
    'outputDir',folder,'makePlots',false);
r=run_paired_tensor(runCfg);
assert(height(r.rankSummary)==2 && all(r.rankSummary.cpEnergyFit<=r.rankSummary.svdEnergyFit+1e-8));
assert(all(r.rankSummary.medianSplitHalfCongruence>=0 & r.rankSummary.medianSplitHalfCongruence<=1+1e-10));
assert(isfile(fullfile(folder,'rank_02_subject.csv')));
assert(numel(r.retainedSubjectIDs)==16);
fprintf('All paired tensor tests passed.\n');
end

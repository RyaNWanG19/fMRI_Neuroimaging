function test_svm_regional_reliance()
% Compare fast linear updates with actual predictions after paired swaps;
% verify OOB integration, missing features and unchanged baseline outputs.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','core'));
rng(91);
n = 80;
X = [randn(n,3)+[2 0 0]; randn(n,3)-[2 0 0]];
y = [ones(n,1); -ones(n,1)];
mdl = fitcsvm(X,y,'KernelFunction','linear','ClassNames',[-1 1]);
rngBefore = rng;
[loss,base] = svm_paired_region_reliance(mdl,X,30,117);
assert(isequal(rngBefore,rng), 'Importance changed the global RNG.');
stream = RandStream('mt19937ar','Seed',117);
reference = zeros(3,1);
for q=1:30
    swaps = rand(stream,n,3)<0.5;
    for j=1:3
        Z = X;
        ids = find(swaps(:,j));
        Z(ids,j) = X(n+ids,j);
        Z(n+ids,j) = X(ids,j);
        yp = predict(mdl,Z);
        reference(j) = reference(j) + mean(yp==y);
    end
end
assert(max(abs(loss-(base-reference/30)))<1e-12);
assert(loss(1)>0.25, 'Known predictive region should show substantial reliance.');
assert(isequal(loss,svm_paired_region_reliance(mdl,X,30,117)));

% Integration fixture: realistic stage dimensions, only three usable features.
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() clean_fixture(folder)); %#ok<NASGU>
oldNboot = getenv('NBOOT');
envCleanup = onCleanup(@() setenv('NBOOT',oldNboot)); %#ok<NASGU>
setenv('NBOOT','2');
Results = zeros(n,200,20,2);
for t=1:20
    Results(:,1:3,t,1) = X(1:n,:);
    Results(:,1:3,t,2) = X(n+1:end,:);
end
Results(:,200,:,:) = NaN;
AlreadyCanonical = true;
input = fullfile(folder,'fixture.mat');
save(input,'Results','AlreadyCanonical');
run_observed_stage_flexible(input,folder,[1 2],{'A','B'},'baseline');
cfg = struct('runImportance',true,'importanceRepeats',4,'importanceSeed',117);
run_observed_stage_flexible(input,folder,[1 2],{'A','B'},'reliance',cfg);
a = load(fullfile(folder,'observed_baseline.mat'));
b = load(fullfile(folder,'observed_reliance.mat'));
for name = {'TA_obs','W_obs','Haufe_obs','PairDiff_obs','TA_subj_oob'}
    assert(isequaln(a.(name{1}),b.(name{1})), 'Baseline output changed: %s',name{1});
end
assert(all(isnan(a.Reliance_obs),'all'));
assert(all(isfinite(b.Reliance_obs(1:3,:,:,:)),'all'));
assert(all(isnan(b.Reliance_obs(4:end,:,:,:)),'all'));
assert(all(b.Reliance_valid_boot(1:3,:,:)==2,'all'));
summary = readtable(fullfile(folder,'regional_reliance_reliance.csv'));
assert(height(summary)==200*20);
assert(all(summary.valid_boot(summary.roi_index==200)==0));
fprintf('SVM regional reliance reference and OOB integration checks passed.\n');
end

function clean_fixture(folder)
% Delete only the known fixture files, never recursively remove a path.
files = {'fixture.mat','observed_baseline.mat','observed_reliance.mat', ...
    'regional_reliance_reliance.csv'};
for j=1:numel(files)
    path = fullfile(folder,files{j});
    if exist(path,'file'), delete(path); end
end
rmdir(folder);
end

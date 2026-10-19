function test_svm_feature_selection()
% Training-only oracle, subject-copy grouping, sparse signal recovery,
% all-missing/constant features, determinism and Stage A persistence checks.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','core'));
rng(125);
n=80;
nuisance=randn(n,2);
P=[2+0.05*randn(n,1),nuisance,zeros(n,1),nan(n,1)];
N=[-2+0.05*randn(n,1),nuisance,zeros(n,1),nan(n,1)];
% Nonadjacent bootstrap copies must stay together.
ids=[(1:n)';(1:15)';(7:12)'];
cfg=struct('innerFolds',3,'featureCounts',[1 2 Inf],'seed',15);
rngBefore=rng;
[mdl,prep,info]=svm_nested_rfe(P(ids,:),N(ids,:),ids,cfg);
assert(isequal(rngBefore,rng),'Nested selection changed global RNG.');
assert(info.valid && info.chosenCount==1 && info.selected(1));
assert(~any(info.selected(4:5)) && ~any(info.eligible(4:5)));
assert(isequal(prep.selected,info.selected));
for id=unique(ids)'
    assert(numel(unique(info.rowFold(ids==id)))==1,'Bootstrap-copy leakage.');
end
assert(numel(mdl.Beta)==1);
[~,~,again]=svm_nested_rfe(P(ids,:),N(ids,:),ids,cfg);
assert(isequaln(info,again),'Selection is not deterministic.');

% Independent all-feature inner CV oracle: validation subjects once, means
% and scales computed only from the duplicated training rows of that fold.
for f=1:3
    tr=info.rowFold~=f;
    va=info.uniqueSubjectIDs(info.subjectFold==f);
    Z=[P(ids(tr),1:3);N(ids(tr),1:3)];
    mu=mean(Z,1); sd=std(Z,0,1);
    m=fitcsvm((Z-mu)./sd,[ones(nnz(tr),1);-ones(nnz(tr),1)], ...
        'KernelFunction','linear','BoxConstraint',1,'ClassNames',[-1 1]);
    y=[ones(numel(va),1);-ones(numel(va),1)];
    acc=mean(predict(m,([P(va,1:3);N(va,1:3)]-mu)./sd)==y);
    assert(abs(info.foldAccuracy(end,f)-acc)<1e-12);
end
% Requested counts beyond availability clip explicitly rather than indexing
% dropped features. All-degenerate data return a marked invalid result.
[~,~,clipped]=svm_nested_rfe(P,N,(1:n)', ...
    struct('innerFolds',2,'featureCounts',[50 Inf]));
assert(clipped.valid && clipped.chosenCount==3);
[empty,~,invalid]=svm_nested_rfe(zeros(n,2),zeros(n,2),(1:n)',cfg);
assert(isempty(empty) && ~invalid.valid && all(isnan(invalid.foldAccuracy),'all'));

folder=tempname; mkdir(folder);
cleanup=onCleanup(@() clean_fixture(folder)); %#ok<NASGU>
oldNboot=getenv('NBOOT');
envCleanup=onCleanup(@() setenv('NBOOT',oldNboot)); %#ok<NASGU>
setenv('NBOOT','2');
Results=zeros(n,200,20,2);
for t=1:20
    Results(:,1:5,t,1)=P;
    Results(:,1:5,t,2)=N;
end
AlreadyCanonical=true;
input=fullfile(folder,'fixture.mat');
save(input,'Results','AlreadyCanonical');
run_observed_stage_flexible(input,folder,[1 2],{'A','B'},'baseline');
stageCfg=struct('runFeatureSelection',true,'featureSelection',cfg, ...
    'runImportance',true,'importanceRepeats',4);
run_observed_stage_flexible(input,folder,[1 2],{'A','B'},'selected',stageCfg);
a=load(fullfile(folder,'observed_baseline.mat'));
b=load(fullfile(folder,'observed_selected.mat'));
for name={'TA_obs','W_obs','Haufe_obs','PairDiff_obs','TA_subj_oob','split_list'}
    assert(isequaln(a.(name{1}),b.(name{1})),'Baseline output changed: %s',name{1});
end
fs=b.FeatureSelection;
assert(fs.enabled && all(fs.valid_boot==2,'all'));
assert(all(fs.selection_frequency(1,:,:)==1,'all'));
assert(all(fs.selection_frequency(4:end,:,:)==0,'all'));
assert(all(fs.accuracy_obs==1,'all'));
for boot=1:2
    q=fs.bootstrap{boot};
    assert(isempty(intersect(q.uniqueSubjectIDs,b.split_list{boot}.test)), ...
        'OOB subjects entered inner selection.');
    assert(all(q.count==1,'all'));
    assert(all(isfinite(q.reliance(1,:,:)),'all'));
    assert(all(isnan(q.reliance(2:end,:,:)),'all'));
end
regional=readtable(fullfile(folder,'regional_selection_selected.csv'));
comparison=readtable(fullfile(folder,'selection_accuracy_selected.csv'));
assert(height(regional)==200*20 && height(comparison)==20);
assert(all(comparison.selected_accuracy==1));
assert(all(comparison.paired_accuracy_difference==0));
fprintf('Nested SVM-RFE oracle, leakage and Stage A integration checks passed.\n');
end

function clean_fixture(folder)
files={'fixture.mat','observed_baseline.mat','observed_selected.mat', ...
    'regional_selection_selected.csv','selection_accuracy_selected.csv', ...
    'regional_reliance_selected.csv'};
for j=1:numel(files)
    path=fullfile(folder,files{j});
    if exist(path,'file'), delete(path); end
end
rmdir(folder);
end

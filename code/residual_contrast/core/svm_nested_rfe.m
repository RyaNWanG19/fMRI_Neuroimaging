function [mdl, prep, info] = svm_nested_rfe(Xpositive, Xnegative, subjectIDs, cfg)
% Training-only nested linear SVM-RFE for paired conditions.
% Inputs contain OUTER TRAINING subjects only, including bootstrap copies.
% Inner folds group both conditions and all copies of a subject. Inner
% validation scores each distinct subject once; training retains multiplicity.
% Imputation, standardization and the entire RFE path are refitted inside
% each inner fold. Maximum mean fold accuracy chooses count, ties favor fewer
% regions. Refit that RFE path on the complete outer training sample.
% No held-out data enter this function. cfg is normalized by svm_rfe_config.

if nargin<4, cfg=struct(); end
cfg=svm_rfe_config(cfg);
validateattributes(Xpositive,{'numeric'},{'2d','nonempty'});
assert(isequal(size(Xpositive),size(Xnegative)), 'Condition matrices must match.');
assert(~any(isinf(Xpositive(:))) && ~any(isinf(Xnegative(:))), ...
    'Infinite input values are not supported; use NaN for missing data.');
ids=subjectIDs(:);
validateattributes(ids,{'numeric'},{'vector','finite','integer','positive'});
assert(numel(ids)==size(Xpositive,1), 'One subject ID per paired training row required.');
[uniqueIDs,firstRow,group]=unique(ids);
assert(numel(uniqueIDs)>=2, 'Need at least two distinct outer training subjects.');
assert(isequaln(Xpositive,Xpositive(firstRow(group),:)) && ...
    isequaln(Xnegative,Xnegative(firstRow(group),:)), ...
    'Bootstrap copies of a subject must have identical condition features.');
nFolds=min(cfg.innerFolds,numel(uniqueIDs));
stream=RandStream('mt19937ar','Seed',cfg.seed);
order=randperm(stream,numel(uniqueIDs));
subjectFold=zeros(numel(uniqueIDs),1);
subjectFold(order)=mod(0:numel(uniqueIDs)-1,nFolds)+1;
rowFold=subjectFold(group);
counts=cfg.featureCounts;
foldAccuracy=nan(numel(counts),nFolds);
foldEffectiveCount=zeros(numel(counts),nFolds);
for f=1:nFolds
    tr=rowFold~=f;
    % Validation excludes duplicate copies; each condition appears once.
    va=firstRow(subjectFold==f);
    Xtr=double([Xpositive(tr,:);Xnegative(tr,:)]);
    ytr=[ones(nnz(tr),1);-ones(nnz(tr),1)];
    [Z,foldPrep]=fit_preprocessing(Xtr);
    if isempty(Z), continue; end
    V=apply_preprocessing(double([Xpositive(va,:);Xnegative(va,:)]),foldPrep);
    yva=[ones(numel(va),1);-ones(numel(va),1)];
    [models,masks]=rfe_path(Z,ytr,counts,cfg.boxConstraint);
    for j=1:numel(counts)
        foldAccuracy(j,f)=mean(predict(models{j},V(:,masks{j}))==yva);
        foldEffectiveCount(j,f)=nnz(masks{j});
    end
end
% A candidate needs evaluable predictions in every planned inner fold.
meanAccuracy=mean(foldAccuracy,2);
eligibleCandidate=all(isfinite(foldAccuracy),2);
mdl=[]; prep=struct();
info=struct('valid',false,'status','no_evaluable_inner_candidate', ...
    'candidateCounts',counts,'foldAccuracy',foldAccuracy, ...
    'foldEffectiveCount',foldEffectiveCount,'meanAccuracy',meanAccuracy, ...
    'uniqueSubjectIDs',uniqueIDs,'subjectFold',subjectFold,'rowFold',rowFold, ...
    'chosenRequestedCount',NaN,'chosenCount',NaN,'innerAccuracy',NaN, ...
    'selected',false(1,size(Xpositive,2)),'eligible',false(1,size(Xpositive,2)));
if ~any(eligibleCandidate), return; end
best=max(meanAccuracy(eligibleCandidate));
chosen=find(eligibleCandidate & meanAccuracy>=best-1e-12,1,'first');
[Z,prep]=fit_preprocessing(double([Xpositive;Xnegative]));
if isempty(Z), info.status='no_eligible_outer_features'; return; end
y=[ones(size(Xpositive,1),1);-ones(size(Xpositive,1),1)];
% Use the SAME predeclared elimination path as inner training, rather than
% taking the top k weights from a single full fit.
[models,masks]=rfe_path(Z,y,counts,cfg.boxConstraint);
mdl=models{chosen};
info.eligible=prep.keep;
original=find(prep.keep);
info.selected(original(masks{chosen}))=true;
prep.selected=info.selected;
info.chosenRequestedCount=counts(chosen);
info.chosenCount=nnz(info.selected);
info.innerAccuracy=meanAccuracy(chosen);
info.valid=true;
info.status='ok';
end

function [models,masks]=rfe_path(Z,y,counts,boxConstraint)
% Batched elimination at grid counts: rank squared standardized coefficients
% and remove the weakest batch, refitting at every remaining count.
effective=min(counts,size(Z,2));
levels=sort(unique(effective),'descend');
% Always start with all eligible features, even if absent from count grid.
levels=unique([size(Z,2),levels]);
levels=sort(levels,'descend');
active=1:size(Z,2);
models=cell(numel(counts),1); masks=cell(numel(counts),1);
for q=1:numel(levels)
    if q>1
        % Deterministic tie break by original standardized feature position.
        [~,rank]=sortrows([-double(mdl.Beta(:)).^2,active(:)],[1 2]);
        active=sort(active(rank(1:levels(q))));
    end
    mdl=fitcsvm(Z(:,active),y,'KernelFunction','linear', ...
        'BoxConstraint',boxConstraint,'ClassNames',[-1 1]);
    mask=false(1,size(Z,2)); mask(active)=true;
    matches=find(effective==levels(q));
    for j=matches(:)'
        models{j}=mdl; masks{j}=mask;
    end
end
end

function [Z,prep]=fit_preprocessing(X)
prep.mu=mean(X,1,'omitnan');
for j=1:size(X,2)
    X(isnan(X(:,j)),j)=prep.mu(j);
end
prep.sd=std(X,0,1);
prep.keep=isfinite(prep.mu) & isfinite(prep.sd) & prep.sd>0;
Z=(X(:,prep.keep)-prep.mu(prep.keep))./prep.sd(prep.keep);
end

function Z=apply_preprocessing(X,prep)
for j=1:size(X,2)
    X(isnan(X(:,j)),j)=prep.mu(j);
end
Z=(X(:,prep.keep)-prep.mu(prep.keep))./prep.sd(prep.keep);
end

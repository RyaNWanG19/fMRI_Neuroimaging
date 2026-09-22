root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
base = fullfile(root,'data','task_residual','jhpce_outputs_wm_population_omnibus_full','Body_LoadDiff_vs_Face_LoadDiff');
build = fileparts(mfilename('fullpath'));
models = {'cHRF','cHRFderiv','sHRF'};
atlas = load_atlas('canlab2018');
labels = string(atlas.labels(:));
assert(numel(labels)==489);
sig=false(489,3); boot=sig; ps=zeros(489,3); stats=ps;
for k=1:3
    f=dir(fullfile(base,models{k},'**','temporal_omnibus_results.mat'));
    assert(numel(f)==1);
    a=load(fullfile(f.folder,f.name),'results'); r=a.results;
    assert(r.config.nPerm==99999 && r.sensitivity.nBootstrap==99999);
    assert(isequal(string(r.roiIDs),string((1:489)')));
    assert(isequal(r.timeSec(:),(0:40)'*.72));
    R{k}=r;
    sig(:,k)=r.significant; boot(:,k)=r.sensitivity.table.bootstrap_significant;
    ps(:,k)=r.p_holm; stats(:,k)=r.S_obs;
    payload.models(k).name=models{k};
    payload.models(k).source=fullfile(f.folder,f.name);
    payload.models(k).significant=sum(sig(:,k));
    payload.models(k).bootstrap=sum(boot(:,k));
    payload.models(k).both=sum(sig(:,k)&boot(:,k));
    payload.models(k).disagreement=sum(sig(:,k)~=boot(:,k));
    payload.models(k).nearMC=sum(r.monteCarloNearThreshold);
    payload.models(k).effectiveN=[min(r.effectiveN) max(r.effectiveN)];
end
payload.timeSec=R{1}.timeSec;
payload.labels=labels;
payload.pHolm=ps;
payload.S=stats;
payload.significant=sig;
payload.sharedAll=sum(all(sig,2));
payload.sharedBothAll=sum(all(sig & boot,2));
payload.any=sum(any(sig,2));
payload.none=sum(~any(sig,2));
patterns=[1 1 1;1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1];
payload.patterns=patterns;
payload.patternCounts=sum(all(reshape(sig,489,1,3)==reshape(patterns,1,7,3),3),1);
% Examples: V1 left as familiar shared ROI, strongest remaining common ROI
% by minimum S across models, strongest cHRF-only ROI by cHRF S.
common=find(all(sig,2));
assert(sig(1,1)&&sig(1,2)&&sig(1,3));
[~,ord]=sort(min(stats(common,:),[],2),'descend');
best=common(ord); best(best==1)=[];
only=find(sig(:,1)&~sig(:,2)&~sig(:,3));
[~,j]=max(stats(only,1));
chosen=[1;best(1);only(j)];
payload.selectionNote='Illustrative post-selection curves: ROI 1 (left V1), remaining shared ROI with largest minimum S across models, cHRF-only ROI with largest cHRF S. No lag-wise testing or model-difference test.';
for j=1:3
    roi=chosen(j);
    payload.examples(j).id=roi;
    payload.examples(j).name=labels(roi);
    payload.examples(j).pHolm=ps(roi,:);
    payload.examples(j).S=stats(roi,:);
    for k=1:3
        payload.examples(j).mean(k,:)=R{k}.mean_D(roi,:);
        payload.examples(j).sd(k,:)=R{k}.sd_D(roi,:);
        payload.examples(j).n(k)=R{k}.effectiveN(roi);
    end
end
payload.notes='N=410 input, 402-410 complete subjects per ROI. One Holm family of 489 ROIs per model. Independence confirmed by data owner; null symmetry assumed. CANlab2018 labels from local atlas matched by ROI index. Stored time axis follows project convention 0:0.72:28.8 sec.';
fid=fopen(fullfile(build,'results.json'),'w'); assert(fid>0); fprintf(fid,'%s',jsonencode(payload)); fclose(fid);
disp(payload.patternCounts); disp(chosen); disp(labels(chosen));
fprintf('common both methods all3=%d\n',payload.sharedBothAll);

% Read-only audit of source data and saved results. Does not rerun inference outputs.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
base=fullfile(root,'data','task_residual');
names={'cHRF','cHRFderiv','sHRF'};
for k=1:3
    f=dir(fullfile(base,'jhpce_outputs_wm_population_omnibus_full', ...
        'Body_LoadDiff_vs_Face_LoadDiff',names{k},'**','temporal_omnibus_results.mat'));
    a=load(fullfile(f(1).folder,f(1).name),'results'); r=a.results;
    data=load(fullfile(base,['WM' names{k} '.mat']),'Results');
    % Independent Holm calculation across the entire planned family.
    [p,ix]=sort(r.p_raw(:)); m=numel(p);
    adj=min(1,cummax(p.*(m:-1:1)')); ph=zeros(m,1); ph(ix)=adj;
    assert(max(abs(ph-r.p_holm))<1e-12);
    assert(isequal(r.significant,ph<=r.config.alpha));
    assert(isequal(string(r.roiIDs(:)),string((1:m)')));
    for roi=[78 211]
        % Explicit condition arithmetic, time x subject -> subject x time.
        v=double(data.Results(roi,:,:,:));
        x=squeeze((v(:,:,:,5)-v(:,:,:,1))-(v(:,:,:,6)-v(:,:,:,2)))';
        keep=all(isfinite(x),2); n=sum(keep); x=x(keep,:);
        assert(isequal(keep,r.subjectMasks(:,roi)));
        mu=mean(x,1); sd=std(x,0,1); tt=mu./(sd/sqrt(n)); ss=sum(tt.^2);
        assert(max(abs(mu-r.mean_D(roi,:)))<1e-12);
        assert(max(abs(sd-r.sd_D(roi,:)))<1e-12);
        assert(max(abs(tt-r.T_obs(roi,:)))<1e-10);
        assert(abs(ss-r.S_obs(roi))<1e-8);
        % Regenerate the full saved sign stream, independent ROI-only arithmetic.
        stream=RandStream('mt19937ar','Seed',r.config.seed);
        exceed=0; maxDirectError=0;
        for first=1:500:r.config.nPerm
            nb=min(500,r.config.nPerm-first+1);
            signs=2*(rand(stream,numel(keep),nb)>=0.5)-1;
            means=(signs(keep,:)'*x)/n;
            vars=(sum(x.^2,1)-n*means.^2)/(n-1);
            sb=sum(means.^2./(vars/n),2);
            if first==1
                for b=1:200
                    xb=x.*signs(keep,b);
                    direct=sum((mean(xb,1)./(std(xb,0,1)/sqrt(n))).^2);
                    maxDirectError=max(maxDirectError,abs(direct-sb(b)));
                end
                assert(maxDirectError<1e-8);
            end
            exceed=exceed+sum(sb>=ss);
        end
        praw=(exceed+1)/(r.config.nPerm+1);
        assert(praw==r.p_raw(roi));
        fprintf('PASS %s ROI %d: source mean/SD/T/S, mask, %d flips; direct error %.3g; p=%.8f adjusted=%.8f MCflag=%d\n', ...
            names{k},roi,r.config.nPerm,maxDirectError,praw,ph(roi),r.monteCarloNearThreshold(roi));
    end
    clear data
end
disp('ALL SIX SOURCE-DATA AUDITS PASSED; all-model Holm and labels also matched.');

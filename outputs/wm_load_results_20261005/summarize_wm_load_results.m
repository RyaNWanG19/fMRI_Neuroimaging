function summary = summarize_wm_load_results(outputDir, makeAnatomy)
% Summarize saved full-run results without modifying or rerunning inference.
% Figures use saved primary decisions, corrected over 1956 tests per model.
% Anatomical labels follow the project's canlab2018 index convention.
if nargin<1 || isempty(outputDir)
    outputDir=fullfile(pwd,'outputs',['wm_load_summary_' datestr(now,'yyyymmddTHHMMSS')]);
end
if nargin<2, makeAnatomy=true; end
assert(~isfolder(outputDir),'Use a new output directory to preserve earlier outputs.');
mkdir(outputDir); figDir=fullfile(outputDir,'figures'); mkdir(figDir);
root=fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
base=fullfile(root,'data','task_residual','wm_load_omnibus','AllFourCategory_LoadDiff');
models={'cHRF','cHRFderiv','sHRF'}; cats={'body','faces','places','tools'};
colors=[22 117 169;212 119 38;106 74 162]/255;
R=cell(1,3); sig=false(489,4,3); rmsMean=zeros(489,4,3);
summary=struct('models',{models},'categories',{cats},'sourceFiles',{{}}, ...
    'correctionFamily','1956 ROI-category hypotheses per model; no across-model correction');
summary.counts=zeros(4,3); summary.mcFlags=zeros(4,3); summary.floorCounts=zeros(4,3);
summary.mcFlaggedSignificant=zeros(4,3); summary.auditMaxErrors=zeros(3,4);
for k=1:3
    f=dir(fullfile(base,models{k},'**','temporal_omnibus_results.mat'));
    assert(numel(f)==1,'Expected exactly one full run for %s.',models{k});
    summary.sourceFiles{k}=fullfile(f.folder,f.name);
    a=load(summary.sourceFiles{k},'results'); r=a.results;
    assert(r.config.nPerm==99999 && strcmp(r.config.runMode,'inference') && r.inferenceValid);
    assert(strcmp(r.config.hrfModelName,models{k}) && r.config.alpha==0.05);
    assert(isequal(r.timeSec(:),(0:40)'*.72));
    assert(isequal(r.testMapping.atlas_roi_index,repmat((1:489)',4,1)));
    assert(isequal(string(r.testMapping.category),repelem(string(cats(:)),489)));
    assert(isequal(r.testMapping.condition_0back,repelem((1:4)',489)));
    assert(isequal(r.testMapping.condition_2back,repelem((5:8)',489)));
    assert(isequal(r.p_holm,holm_adjust(r.p_raw)));
    assert(isequal(r.significant,r.p_holm<=.05 & r.effectiveN>=2));
    assert(~isfield(r,'sensitivity'),'Update report for bootstrap results before proceeding.');
    % Check all observed summaries against original subject-level residuals.
    raw=load(fullfile(root,'data','task_residual',['WM' models{k} '.mat']),'Results');
    for c=1:4
        for roi=1:489
            ix=(c-1)*489+roi;
            x=squeeze(double(raw.Results(roi,:,:,c+4))-double(raw.Results(roi,:,:,c)))';
            keep=all(isfinite(x),2); x=x(keep,:);
            assert(isequal(keep,r.subjectMasks(:,ix)) && sum(keep)==r.effectiveN(ix));
            mu=mean(x,1); sd=std(x,0,1); tt=mu./(sd/sqrt(size(x,1))); ss=sum(tt.^2);
            assert(all(isfinite(tt)),'Unexpected degenerate case needs review.');
            e=[max(abs(mu-r.mean_D(ix,:))),max(abs(sd-r.sd_D(ix,:))), ...
                max(abs(tt-r.T_obs(ix,:))),abs(ss-r.S_obs(ix))];
            summary.auditMaxErrors(k,:)=max(summary.auditMaxErrors(k,:),e);
        end
    end
    assert(all(summary.auditMaxErrors(k,:)<1e-8)); clear raw
    R{k}=r; sig(:,:,k)=reshape(r.significant,489,4);
    summary.counts(:,k)=sum(sig(:,:,k),1)';
    summary.mcFlags(:,k)=sum(reshape(r.monteCarloNearThreshold,489,4),1)';
    summary.mcFlaggedSignificant(:,k)=sum(reshape(r.monteCarloNearThreshold&r.significant,489,4),1)';
    summary.floorCounts(:,k)=sum(reshape(r.p_raw==1/100000,489,4),1)';
    summary.effectiveN(k,:)=[min(r.effectiveN) max(r.effectiveN)];
    rmsMean(:,:,k)=reshape(sqrt(mean(r.mean_D.^2,2)),489,4);
    fprintf('SOURCE AUDIT PASS %s; max mean/SD/T/S errors: %s\n',models{k},mat2str(summary.auditMaxErrors(k,:),4));
end
patterns=[1 1 1;1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
patternNames={'All three','cHRF + deriv','cHRF + sHRF','deriv + sHRF', ...
    'cHRF only','deriv only','sHRF only','None'};
summary.patternNames=patternNames; summary.patternCounts=zeros(4,8);
for c=1:4
    for p=1:8
        summary.patternCounts(c,p)=sum(all(reshape(sig(:,c,:),489,3)==patterns(p,:),2));
    end
end
assert(all(sum(summary.patternCounts,2)==489));
summary.sharedAll=summary.patternCounts(:,1);
summary.anyModel=489-summary.patternCounts(:,8);
summary.significant=sig;
summary.sharedEveryCategoryAndModel=sum(all(all(sig,3),2));
summary.sharedCategoryCountByROI=sum(all(sig,3),2);
summary.medianRMS=squeeze(median(rmsMean,1));
summary.rmsMean=rmsMean;
summary.timeSec=R{1}.timeSec;
atlasObj=load_atlas('canlab2018'); labels=string(atlasObj.labels(:));
assert(numel(labels)==489); summary.atlasLabels=labels;

% Counts with numeric labels; no correction terminology in the title.
f=newfig([1150 640]); ax=axes(f); b=bar(ax,summary.counts,'grouped');
for k=1:3
    b(k).FaceColor=colors(k,:);
    text(ax,b(k).XEndPoints,b(k).YEndPoints+8,string(summary.counts(:,k)), ...
        'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
end
set(ax,'XTick',1:4,'XTickLabel',cats,'FontSize',16,'Box','off');
ylim(ax,[0 489]); ylabel(ax,'Number of regions');
title(ax,'Number of significant regions','FontSize',23);
subtitle(ax,'Memory load difference by stimulus category');
legend(ax,models,'Location','northoutside','Orientation','horizontal');
savefigpng(f,figDir,'significant_region_counts');

f=newfig([1400 850]); tl=tiledlayout(f,2,2,'TileSpacing','compact');
for c=1:4
    ax=nexttile(tl); vals=summary.patternCounts(c,:);
    b=bar(ax,vals,'FaceColor','flat'); b.CData=repmat([.57 .66 .71],8,1); b.CData(1,:)=colors(1,:);
    text(ax,1:8,vals+8,string(vals),'HorizontalAlignment','center','FontSize',12);
    set(ax,'XTick',1:8,'XTickLabel',patternNames,'FontSize',12,'Box','off'); xtickangle(ax,28);
    ylim(ax,[0 489]); ylabel(ax,'Number of regions'); title(ax,cats{c},'FontSize',19);
end
title(tl,'Overlap of significant regions across HRF models','FontSize',23);
subtitle(tl,'Mutually exclusive groups, including regions significant in none');
savefigpng(f,figDir,'model_overlap');

f=newfig([1250 720]); tl=tiledlayout(f,2,2,'TileSpacing','compact');
for c=1:4
    ax=nexttile(tl); hold(ax,'on');
    for k=1:3
        plot(ax,sort(rmsMean(:,c,k)),(1:489)/489,'Color',colors(k,:),'LineWidth',2);
    end
    title(ax,cats{c}); xlabel(ax,'RMS of population-average residual contrast');
    ylabel(ax,'Fraction of regions'); set(ax,'FontSize',12,'Box','off');
    if c==1, lg=legend(ax,models,'Orientation','horizontal'); lg.Layout.Tile='south'; end
end
title(tl,'Residual contrast magnitude across all 489 regions','FontSize',22);
savefigpng(f,figDir,'residual_magnitude_distributions');

% First atlas index in each exclusive pattern: never substitute another group.
examplePatterns=[1 5 6 7 8]; summary.examples=struct([]);
for c=1:4
    selected=nan(1,5); flags=reshape(sig(:,c,:),489,3);
    for g=1:5
        candidates=find(all(flags==patterns(examplePatterns(g),:),2));
        if ~isempty(candidates), selected(g)=candidates(1); end
    end
    summary.examples(c).category=cats{c}; summary.examples(c).roi=selected;
    summary.examples(c).pHolm=nan(5,3);
    limMean=0; limT=0;
    for roi=selected(isfinite(selected))
        ix=(c-1)*489+roi;
        for k=1:3
            limMean=max(limMean,max(abs(R{k}.mean_D(ix,:))));
            limT=max(limT,max(abs(R{k}.T_obs(ix,:))));
        end
    end
    for mode=1:2
        f=newfig([1450 850]); tl=tiledlayout(f,2,3,'TileSpacing','compact');
        for g=1:5
            ax=nexttile(tl); hold(ax,'on'); roi=selected(g);
            titleText=patternNames{examplePatterns(g)};
            if isnan(roi)
                axis(ax,'off'); text(ax,.1,.5,'No regions in this group','FontSize',14);
                title(ax,titleText); continue
            end
            ix=(c-1)*489+roi; h=gobjects(1,3);
            for k=1:3
                if mode==1, y=R{k}.mean_D(ix,:); else, y=R{k}.T_obs(ix,:); end
                h(k)=plot(ax,R{k}.timeSec,y,'Color',colors(k,:),'LineWidth',1.8);
                summary.examples(c).pHolm(g,k)=R{k}.p_holm(ix);
            end
            title(ax,{titleText,sprintf('%s (ROI %d)',labels(roi),roi)},'Interpreter','none','FontSize',13);
            yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[0 28.8]); xticks(ax,0:7.2:28.8);
            xlabel(ax,'Time (s)'); set(ax,'FontSize',11,'Box','off');
            if mode==1, ylabel(ax,'Mean residual contrast'); ylim(ax,[-1 1]*max(limMean*1.1,eps));
            else, ylabel(ax,'T statistic'); ylim(ax,[-1 1]*max(limT*1.1,eps)); end
            if g==1, lg=legend(ax,h,models,'Orientation','horizontal'); lg.Layout.Tile='south'; end
        end
        if mode==1, what='Mean residual contrast trajectories'; suffix='mean';
        else, what='Standardized residual contrast trajectories'; suffix='T'; end
        title(tl,sprintf('%s: %s',what,cats{c}),'FontSize',22);
        subtitle(tl,'First ROI in each group, selected after testing. No individual-lag significance implied.');
        savefigpng(f,figDir,[cats{c} '_examples_' suffix]);
    end
end
save(fullfile(outputDir,'summary_data.mat'),'summary','-v7.3');
writeReport(summary,outputDir);

if makeAnatomy
    for c=1:4
        for k=1:4
            if k<=3, selected=sig(:,c,k); tag=models{k}; color=colors(k,:);
            else, selected=all(sig(:,c,:),3); tag='shared'; color=[35 140 110]/255; end
            drawAnatomy(atlasObj,selected,cats{c},tag,color, ...
                fullfile(figDir,[cats{c} '_' tag '_orthviews.png']));
        end
    end
end
fprintf('REPORT SAVED: %s\n',outputDir);
disp(summary.counts); disp(summary.patternCounts);
end

function f=newfig(sz)
f=figure('Visible','off','Color','w','Position',[50 50 sz]);
end
function savefigpng(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',220); close(f);
end
function drawAnatomy(atlasObj,mask,category,tag,color,file)
if ~any(mask), return; end
sub=select_atlas_subset(atlasObj,find(mask),'deterministic'); r=atlas2region(sub);
orthviews(r,repmat({color},1,numel(r)),'solid'); xyz=[0 -20 20];
spm_orthviews('Reposition',xyz'); drawnow;
fig=spm_figure('FindWin','Graphics');
txt=findall(fig,'Type','text');
for j=1:numel(txt)
    value=get(txt(j),'String');
    if ischar(value)
        tok=regexp(value,'^\s*([xyz])\s*=','tokens','once');
        if ~isempty(tok), set(txt(j),'String',sprintf('%s = %g',tok{1},xyz('xyz'==tok{1}))); end
    end
end
annotation(fig,'textbox',[.03 .91 .94 .07],'String', ...
    sprintf('%s: %s, %d significant regions',category,tag,sum(mask)), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',14,'FontWeight','bold', ...
    'Interpreter','none','Color','k','BackgroundColor','w');
annotation(fig,'textbox',[.03 .01 .94 .065],'String', ...
    sprintf('2-back minus 0-back. MNI [0 -20 20] mm.\nColor shows ROI membership, not magnitude or direction. Only selected slices are visible.'), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9,'Color','k','BackgroundColor','w');
drawnow; exportgraphics(fig,file,'Resolution',200);
fprintf('ANATOMY SAVED %s %s\n',category,tag);
end

function writeReport(s,out)
fid=fopen(fullfile(out,'results_report.md'),'w'); assert(fid>0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# WM load-difference residual results\n\n');
fprintf(fid,['Each test assesses whether the population-average 2-back minus 0-back residual curve departs from zero somewhere across 41 sampled lags (0-28.8 s). ' ...
    'The four categories are body, faces, places, and tools. Full runs used 99,999 subject sign flips, alpha 0.05, and Holm correction across 1,956 ROI-category tests separately within each HRF model.\n\n']);
fprintf(fid,'## Significant regions\n\n| Category | cHRF | cHRFderiv | sHRF | Shared by all three | At least one model | None |\n|---|---:|---:|---:|---:|---:|---:|\n');
for c=1:4
    fprintf(fid,'| %s | %d (%.1f%%) | %d (%.1f%%) | %d (%.1f%%) | %d | %d | %d |\n', ...
        s.categories{c},s.counts(c,1),s.counts(c,1)/489*100,s.counts(c,2),s.counts(c,2)/489*100, ...
        s.counts(c,3),s.counts(c,3)/489*100,s.sharedAll(c),s.anyModel(c),s.patternCounts(c,8));
end
fprintf(fid,['\nThe saved decisions show systematic load-dependent residual contrasts in many regions in every stimulus category. ' ...
    'sHRF has fewer significant regions than either canonical model in all four categories. This is a descriptive comparison: counts alone do not establish that sHRF fits better or that the models differ statistically.\n\n' ...
    '%d regions are significant in every category under every model. Per-category overlaps describe agreement between decisions, not independent replication.\n\n'],s.sharedEveryCategoryAndModel);
fprintf(fid,'## Exclusive model overlap\n\n| Category | All three | cHRF + deriv only | cHRF + sHRF only | deriv + sHRF only | cHRF only | deriv only | sHRF only | None |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for c=1:4, fprintf(fid,'| %s | %d | %d | %d | %d | %d | %d | %d | %d |\n',s.categories{c},s.patternCounts(c,:)); end
fprintf(fid,'\n## Monte Carlo resolution and robustness\n\n');
fprintf(fid,['The minimum raw p-value is 0.00001. At the minimum, 1,956-test Holm adjustment yields 0.01956. ' ...
    'With zero exceedances, the code estimates a raw-p Monte Carlo SE near 0.00001. Its upper two-SE input is then near 0.00003, which can produce a Holm upper sensitivity value near 0.05868. ' ...
    'Because Holm uses cumulative maxima, this early sorted-p bound can propagate to many later hypotheses. The broad flags therefore do not mean each reported p-value is individually close to 0.05. ' ...
    'This is a heuristic two-SE sensitivity screen, not a simultaneous confidence interval. It flags possible finite-permutation sensitivity, rather than demonstrating that all flagged decisions would change.\n\n' ...
    '| Category | Flagged significant cHRF | Flagged significant deriv | Flagged significant sHRF |\n|---|---:|---:|---:|\n']);
for c=1:4, fprintf(fid,'| %s | %d | %d | %d |\n',s.categories{c},s.mcFlaggedSignificant(c,:)); end
fprintf(fid,['\nThese runs do not contain centered-bootstrap sensitivity results. The earlier interaction bootstrap findings cannot be transferred to these tests. ' ...
    'Before treating borderline classifications as stable, increase the permutation count or use a statistically justified Monte Carlo precision assessment. No additional inference was run for this report.\n\n']);
fprintf(fid,'## Residual magnitude\n\nMedian RMS of each ROI population-average curve, over all 489 ROIs (native residual units):\n\n| Category | cHRF | cHRFderiv | sHRF |\n|---|---:|---:|---:|\n');
for c=1:4, fprintf(fid,'| %s | %.6g | %.6g | %.6g |\n',s.categories{c},s.medianRMS(c,:)); end
fprintf(fid,['\nRMS describes curve magnitude and uses no uncertainty weighting. Significance instead depends on the time-specific mean-to-SE ratios and the permutation distribution. ' ...
    'Magnitude distributions are descriptive and have no between-model p-value.\n\n## Example trajectories\n\n' ...
    'Each category has mean and T-statistic figures for five patterns: all three models, each model only, and none. ' ...
    'Examples are the first atlas index in each pattern, chosen after testing. Empty patterns remain explicitly empty. ' ...
    'The example curves are not representative estimates for their groups. Each category figure uses common y limits across its panels.\n\n']);
for c=1:4
    fprintf(fid,'### %s\n\n| Pattern | ROI | Label | cHRF adjusted p | deriv adjusted p | sHRF adjusted p |\n|---|---:|---|---:|---:|---:|\n',s.categories{c});
    pi=[1 5 6 7 8];
    for g=1:5
        roi=s.examples(c).roi(g);
        if isnan(roi), fprintf(fid,'| %s | - | No regions | - | - | - |\n',s.patternNames{pi(g)});
        else, fprintf(fid,'| %s | %d | %s | %.6g | %.6g | %.6g |\n',s.patternNames{pi(g)},roi,s.atlasLabels(roi),s.examples(c).pHolm(g,:)); end
    end
    fprintf(fid,'\n');
end
fprintf(fid,['## Anatomy and interpretation boundaries\n\n' ...
    '- Orthviews use deterministic canlab2018 parcel boundaries and the project atlas-index convention. Historical atlas-version identity cannot be independently established from ordinal IDs alone.\n' ...
    '- All anatomical maps use the same MNI coordinates [0 -20 20] mm. Three slices do not display every significant region. Color denotes ROI membership, not effect sign, magnitude, or voxelwise significance.\n' ...
    '- A significant omnibus test provides evidence of a nonzero population-average load contrast somewhere in the sampled window. It does not locate significant lags.\n' ...
    '- Shared mismatch under 2-back and 0-back can still cancel. These are not absolute condition-specific residual tests.\n' ...
    '- The sign-flip inference is conditional on the documented independence and joint null-symmetry assumptions. Numerical agreement does not verify those assumptions.\n' ...
    '- Correction covers four categories within each model, not all models or the earlier interaction analysis together. Different correction families prevent direct interpretation of count changes as effect changes.\n\n' ...
    '## Checks and provenance\n\nAll 5,868 observed ROI-category/model means, SDs, T curves, S statistics and subject masks were recomputed from the local source residual arrays and matched within numerical tolerance. ' ...
    'Saved adjusted p-values and labels match the 1,956-test Holm calculation. This check does not independently regenerate the 99,999 permutations or audit upstream HRF fitting.\n\n']);
for k=1:3, fprintf(fid,'- %s: N %d-%d, source `%s`\n',s.models{k},s.effectiveN(k,:),s.sourceFiles{k}); end
fprintf(fid,'\n## Slide-ready interpretation\n\n');
fprintf(fid,['Across all four stimulus categories, the population-average 2-back minus 0-back residual contrast remained significant in many regions. ' ...
    'The sHRF model yielded fewer significant regions than cHRF or cHRFderiv, while a subset of regions passed under all three models. ' ...
    'These results identify load-dependent residual structure. They do not establish superior fit from counts alone or exclude mismatch shared across loads. ' ...
    'Finite-permutation uncertainty warrants caution when interpreting exact significance boundaries.\n']);
% Complete shared anatomical sets, not restricted to displayed slices.
fid2=fopen(fullfile(out,'shared_region_lists.md'),'w'); assert(fid2>0); close2=onCleanup(@() fclose(fid2)); %#ok<NASGU>
fprintf(fid2,'# Regions significant in all three models\n\nAtlas labels follow local canlab2018 index order. Lists use saved primary decisions.\n\n');
for c=1:4
    fprintf(fid2,'## %s\n\n| Atlas ROI | Label |\n|---:|---|\n',s.categories{c});
    for roi=find(all(s.significant(:,c,:),3))', fprintf(fid2,'| %d | %s |\n',roi,s.atlasLabels(roi)); end
    fprintf(fid2,'\n');
end
end

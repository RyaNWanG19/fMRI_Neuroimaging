function analysis = compare_wm_equal_size_overlap(outputDir, summaryFile)
% Equal-size category overlap for RMS magnitude and S=sum(T.^2) rankings.
% All 489 ROIs are ranked without a significance filter. Uses saved means/T/S
% and does not rerun inference. Exact score ties use ascending atlas index.
root=fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
if nargin<1 || isempty(outputDir)
    outputDir=fullfile(root,'outputs',['wm_equal_size_overlap_' datestr(now,'yyyymmddTHHMMSS')]);
end
if nargin<2 || isempty(summaryFile)
    summaryFile=fullfile(root,'outputs','wm_load_results_20261005','summary_data.mat');
end
assert(~isfolder(outputDir),'Choose a new output directory.');
a=load(summaryFile,'summary'); s=a.summary; M=489;
assert(isequal(s.models,{'cHRF','cHRFderiv','sHRF'}));
assert(isequal(s.categories,{'body','faces','places','tools'}));
scores=zeros(M,4,3,2); sources=cell(3,1);
for model=1:3
    sources{model}=s.sourceFiles{model}; r=load(sources{model},'results'); r=r.results;
    assert(isequal(r.testMapping.atlas_roi_index,repmat((1:M)',4,1)));
    assert(isequal(string(r.testMapping.category),repelem(string(s.categories(:)),M)));
    assert(isequal(size(r.mean_D),[1956 41])&&isequal(size(r.T_obs),[1956 41]));
    assert(all(r.effectiveN>=2)&&all(isfinite(r.mean_D),'all')&&all(isfinite(r.S_obs)));
    assert(max(abs(sum(r.T_obs.^2,2)-r.S_obs))<1e-8);
    assert(strcmp(r.config.hrfModelName,s.models{model}));
    scores(:,:,model,1)=reshape(sqrt(mean(r.mean_D.^2,2)),M,4);
    scores(:,:,model,2)=reshape(r.S_obs,M,4);
    assert(max(abs(scores(:,:,model,1)-s.rmsMean(:,:,model)),[],'all')<1e-12);
end
mkdir(outputDir); figDir=fullfile(outputDir,'figures'); mkdir(figDir);
analysis=struct('models',{s.models},'categories',{s.categories}, ...
    'rankings',{{'RMS magnitude','Standardized statistic S'}}, ...
    'sourceFiles',{sources},'summaryFile',summaryFile,'K',(1:M)', ...
    'reportK',[50 75 100 150 200 250 300], ...
    'tieRule','Exact score ties broken by ascending atlas index', ...
    'selection','Top K of all 489 ROIs, without a significance filter');
analysis.scores=scores; analysis.rankOrder=zeros(M,4,3,2);
analysis.commonCount=zeros(M,3,2); analysis.commonFraction=zeros(M,3,2);
analysis.pairwiseJaccard=zeros(M,6,3,2);
analysis.boundaryTieSize=zeros(numel(analysis.reportK),4,3,2);
analysis.commonIDs=cell(numel(analysis.reportK),3,2);
pairs=nchoosek(1:4,2); analysis.categoryPairs=pairs;
for method=1:2
    for model=1:3
        ranks=zeros(M,4);
        for c=1:4
            [~,order]=sortrows([-scores(:,c,model,method),(1:M)'],[1 2]);
            analysis.rankOrder(:,c,model,method)=order;
            ranks(order,c)=(1:M)';
            assert(all(diff(scores(order,c,model,method))<=0));
            for j=1:numel(analysis.reportK)
                value=scores(order(analysis.reportK(j)),c,model,method);
                analysis.boundaryTieSize(j,c,model,method)=sum(scores(:,c,model,method)==value);
            end
        end
        worst=max(ranks,[],2);
        counts=cumsum(accumarray(worst,1,[M 1]));
        analysis.commonCount(:,model,method)=counts;
        analysis.commonFraction(:,model,method)=counts./(1:M)';
        for p=1:6
            jointRank=max(ranks(:,pairs(p,:)),[],2);
            intersection=cumsum(accumarray(jointRank,1,[M 1]));
            analysis.pairwiseJaccard(:,p,model,method)=intersection./(2*(1:M)'-intersection);
        end
        for j=1:numel(analysis.reportK)
            K=analysis.reportK(j); chosen=all(ranks<=K,2);
            ids=find(chosen); analysis.commonIDs{j,model,method}=ids;
            % Check against explicit set intersections at every reported K.
            direct=analysis.rankOrder(1:K,1,model,method);
            for c=2:4, direct=intersect(direct,analysis.rankOrder(1:K,c,model,method)); end
            assert(isequal(sort(direct),ids)&&counts(K)==numel(ids));
        end
        assert(counts(M)==M && all(diff(counts)>=0));
        assert(analysis.commonFraction(M,model,method)==1);
        assert(all(analysis.pairwiseJaccard(M,:,model,method)==1));
    end
end
analysis.meanPairwiseJaccard=reshape(mean(analysis.pairwiseJaccard,2),M,3,2);
analysis.uniformIndependentExpectedCount=(1:M)'.^4/M^3;
analysis.uniformIndependentExpectedFraction=((1:M)'/M).^3;
save(fullfile(outputDir,'equal_size_overlap_data.mat'),'analysis','-v7.3');
colors=[22 117 169;212 119 38;106 74 162]/255;
titles={'Common top-ranked regions across categories', ...
    'Fraction of top-ranked regions common to all four categories', ...
    'Pairwise similarity of top-ranked regions across categories'};
ylabs={'Number of regions common to all four','Common regions / K','Mean pairwise Jaccard similarity'};
fileNames={'equal_size_common_counts','equal_size_common_fraction','equal_size_pairwise_similarity'};
values={analysis.commonCount,analysis.commonFraction,analysis.meanPairwiseJaccard};
shown=50:300;
for display=1:3
    f=figure('Visible','off','Color','w','Position',[50 50 1300 620]);
    tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    for method=1:2
        ax=nexttile(tl); hold(ax,'on'); h=gobjects(1,3);
        for model=1:3
            h(model)=plot(ax,shown,values{display}(shown,model,method), ...
                'Color',colors(model,:),'LineWidth',2.3);
            plot(ax,analysis.reportK,values{display}(analysis.reportK,model,method), ...
                'o','MarkerSize',5,'Color',colors(model,:),'HandleVisibility','off');
        end
        title(ax,analysis.rankings{method},'FontSize',19);
        xlabel(ax,'Number of top-ranked regions selected per category (K)');
        ylabel(ax,ylabs{display}); xlim(ax,[50 300]); xticks(ax,[50 100 150 200 250 300]);
        if display==1, ylim(ax,[0 300]); else, ylim(ax,[0 1]); end
        set(ax,'FontSize',13,'Box','off','TickDir','out'); grid(ax,'on');
        if method==1, lg=legend(ax,h,s.models,'Orientation','horizontal'); lg.Layout.Tile='south'; end
    end
    title(tl,titles{display},'FontSize',22);
    exportgraphics(f,fullfile(figDir,[fileNames{display} '.png']),'Resolution',230); close(f);
end
writeReport(analysis,s.atlasLabels,outputDir);
for method=1:2
    fprintf('%s: common counts at K=100: %s\n',analysis.rankings{method}, ...
        mat2str(analysis.commonCount(100,:,method)));
    fprintf('Common fractions at K=50,100,200,300:\n'); disp(analysis.commonFraction([50 100 200 300],:,method));
end
fprintf('Boundary ties greater than 1: %d\n',sum(analysis.boundaryTieSize>1,'all'));
fprintf('EQUAL-SIZE COMPARISON SAVED: %s\n',outputDir);
end

function writeReport(a,labels,out)
fid=fopen(fullfile(out,'equal_size_overlap_report.md'),'w'); assert(fid>0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Equal-sized region-set overlap across WM categories\n\n');
fprintf(fid,['For each category and model, rank all 489 regions and select exactly K. Compare the intersection across body, faces, places and tools. ' ...
    'Selections use no significance filter. RMS ranks the magnitude of each population-average residual contrast curve. ' ...
    'S=sum(T^2) ranks standardized departures from zero across 41 times. S is not a p-value, and differing null distributions can change its relationship to p-value rankings.\n\n' ...
    'Three figures show four-category common counts, common count divided by K, and mean Jaccard similarity across the six category pairs. ' ...
    'Displayed K ranges from 50 to 300; the MAT file contains every K from 1 to 489. All models share the same K and eligible ROI universe.\n\n']);
for method=1:2
    fprintf(fid,'## %s\n\n| K | cHRF common count | cHRFderiv common count | sHRF common count | cHRF fraction | deriv fraction | sHRF fraction |\n|---:|---:|---:|---:|---:|---:|---:|\n',a.rankings{method});
    for K=a.reportK
        fprintf(fid,'| %d | %d | %d | %d | %.3f | %.3f | %.3f |\n',K, ...
            a.commonCount(K,:,method),a.commonFraction(K,:,method));
    end
    fprintf(fid,'\n| K | Mean pairwise Jaccard cHRF | cHRFderiv | sHRF |\n|---:|---:|---:|---:|\n');
    for K=a.reportK, fprintf(fid,'| %d | %.3f | %.3f | %.3f |\n',K,a.meanPairwiseJaccard(K,:,method)); end
    fprintf(fid,'\n');
end
fprintf(fid,['## Interpretation limits\n\n' ...
    '- Equal K controls set size, while preserving each model/category spatial ranking. It does not equalize sampling uncertainty or power.\n' ...
    '- Uniform positive scaling of every ROI RMS in a category preserves its RMS ranking. A lower overlap can reflect relative spatial ordering or noisier ranks, beyond uniform magnitude reduction.\n' ...
    '- Common RMS-ranked membership does not establish matching curve signs or shapes. RMS ignores sign, and S squares the standardized contrast.\n' ...
    '- The comparison is descriptive. No overlap p-values, model-difference tests, subject bootstrap or spatial-null analysis were run.\n' ...
    '- Exact score ties use ascending atlas index. There are %d tied cutoffs at the reported K values across all rankings/categories/models.\n' ...
    '- Permutation p-value ties and their Monte Carlo precision do not enter these rankings. Sampling uncertainty in observed means and standard errors remains.\n' ...
    '- Region identities and rank score order were checked against mapping in the saved source results. Explicit set intersections agree with the cumulative-rank method.\n\n'],sum(a.boundaryTieSize>1,'all'));
fprintf(fid,['Under independently uniform random K-region sets, expected common count is K^4/489^3 and expected common fraction is (K/489)^3. ' ...
    'This is an illustrative size benchmark, not an anatomical significance test. At K=100 the expected count is %.3f. ' ...
    'Spatial dependence and shared subjects invalidate treating it as the study null.\n\n'],a.uniformIndependentExpectedCount(100));
fprintf(fid,'## Sources\n\n');
for k=1:3, fprintf(fid,'- %s: `%s`\n',a.models{k},a.sourceFiles{k}); end
fid2=fopen(fullfile(out,'common_top100_region_lists.md'),'w'); assert(fid2>0);
cleanup2=onCleanup(@() fclose(fid2)); %#ok<NASGU>
fprintf(fid2,'# Regions common across all four top-100 category sets\n\nThese are ranked selections, without a significance requirement.\n\n');
j=find(a.reportK==100);
for method=1:2
    for model=1:3
        ids=a.commonIDs{j,model,method};
        fprintf(fid2,'## %s, %s (%d regions)\n\n| Atlas ROI | Label |\n|---:|---|\n', ...
            a.rankings{method},a.models{model},numel(ids));
        for roi=ids', fprintf(fid2,'| %d | %s |\n',roi,labels(roi)); end
        fprintf(fid2,'\n');
    end
end
end

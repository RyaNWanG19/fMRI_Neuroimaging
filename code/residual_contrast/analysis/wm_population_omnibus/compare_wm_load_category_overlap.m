function overlap = compare_wm_load_category_overlap(summaryFile, outputDir)
% Descriptive overlap across body/faces/places/tools, separately per model.
% Uses saved primary decisions. Does not rerun tests or change multiplicity.
root=fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
if nargin<1 || isempty(summaryFile)
    summaryFile=fullfile(root,'outputs','wm_load_results_20261005','summary_data.mat');
end
if nargin<2 || isempty(outputDir)
    outputDir=fullfile(root,'outputs',['wm_category_overlap_' datestr(now,'yyyymmddTHHMMSS')]);
end
assert(~isfolder(outputDir),'Choose a new output directory.');
a=load(summaryFile,'summary'); s=a.summary;
assert(isequal(size(s.significant),[489 4 3]));
assert(isequal(s.categories,{'body','faces','places','tools'}));
assert(isequal(s.models,{'cHRF','cHRFderiv','sHRF'}));
% Ensure the cached summary still agrees with the actual source result files.
for k=1:3
    r=load(s.sourceFiles{k},'results');
    assert(isequal(reshape(r.results.significant,489,4),s.significant(:,:,k)), ...
        'Cached summary differs from source decisions.');
end
mkdir(outputDir); figures=fullfile(outputDir,'figures'); mkdir(figures);
colors=[22 117 169;212 119 38;106 74 162]/255;
% Category bits are body=1, faces=2, places=4, tools=8. Patterns are exclusive.
patterns=false(16,4);
for c=1:4, patterns(:,c)=logical(bitget((0:15)',c)); end
patternLabels=cell(16,1);
for p=1:16
    if p==1, patternLabels{p}='None';
    elseif p==16, patternLabels{p}='All four';
    else, patternLabels{p}=strjoin(s.categories(patterns(p,:)),' + '); end
end
overlap=struct('models',{s.models},'categories',{s.categories}, ...
    'summarySource',summaryFile,'resultSources',{s.sourceFiles}, ...
    'patternMatrix',patterns,'patternLabels',{patternLabels}, ...
    'note','Descriptive overlap of saved corrected decisions. No overlap significance test.');
overlap.patternCounts=zeros(16,3); overlap.pairwiseCounts=zeros(4,4,3);
overlap.jaccard=zeros(4,4,3); overlap.categoryNumberCounts=zeros(5,3);
overlap.allFourIDs=cell(3,1); overlap.categoryCounts=zeros(4,3);
for k=1:3
    mask=s.significant(:,:,k); code=double(mask)*[1;2;4;8];
    overlap.patternCounts(:,k)=accumarray(code+1,1,[16 1]);
    overlap.categoryNumberCounts(:,k)=accumarray(sum(mask,2)+1,1,[5 1]);
    overlap.categoryCounts(:,k)=sum(mask,1)';
    overlap.allFourIDs{k}=find(all(mask,2));
    for c=1:4
        for d=1:4
            intersection=sum(mask(:,c)&mask(:,d)); union=sum(mask(:,c)|mask(:,d));
            overlap.pairwiseCounts(c,d,k)=intersection;
            if union>0, overlap.jaccard(c,d,k)=intersection/union;
            else, overlap.jaccard(c,d,k)=NaN; end
        end
    end
    assert(sum(overlap.patternCounts(:,k))==489);
    assert(isequal(double(patterns)'*overlap.patternCounts(:,k),overlap.categoryCounts(:,k)));
    assert(isequal(diag(overlap.pairwiseCounts(:,:,k)),overlap.categoryCounts(:,k)));
    assert(overlap.patternCounts(16,k)==numel(overlap.allFourIDs{k}));
    % Put all-four first, followed by other patterns in descending count.
    [~,order]=sort(overlap.patternCounts(2:15,k),'descend'); order=[16;order+1;1];
    vals=overlap.patternCounts(order,k);
    f=figure('Visible','off','Color','w','Position',[50 50 1150 800]); ax=axes(f);
    b=barh(ax,vals,'FaceColor','flat','EdgeColor','none');
    b.CData=repmat([.65 .71 .75],16,1); b.CData(1,:)=colors(k,:);
    set(ax,'YDir','reverse','YTick',1:16,'YTickLabel',patternLabels(order), ...
        'FontSize',13,'Box','off','TickDir','out');
    xlim(ax,[0 max(vals)*1.18+1]); xlabel(ax,'Number of regions');
    text(ax,vals+max(vals)*.02,1:16,string(vals),'FontSize',13,'VerticalAlignment','middle');
    title(ax,sprintf('Category overlap of significant regions: %s',s.models{k}),'FontSize',21);
    subtitle(ax,'Each row is an exact combination, excluding unlisted categories');
    exportgraphics(f,fullfile(figures,[s.models{k} '_exclusive_category_overlap.png']),'Resolution',220); close(f);
end
overlap.allFourCounts=overlap.patternCounts(16,:);
overlap.noneCounts=overlap.patternCounts(1,:);
overlap.anyCounts=489-overlap.noneCounts;
% Show counts and normalized overlap side by side, with matching color scales.
for kind=1:2
    f=figure('Visible','off','Color','w','Position',[50 50 1450 510]);
    tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    if kind==1, matrices=overlap.pairwiseCounts; maxColor=489;
    else, matrices=overlap.jaccard; maxColor=1; end
    for k=1:3
        ax=nexttile(tl); imagesc(ax,matrices(:,:,k)); axis(ax,'image'); clim(ax,[0 maxColor]);
        colormap(ax,parula); set(ax,'XTick',1:4,'YTick',1:4, ...
            'XTickLabel',s.categories,'YTickLabel',s.categories,'FontSize',14);
        title(ax,s.models{k},'FontSize',20); colorbar(ax);
        for c=1:4
            for d=1:4
                v=matrices(c,d,k);
                if kind==1, label=sprintf('%d',v); else, label=sprintf('%.2f',v); end
                if v<.5*maxColor, color='w'; else, color='k'; end
                text(ax,d,c,label,'HorizontalAlignment','center', ...
                    'FontSize',17,'FontWeight','bold','Color',color);
            end
        end
    end
    if kind==1
        title(tl,'Regions significant in both stimulus categories','FontSize',23);
        name='pairwise_category_overlap_counts';
    else
        title(tl,'Similarity of significant-region sets across categories','FontSize',23);
        subtitle(tl,'Jaccard similarity = intersection / union');
        name='pairwise_category_overlap_jaccard';
    end
    exportgraphics(f,fullfile(figures,[name '.png']),'Resolution',220); close(f);
end
f=figure('Visible','off','Color','w','Position',[50 50 1050 650]); ax=axes(f);
b=bar(ax,overlap.categoryNumberCounts,'grouped');
for k=1:3
    b(k).FaceColor=colors(k,:);
    text(ax,b(k).XEndPoints,b(k).YEndPoints+5,string(overlap.categoryNumberCounts(:,k)), ...
        'HorizontalAlignment','center','FontSize',12);
end
set(ax,'XTick',1:5,'XTickLabel',string(0:4),'FontSize',15,'Box','off');
xlabel(ax,'Number of categories in which a region is significant'); ylabel(ax,'Number of regions');
ylim(ax,[0 max(overlap.categoryNumberCounts,[],'all')+35]);
legend(ax,s.models,'Location','northoutside','Orientation','horizontal');
title(ax,'Common significant regions across category','FontSize',21);
exportgraphics(f,fullfile(figures,'number_of_significant_categories.png'),'Resolution',220); close(f);
save(fullfile(outputDir,'category_overlap_data.mat'),'overlap','-v7.3');
fid=fopen(fullfile(outputDir,'category_overlap_report.md'),'w'); assert(fid>0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Significant-region overlap across stimulus categories\n\n');
fprintf(fid,['These comparisons use the saved 2-back minus 0-back omnibus decisions for body, faces, places and tools. ' ...
    'Holm correction covers all 1,956 ROI-category hypotheses separately within each model. The overlap calculations add no new hypothesis tests.\n\n']);
fprintf(fid,'| Model | Body | Faces | Places | Tools | All four | At least one | None |\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for k=1:3
    fprintf(fid,'| %s | %d | %d | %d | %d | %d | %d | %d |\n',s.models{k}, ...
        overlap.categoryCounts(:,k)',overlap.allFourCounts(k),overlap.anyCounts(k),overlap.noneCounts(k));
end
fprintf(fid,['\n“All four” means the same atlas region passes the load-difference test for every stimulus category under that model. ' ...
    'It supports broadly recurring load-dependent residual structure across categories. It does not establish equal curve shapes, equal directions or a common physiological cause. ' ...
    'Threshold patterns can change even when underlying effects are similar. Shared subjects and stimuli make these descriptive comparisons dependent.\n\n']);
for k=1:3
    fprintf(fid,'## %s\n\n',s.models{k});
    fprintf(fid,'%d regions (%.1f%% of all 489 regions; %.1f%% of the union) pass in all four categories.\n\n', ...
        overlap.allFourCounts(k),100*overlap.allFourCounts(k)/489,100*overlap.allFourCounts(k)/overlap.anyCounts(k));
    fprintf(fid,'### Pairwise overlap\n\n| Categories | Intersection | Union | Jaccard similarity |\n|---|---:|---:|---:|\n');
    for c=1:3
        for d=c+1:4
            mask=s.significant(:,:,k); union=sum(mask(:,c)|mask(:,d));
            fprintf(fid,'| %s + %s | %d | %d | %.3f |\n',s.categories{c},s.categories{d}, ...
                overlap.pairwiseCounts(c,d,k),union,overlap.jaccard(c,d,k));
        end
    end
    fprintf(fid,'\n### Exact category combinations\n\n| Combination | Regions |\n|---|---:|\n');
    for p=16:-1:1, fprintf(fid,'| %s | %d |\n',patternLabels{p},overlap.patternCounts(p,k)); end
    fprintf(fid,'\n### Regions significant in all four categories\n\n| Atlas ROI | Label |\n|---:|---|\n');
    for roi=overlap.allFourIDs{k}', fprintf(fid,'| %d | %s |\n',roi,s.atlasLabels(roi)); end
    fprintf(fid,'\n');
end
fprintf(fid,['## Figure interpretation\n\n' ...
    '- Pairwise count matrices: diagonal values are category totals; off-diagonal values count regions significant in both categories. A region can contribute to several pairs.\n' ...
    '- Jaccard matrices: intersection divided by union. Zero means no overlap and one means identical sets. This normalizes for differing set sizes.\n' ...
    '- Exact-combination bar charts: every region contributes to exactly one row. “body + faces” excludes places and tools.\n' ...
    '- Breadth chart: regions significant in exactly zero, one, two, three or four categories.\n\n' ...
    'The original Monte Carlo precision qualification still applies to the saved classifications. No bootstrap stability or overlap significance analysis is claimed. ' ...
    'Smaller sets can have lower Jaccard overlap due to thresholding and power. This does not independently establish greater biological specificity or better fit.\n\n']);
fprintf('All-four counts: %s\nAny-category counts: %s\nNone counts: %s\n', ...
    mat2str(overlap.allFourCounts),mat2str(overlap.anyCounts),mat2str(overlap.noneCounts));
fprintf('Category overlap output: %s\n',outputDir);
end

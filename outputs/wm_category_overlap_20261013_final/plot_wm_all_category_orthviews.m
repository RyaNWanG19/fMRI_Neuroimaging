function roiIndices = plot_wm_all_category_orthviews(model, xyz, outputFile, summaryFile)
% Map ROIs significant in all four category-specific 2-back minus 0-back tests.
% Example (CANlab and SPM on path):
%   plot_wm_all_category_orthviews('cHRF',[0 -20 20],'cHRF_all_four.png');
% Whole-parcel membership overlay, not a voxelwise or effect-direction map.
root=fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
if nargin<1 || isempty(model), model='cHRF'; end
if nargin<2 || isempty(xyz), xyz=[0 -20 20]; end
if nargin<3, outputFile=''; end
if nargin<4 || isempty(summaryFile)
    summaryFile=fullfile(root,'outputs','wm_load_results_20261005','summary_data.mat');
end
validateattributes(xyz,{'numeric'},{'real','finite','numel',3});
if ~isempty(outputFile), assert(~isfile(outputFile),'Refusing to overwrite existing output.'); end
a=load(summaryFile,'summary'); s=a.summary;
k=find(strcmpi(model,s.models),1); assert(~isempty(k),'Unknown HRF model.');
model=s.models{k};
r=load(s.sourceFiles{k},'results'); r=r.results;
assert(isequal(r.testMapping.atlas_roi_index,repmat((1:489)',4,1)));
assert(isequal(string(r.testMapping.category),repelem(string(s.categories(:)),489)));
assert(isequal(reshape(r.significant,489,4),s.significant(:,:,k)), ...
    'Cached summary differs from source significance decisions.');
selected=all(reshape(r.significant,489,4),2);
roiIndices=find(selected);
assert(~isempty(roiIndices),'No regions pass in all four categories.');
atlasObj=load_atlas('canlab2018');
assert(isequal(string(atlasObj.labels(:)),s.atlasLabels(:)), ...
    'Current atlas label order differs from summary.');
sub=select_atlas_subset(atlasObj,roiIndices,'deterministic'); regions=atlas2region(sub);
colors=[22 117 169;212 119 38;106 74 162]/255;
orthviews(regions,repmat({colors(k,:)},1,numel(regions)),'solid');
spm_orthviews('Reposition',xyz(:)); drawnow;
fig=spm_figure('FindWin','Graphics');
assert(~isempty(fig)&&isgraphics(fig),'No SPM Graphics figure found.');
% Correct static CANlab coordinate labels after setting the common slices.
txt=findall(fig,'Type','text');
for j=1:numel(txt)
    value=get(txt(j),'String');
    if ischar(value)
        tok=regexp(value,'^\s*([xyz])\s*=','tokens','once');
        if ~isempty(tok)
            set(txt(j),'String',sprintf('%s = %g',tok{1},xyz('xyz'==tok{1})));
        end
    end
end
annotation(fig,'textbox',[.03 .90 .94 .09],'String', ...
    sprintf('Significant in all four stimulus categories\n%s (n = %d regions)',model,numel(roiIndices)), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',14,'FontWeight','bold', ...
    'Interpreter','none','Color','k','BackgroundColor','w');
annotation(fig,'textbox',[.03 .01 .94 .075],'String', ...
    sprintf(['Memory load difference: body, faces, places and tools\n' ...
    'MNI [%g %g %g] mm. Color shows region membership.\n' ...
    'Displayed slices show only part of the significant-region set.'],xyz), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9, ...
    'Color','k','BackgroundColor','w');
drawnow;
if ~isempty(outputFile), exportgraphics(fig,outputFile,'Resolution',250); end
fprintf('%s: %d regions significant in all four categories.\n',model,numel(roiIndices));
end

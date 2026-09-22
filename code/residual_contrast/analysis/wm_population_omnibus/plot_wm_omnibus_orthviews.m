function roiIndices = plot_wm_omnibus_orthviews(mapName, xyz, outputFile, baseDir)
%PLOT_WM_OMNIBUS_ORTHVIEWS Anatomical display of saved significant ROI sets.
% Requires CANlab and SPM. Run from MATLAB with those toolboxes on the path.
% Examples:
%   plot_wm_omnibus_orthviews('shared');
%   plot_wm_omnibus_orthviews('cHRF', [0 -20 20]);
%   plot_wm_omnibus_orthviews('sHRF', [0 -20 20], 'sHRF_map.png');
% Colors indicate membership, not effect size, direction, or voxelwise tests.
% xyz is the common MNI crosshair position (mm); change it interactively or
% with spm_orthviews('Reposition', [x y z]). Three slices are not whole-brain coverage.
% ROI mapping follows the extraction script's canlab2018 index convention;
% ordinal IDs cannot independently establish historical atlas-version identity.

if nargin < 1 || isempty(mapName), mapName = 'shared'; end
if nargin < 2 || isempty(xyz), xyz = [0 -20 20]; end
if nargin < 3, outputFile = ''; end
if nargin < 4 || isempty(baseDir)
    root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
    baseDir = fullfile(root, 'data', 'task_residual', ...
        'jhpce_outputs_wm_population_omnibus_full', 'Body_LoadDiff_vs_Face_LoadDiff');
end
validateattributes(xyz, {'numeric'}, {'real','finite','numel',3});
if ~isempty(outputFile) && isfile(outputFile)
    error('Refusing to overwrite existing file: %s', outputFile);
end
models = {'cHRF','cHRFderiv','sHRF'};
colors = [22 117 169; 212 119 38; 106 74 162; 35 140 110]/255;
atlasObj = load_atlas('canlab2018');
nROI = numel(atlasObj.labels);
assert(nROI == 489, 'Expected the 489-region canlab2018 atlas.');
sig = false(nROI,3);
for k = 1:3
    files = dir(fullfile(baseDir,models{k},'**','temporal_omnibus_results.mat'));
    assert(numel(files)==1, 'Expected exactly one results file for %s.', models{k});
    a = load(fullfile(files(1).folder,files(1).name),'results');
    r = a.results;
    assert(isequal(string(r.roiIDs(:)),string((1:nROI)')), 'Unexpected ROI index mapping.');
    assert(strcmpi(r.config.hrfModelName,models{k}), 'Model label mismatch.');
    assert(strcmpi(r.config.contrastName,'Body_LoadDiff_vs_Face_LoadDiff'), ...
        'Unexpected contrast.');
    assert(isequal(logical(r.significant(:)),r.p_holm(:)<=r.config.alpha), ...
        'Saved significance labels do not match adjusted p-values.');
    sig(:,k) = r.significant(:);
end
if strcmpi(mapName,'shared')
    selected = all(sig,2); color = colors(4,:);
    label = 'Significant regions shared by all three HRF models';
else
    k = find(strcmpi(mapName,models),1);
    assert(~isempty(k), 'mapName must be cHRF, cHRFderiv, sHRF, or shared.');
    selected = sig(:,k); color = colors(k,:);
    label = sprintf('Significant regions: %s',models{k});
end
roiIndices = find(selected);
assert(~isempty(roiIndices), 'No significant regions in the requested set.');

% Deterministic atlas boundaries prevent probabilistic parcel expansion.
subset = select_atlas_subset(atlasObj,roiIndices,'deterministic');
regions = atlas2region(subset);
orthviews(regions,repmat({color},1,numel(regions)),'solid');
spm_orthviews('Reposition',xyz(:));
drawnow;
fig = spm_figure('FindWin','Graphics');
assert(~isempty(fig) && isgraphics(fig), 'Cannot find the SPM Graphics window.');
% CANlab may leave static coordinate labels from its initial peak location.
% Update those labels after repositioning so they match the displayed slices.
txt = findall(fig,'Type','text');
for j = 1:numel(txt)
    value = get(txt(j),'String');
    if ischar(value)
        token = regexp(value,'^\s*([xyz])\s*=','tokens','once');
        if ~isempty(token)
            axisIndex = find('xyz'==token{1});
            set(txt(j),'String',sprintf('%s = %g',token{1},xyz(axisIndex)));
        end
    end
end
set(fig,'Name',label);
annotation(fig,'textbox',[0.03 0.91 0.94 0.07], ...
    'String',sprintf('%s (n = %d)',label,numel(roiIndices)), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',13, ...
    'FontWeight','bold','Interpreter','none','Color','k','BackgroundColor','w');
annotation(fig,'textbox',[0.03 0.01 0.94 0.065], ...
    'String',sprintf(['Body load difference minus face load difference | MNI [%g %g %g] mm\n' ...
    'Color = ROI significance membership; displayed slices do not show all regions.'],xyz), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9, ...
    'Interpreter','none','Color','k','BackgroundColor','w');
drawnow;
if ~isempty(outputFile)
    exportgraphics(fig,outputFile,'Resolution',250);
end
fprintf('%s: %d of %d ROIs.\n',label,numel(roiIndices),nROI);
end

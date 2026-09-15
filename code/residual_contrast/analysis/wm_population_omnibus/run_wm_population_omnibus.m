function results = run_wm_population_omnibus(cfg)
%RUN_WM_POPULATION_OMNIBUS Project driver for WM residual contrasts.
%
% The source Results arrays are [ROI x time x subject x condition]. This
% driver constructs, without temporal averaging,
%   D = (2-back category A - 0-back category A) ...
%       - (2-back category B - 0-back category B),
% then temporal_omnibus_signflip maps [ROI x time x subject] to canonical
% [subject x ROI x time].
%
% Full inference requires real IDs plus verified independent rows:
%   cfg.subjectIDs = ...; % aligned with Results dimension 3
%   cfg.resampling.independentSubjectsVerified = true;
%   cfg.resampling.nullSymmetryVerified = true;
%   cfg.resampling.verificationNote = 'How independence and symmetry were justified';
%
% Set cfg.smokeTest=true for a non-inferential small run.

if nargin < 1
    cfg = struct();
end
if ~isstruct(cfg)
    error('cfg must be a struct.');
end
moduleDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(moduleDir, '..', '..', '..', '..');
cfg = set_default(cfg, 'repoRoot', repoRoot);
cfg = set_default(cfg, 'hrfModelName', 'cHRF');
cfg = set_default(cfg, 'contrastName', 'Body_LoadDiff_vs_Face_LoadDiff');
cfg = set_default(cfg, 'smokeTest', false);
cfg = set_default(cfg, 'saveOutputs', true);
cfg = set_default(cfg, 'TR', 0.72);

conditionMap = verified_condition_map(cfg.repoRoot);
[dataFile, dataVariable] = resolve_data(cfg);
payload = load(dataFile, dataVariable);
residuals = payload.(dataVariable);
sz = size(residuals);
sz(end+1:4) = 1;
if ~(sz(1) == 489 && sz(2) == 41 && sz(4) == 8)
    error('Expected Results [489 x 41 x subject x 8], found [%s].', num2str(sz));
end
nSubject = sz(3);

[pair, contrastLabel] = resolve_contrast(cfg.contrastName, conditionMap);
% pair = [2back_A, 0back_A, 2back_B, 0back_B]
D = (residuals(:, :, :, pair(1)) - residuals(:, :, :, pair(2))) - ...
    (residuals(:, :, :, pair(3)) - residuals(:, :, :, pair(4)));
clear residuals payload

if isfield(cfg, 'subjectIDs') && ~isempty(cfg.subjectIDs)
    subjectIDs = cfg.subjectIDs;
    subjectIDSource = 'cfg.subjectIDs';
else
    subjectIDs = (1:nSubject)';
    subjectIDSource = ['ordinal row indices; source MAT file contains no real HCP IDs'];
end
if numel(subjectIDs) ~= nSubject
    error('cfg.subjectIDs has %d entries; data has %d subjects.', numel(subjectIDs), nSubject);
end

if isfield(cfg, 'roiIDs') && ~isempty(cfg.roiIDs)
    roiIDs = cfg.roiIDs;
else
    roiIDs = (1:sz(1))';
end
if isfield(cfg, 'roiNames') && ~isempty(cfg.roiNames)
    roiNames = cfg.roiNames;
    roiNameSource = 'cfg.roiNames';
else
    [roiNames, roiNameSource] = resolve_roi_names(sz(1));
end

% No time vector is stored in the MAT files. Fit_HRF_Residuals_HCP.m calls
% Fit_sFIR_all with TR=0.72 and existing project plots explicitly use
% (0:T-1)*TR for these 41 coefficients, so reconstruct that established
% extraction/plotting convention rather than guessing from the array alone.
timeSec = (0:(sz(2)-1))' .* cfg.TR;

cfg.inputDimensionOrder = 'roi-time-subject';
cfg.contrastName = contrastLabel;
cfg = set_default(cfg, 'seed', 20260914);
cfg = set_default(cfg, 'alpha', 0.05);
if cfg.smokeTest
    cfg.runMode = 'smoke-test';
    cfg = set_default(cfg, 'nPerm', 999);
else
    cfg.runMode = 'inference';
    cfg = set_default(cfg, 'nPerm', 99999);
    if strcmp(subjectIDSource(1:7), 'ordinal')
        error('run_wm_population_omnibus:RealSubjectIDsRequired', ...
            ['The HRF-model MAT files do not store real subject IDs. Supply ', ...
             'cfg.subjectIDs aligned with Results dimension 3 before inference.']);
    end
end
if ~isfield(cfg, 'resampling') || ~isstruct(cfg.resampling)
    cfg.resampling = struct();
end
cfg.resampling = set_default(cfg.resampling, 'scheme', 'independent-subject-signflip');
cfg.resampling = set_default(cfg.resampling, 'independentSubjectsVerified', false);
cfg.resampling = set_default(cfg.resampling, 'nullSymmetryVerified', false);
cfg.resampling = set_default(cfg.resampling, 'verificationNote', '');

if cfg.saveOutputs && (~isfield(cfg, 'outputDir') || isempty(cfg.outputDir))
    stamp = datestr(now, 'yyyymmddTHHMMSSFFF');
    modeSuffix = '';
    if cfg.smokeTest, modeSuffix = '_SMOKE'; end
    cfg.outputDir = fullfile(cfg.repoRoot, 'data', 'task_residual', ...
        'wm_population_omnibus', safe_name(contrastLabel), cfg.hrfModelName, ...
        [stamp modeSuffix]);
end

cfg.dataMetadata = struct( ...
    'dataFile', dataFile, ...
    'dataVariable', dataVariable, ...
    'sourceSize', sz, ...
    'sourceDimensionOrder', 'ROI x time x subject x condition', ...
    'DSourceDimensionOrder', 'ROI x time x subject', ...
    'DCanonicalPermutation', [3 1 2], ...
    'subjectIDSource', subjectIDSource, ...
    'roiNameSource', roiNameSource, ...
    'conditionIndices', pair, ...
    'conditionLabels', {conditionMap.conditionLabels}, ...
    'conditionDefinitionSources', {conditionMap.sourceFiles}, ...
    'timeVectorSource', ['Reconstructed as (0:40)*0.72 s from Fit_HRF_Residuals_HCP.m ' ...
        'and established WM plotting code.']);

results = temporal_omnibus_signflip(D, subjectIDs, roiIDs, roiNames, timeSec, cfg);
end

function map = verified_condition_map(repoRoot)
sourceFiles = { ...
    fullfile(repoRoot, 'functional_clustering', 'build_wm_contrast.m'), ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'data_extract', 'GenerateOnsetsHCP.m'), ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'data_extract', 'PredictContrast.m'), ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'core', 'PredictContrast_WM.m')};
for i = 1:numel(sourceFiles)
    if ~isfile(sourceFiles{i})
        error('Required condition-definition source is missing: %s', sourceFiles{i});
    end
end
txt = lower(strjoin(cellfun(@fileread, sourceFiles, 'UniformOutput', false), newline));
required = {'0bk_body','0bk_faces','0bk_places','0bk_tools', ...
    '2bk_body','2bk_faces','2bk_places','2bk_tools'};
if ~all(cellfun(@(x) contains(txt, x), required)) || ...
        ~contains(txt, 'idx_0 = [1 2 3 4]') || ~contains(txt, 'idx_2 = [5 6 7 8]')
    error('Could not verify WM condition order from repository source files.');
end
map.sourceFiles = sourceFiles;
map.conditionLabels = required;
map.pairs = [5 1 6 2; 5 1 7 3; 5 1 8 4; 6 2 7 3; 6 2 8 4; 7 3 8 4];
map.labels = {'Body_LoadDiff_vs_Face_LoadDiff'; ...
    'Body_LoadDiff_vs_Place_LoadDiff'; 'Body_LoadDiff_vs_Tool_LoadDiff'; ...
    'Face_LoadDiff_vs_Place_LoadDiff'; 'Face_LoadDiff_vs_Tool_LoadDiff'; ...
    'Place_LoadDiff_vs_Tool_LoadDiff'};
end

function [pair, label] = resolve_contrast(requested, map)
idx = find(strcmpi(char(requested), map.labels), 1);
if isempty(idx)
    error('Unknown WM contrast %s. Allowed values:\n  %s', char(requested), ...
        strjoin(map.labels, '\n  '));
end
pair = map.pairs(idx, :);
label = map.labels{idx};
end

function [dataFile, dataVariable] = resolve_data(cfg)
if isfield(cfg, 'dataFile') && ~isempty(cfg.dataFile)
    dataFile = cfg.dataFile;
else
    dataDir = fullfile(cfg.repoRoot, 'data', 'task_residual');
    switch lower(char(cfg.hrfModelName))
        case 'chrf'
            fileName = 'WMcHRF.mat';
        case 'chrfderiv'
            fileName = 'WMcHRFderiv.mat';
        case 'shrf'
            fileName = 'WMsHRF.mat';
        otherwise
            error('hrfModelName must be cHRF, cHRFderiv, or sHRF.');
    end
    dataFile = fullfile(dataDir, fileName);
end
if ~isfile(dataFile)
    error('WM data file not found: %s', dataFile);
end
dataVariable = 'Results';
vars = whos('-file', dataFile);
if ~ismember(dataVariable, {vars.name})
    error('%s does not contain variable Results.', dataFile);
end
end

function [roiNames, source] = resolve_roi_names(nROI)
roiNames = [];
source = 'generic Parcel_### labels';
if exist('load_atlas', 'file') ~= 2
    return
end
try
    atlas = load_atlas('canlab2018');
    if isobject(atlas) && isprop(atlas, 'labels')
        labels = atlas.labels;
    elseif isstruct(atlas) && isfield(atlas, 'labels')
        labels = atlas.labels;
    else
        labels = [];
    end
    if numel(labels) == nROI
        roiNames = labels(:);
        source = 'load_atlas(''canlab2018'').labels';
    end
catch ME
    warning('run_wm_population_omnibus:AtlasLabelsUnavailable', ...
        'Could not load canlab2018 ROI labels (%s); using generic labels.', ME.message);
end
end

function s = safe_name(s)
s = regexprep(char(s), '[^A-Za-z0-9_+-]+', '_');
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

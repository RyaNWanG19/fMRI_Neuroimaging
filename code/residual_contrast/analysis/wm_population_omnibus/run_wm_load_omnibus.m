function results = run_wm_load_omnibus(cfg)
% Four category-specific load tests, jointly Holm-corrected within one model.
% Uses existing sign-flip engine with 1956 ROI-category hypotheses. Each
% subject's sign is shared across ALL categories, ROIs and times. Missingness
% is fixed independently per ROI-category across all 41 times. Optional
% centered bootstrap uses the same 1956-test correction family.
% This tests differential load mismatch, NOT absolute residuals per condition.
if nargin<1, cfg=struct(); end
moduleDir=fileparts(mfilename('fullpath'));
cfg=defaults(cfg,'repoRoot',fileparts(fileparts(fileparts(fileparts(moduleDir)))));
cfg=defaults(cfg,'hrfModelName','cHRF');
models={'cHRF','cHRFderiv','sHRF'};
k=find(strcmpi(cfg.hrfModelName,models),1);
assert(~isempty(k),'Unknown HRF model.'); cfg.hrfModelName=models{k};
cfg=defaults(cfg,'dataFile',fullfile(cfg.repoRoot,'data','task_residual',['WM' models{k} '.mat']));
cfg=defaults(cfg,'saveOutputs',true);
cfg=defaults(cfg,'smokeTest',false);
cfg=defaults(cfg,'TR',0.72);
validateattributes(cfg.TR,{'numeric'},{'scalar','finite','positive'});
if cfg.smokeTest
    cfg.runMode='smoke-test'; cfg=defaults(cfg,'nPerm',999);
else
    cfg.runMode='inference'; cfg=defaults(cfg,'nPerm',99999);
end
cfg.contrastName='AllFourCategory_LoadDiff';
cfg.inputDimensionOrder='roi-time-subject';
cfg=defaults(cfg,'outputDir',fullfile(cfg.repoRoot,'data','task_residual', ...
    'wm_load_omnibus',cfg.hrfModelName,datestr(now,'yyyymmddTHHMMSSFFF')));
saveOutputs=cfg.saveOutputs;
matFile=fullfile(cfg.outputDir,'temporal_omnibus_results.mat');
csvFile=fullfile(cfg.outputDir,'temporal_omnibus_roi_summary.csv');
bootFile=fullfile(cfg.outputDir,'omnibus_sensitivity.csv');
if saveOutputs
    assert(~isfile(matFile)&&~isfile(csvFile)&&~isfile(bootFile), ...
        'Refusing to overwrite existing analysis outputs. Choose a new outputDir.');
end
a=load(cfg.dataFile,'Results');
assert(isfield(a,'Results'),'Data file must contain Results.');
sz=size(a.Results); sz(end+1:4)=1;
assert(sz(1)==489 && sz(2)==41 && sz(4)==8,'Expected 489 x 41 x subject x 8.');
[D,mapping]=build_wm_load_differences(a.Results); clear a
n=sz(3);
cfg=defaults(cfg,'subjectIDs',compose('row_%03d',(1:n)'));
assert(numel(cfg.subjectIDs)==n,'Subject ID count does not match source data.');
cfg.dataMetadata=struct('dataFile',cfg.dataFile,'dataVariable','Results', ...
    'sourceSize',sz,'conditionOrder',{{'0bk_body','0bk_faces','0bk_places','0bk_tools', ...
    '2bk_body','2bk_faces','2bk_places','2bk_tools'}}, ...
    'conditionPairs',[5 1;6 2;7 3;8 4], ...
    'mappingConvention','canlab2018 ROI index; category-major flattened tests');
% Suppress generic per-ROI persistence/plots until category metadata is attached.
cfg.saveOutputs=false; cfg.makePlots=false;
results=temporal_omnibus_signflip(D,cfg.subjectIDs,mapping.test_id, ...
    mapping.test_id,(0:40)'*cfg.TR,cfg);
results.config.saveOutputs=saveOutputs;
results.analysisType='category-specific 2-back minus 0-back residual omnibus';
results.testMapping=mapping;
results.correctionFamily=sprintf(['All 489 ROIs x 4 categories = 1956 tests within %s. ' ...
    'No correction across HRF models. Untestable hypotheses retained with p=1.'],models{k});
results.summaryTable.atlas_roi_index=mapping.atlas_roi_index;
results.summaryTable.category=mapping.category;
results.summaryTable.correction_family=repmat(string(results.correctionFamily),height(mapping),1);
results.categoryCounts=table(["body";"faces";"places";"tools"], ...
    sum(reshape(results.significant,489,4),1)', ...
    'VariableNames',{'category','significant_regions'});
if isfield(results,'sensitivity')
    results.sensitivity.table.atlas_roi_index=mapping.atlas_roi_index;
    results.sensitivity.table.category=mapping.category;
end
if saveOutputs
    if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
    results.outputFiles.mat=matFile; results.outputFiles.csv=csvFile;
    writetable(results.summaryTable,csvFile);
    if isfield(results,'sensitivity')
        results.outputFiles.sensitivityCSV=bootFile;
        writetable(results.sensitivity.table,bootFile);
    end
    save(matFile,'results','-v7.3');
end
disp(results.categoryCounts);
end

function cfg=defaults(cfg,name,value)
if ~isfield(cfg,name)||isempty(cfg.(name)), cfg.(name)=value; end
end

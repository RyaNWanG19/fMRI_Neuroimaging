function test_wm_load_omnibus()
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','analysis','wm_population_omnibus'));
rng(704);
x=randn(3,5,7,8);
[d,map]=build_wm_load_differences(x);
assert(isequal(size(d),[12 5 7]));
for c=1:4
    assert(isequal(d((c-1)*3+(1:3),:,:),x(:,:,:,c+4)-x(:,:,:,c)));
end
assert(isequal(map.atlas_roi_index,repmat((1:3)',4,1)));
% Synthetic full-shape driver run: no real-data inference or output overwrites.
scratch=tempname; mkdir(scratch); cleanup=onCleanup(@() rmdir(scratch,'s')); %#ok<NASGU>
Results=randn(489,41,7,8); Results(2,3,1,5)=NaN;
save(fullfile(scratch,'source.mat'),'Results');
cfg=struct('dataFile',fullfile(scratch,'source.mat'),'outputDir',fullfile(scratch,'out'), ...
    'smokeTest',true,'nPerm',17,'seed',129,'verbose',false,'saveOutputs',true, ...
    'runDiagnostics',true,'nBootstrap',7);
r=run_wm_load_omnibus(cfg);
assert(numel(r.p_raw)==1956);
assert(isequal(r.p_holm,holm_adjust(r.p_raw)));
assert(r.effectiveN(2)==6 && r.effectiveN(489+2)==7);
for c=1:4
    roi=4; idx=(c-1)*489+roi;
    z=squeeze(Results(roi,:,:,c+4)-Results(roi,:,:,c))';
    assert(max(abs(r.mean_D(idx,:)-mean(z,1)))<1e-12);
    assert(max(abs(r.T_obs(idx,:)-mean(z,1)./(std(z,0,1)/sqrt(7))))<1e-12);
end
assert(isequal(r.sensitivity.table.bootstrap_p_holm, ...
    holm_adjust(r.sensitivity.table.bootstrap_p_raw)));
a=load(r.outputFiles.mat,'results'); assert(isequal(a.results.testMapping,r.testMapping));
assert(isfile(r.outputFiles.csv));
blocked=false;
try, run_wm_load_omnibus(cfg); catch ME, blocked=contains(ME.message,'overwrite'); end
assert(blocked);
cfg.saveOutputs=false; cfg.smokeTest=false;
blocked=false;
try, run_wm_load_omnibus(cfg); catch ME, blocked=strcmp(ME.identifier,'temporal_omnibus_signflip:IndependenceUnresolved'); end
assert(blocked);
fprintf('All four-load-difference tests passed.\n');
end

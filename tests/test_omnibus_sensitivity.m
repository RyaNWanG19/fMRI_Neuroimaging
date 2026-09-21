function test_omnibus_sensitivity()
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','analysis','wm_population_omnibus'));
rng(42); D = exp(randn(13,3,4)); D(1,2,1)=NaN;
c = struct('runMode','inference','inputDimensionOrder','srt', ...
    'nPerm',31,'nBootstrap',37,'seed',13,'batchSize',11, ...
    'runDiagnostics',true,'saveOutputs',false,'verbose',false, ...
    'storePermutationStatistics',true);
c.resampling = struct('independentSubjectsVerified',true, ...
    'nullSymmetryAssumed',true,'verificationNote','synthetic independent draws');
r = temporal_omnibus_signflip(D,(1:13)',[],[],1:4,c);
assert(r.inferenceStatus == "conditional_on_independence_and_null_symmetry");
assert(~r.config.resampling.nullSymmetryVerified);
% Independently reconstruct each bootstrap draw, direct mean/sample SD.
for roi = 1:3
    x = reshape(D(r.subjectMasks(:,roi),roi,:),[],4);
    n = size(x,1); y = x-mean(x,1);
    stream = RandStream('mt19937ar','Seed',c.seed+104729+roi);
    direct = zeros(c.nBootstrap,1);
    for b = 1:c.nBootstrap
        xb = y(randi(stream,n,n,1),:);
        direct(b) = sum((mean(xb,1)./(std(xb,0,1)/sqrt(n))).^2);
    end
    assert(max(abs(direct-r.sensitivity.S_boot(:,roi)))<1e-8);
    assert(r.sensitivity.table.bootstrap_p_raw(roi)== ...
        (1+sum(direct>=r.S_obs(roi)))/(c.nBootstrap+1));
end
c.batchSize=7;
again=temporal_omnibus_signflip(D,(1:13)',[],[],1:4,c);
assert(max(abs(again.sensitivity.S_boot-r.sensitivity.S_boot),[],'all')<1e-8);
% Degenerate observations must not become bootstrap discoveries.
D(:,1,:)=2;
d=temporal_omnibus_signflip(D,(1:13)',[],[],1:4,c);
assert(d.sensitivity.table.bootstrap_p_raw(1)==1);
assert(contains(d.sensitivity.table.bootstrap_status(1),'degenerate'));
fprintf('Sensitivity direct-loop, asymmetric-input, reproducibility and degeneracy checks passed.\n');
end

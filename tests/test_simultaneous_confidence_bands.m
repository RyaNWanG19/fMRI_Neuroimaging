function test_simultaneous_confidence_bands()
% Numerical reference, joint dependence, degeneracy and persistence checks.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','analysis','wm_population_omnibus'));
cfg=struct('runMode','verification','inputDimensionOrder','subject-roi-time', ...
    'nPerm',7,'seed',991,'alpha',0.05,'batchSize',23,'roiBatchSize',2, ...
    'saveOutputs',false,'makePlots',false,'verbose',false, ...
    'runSimultaneousBands',true,'scb',struct('nBootstrap',399,'minSubjects',2, ...
    'method','gaussian-multiplier'));
rng(1701);
D=randn(32,3,7);
D(:,2,:)=D(:,1,:); % exact duplication must not count as independent tests
D(1,3,2)=NaN;
r=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfg);
b=r.simultaneousBands;
assert(isequal(size(b.lower),[3 7]) && b.plannedCellCount==21);
assert(all(b.eligible,'all') && ~b.inferenceValid);
assert(max(abs(b.standardError-r.sd_D./sqrt(r.effectiveN)),[],'all')<1e-12);

% Direct whole-family reference with shared subject multipliers and masks.
stream=RandStream('mt19937ar','Seed',b.config.seed);
stream.NormalTransform='Inversion';
G=randn(stream,32,b.config.nBootstrap)';
reference=zeros(b.config.nBootstrap,1);
for j=1:3
    keep=r.subjectMasks(:,j);
    x=reshape(D(keep,j,:),nnz(keep),7);
    e=x-mean(x,1);
    z=G(:,keep)*(e./vecnorm(e,2,1));
    reference=max(reference,max(abs(z),[],2));
end
assert(max(abs(b.bootstrapMaxima-reference))<1e-12);
ordered=sort(reference);
assert(b.criticalValue==ordered(ceil((399+1)*0.95)));
assert(max(abs(b.lower-(r.mean_D-b.criticalValue*b.standardError)),[],'all')<1e-12);
assert(isequal(b.positive,b.lower>0) && isequal(b.negative,b.upper<0));

% Changing batching never changes RNG, family calibration or omnibus output.
cfg2=cfg; cfg2.scb.batchSize=71; cfg2.scb.testBatchSize=1;
r2=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfg2);
assert(max(abs(b.bootstrapMaxima-r2.simultaneousBands.bootstrapMaxima))<1e-12);
assert(abs(b.criticalValue-r2.simultaneousBands.criticalValue)<1e-12);
cfgOff=cfg; cfgOff.runSimultaneousBands=false;
off=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfgOff);
assert(isequal(r.p_raw,off.p_raw) && isequal(r.p_holm,off.p_holm));
assert(~isfield(off,'simultaneousBands'));

% Translation does not change uncertainty; sign/scale equivariance.
shifted=temporal_omnibus_signflip(D+4,(1:32)',[],[],(0:6)',cfg);
assert(abs(b.criticalValue-shifted.simultaneousBands.criticalValue)<1e-12);
assert(max(abs(shifted.simultaneousBands.lower-(b.lower+4)),[],'all')<1e-12);
neg=temporal_omnibus_signflip(-2*D,(1:32)',[],[],(0:6)',cfg);
assert(max(abs(neg.simultaneousBands.lower+2*b.upper),[],'all')<1e-12);
assert(isequal(neg.simultaneousBands.negative,b.positive));

% Perfect time/test dependence contributes only one effective coordinate.
v=randn(32,1);
duplicated=repmat(reshape(v,32,1,1),1,5,7);
dup=temporal_omnibus_signflip(duplicated,(1:32)',[],[],(0:6)',cfg);
single=temporal_omnibus_signflip(reshape(v,32,1,1),(1:32)',[],[],0,cfg);
assert(abs(dup.simultaneousBands.criticalValue-single.simultaneousBands.criticalValue)<1e-12);

% Never report certainty from constant cells or too few complete subjects.
q=D; q(:,1,2)=5; q(:,2,:)=0; q(3:end,3,:)=NaN;
cfgGuard=cfg; cfgGuard.scb.minSubjects=10;
bad=temporal_omnibus_signflip(q,(1:32)',[],[],(0:6)',cfgGuard);
bb=bad.simultaneousBands;
assert(~bb.eligible(1,2) && ~any(bb.eligible(2:3,:),'all'));
assert(all(bb.lower(~bb.eligible)==-Inf) && all(bb.upper(~bb.eligible)==Inf));
assert(~any(bb.positive(~bb.eligible)) && ~any(bb.negative(~bb.eligible)));
assert(height(bb.table)==21 && bb.table.mean_difference(2)==5);
cfgTiny=cfg; cfgTiny.scb.nBootstrap=3;
tiny=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfgTiny);
assert(isinf(tiny.simultaneousBands.criticalValue));
assert(~any(tiny.simultaneousBands.positive,'all'));
none=temporal_omnibus_signflip(zeros(32,2,7),(1:32)',[],[],(0:6)',cfg);
assert(~any(none.simultaneousBands.eligible,'all') && isinf(none.simultaneousBands.criticalValue));

% Inferential eligibility is conditional; band calculation needs no null
% symmetry (the accompanying omnibus engine retains its symmetry gate).
primary=r; primary.config.runMode='inference';
primary.config.resampling.independentSubjectsVerified=true;
primary.config.resampling.verificationNote='Synthetic independent subjects';
primary.config.resampling.nullSymmetryAssumed=false;
inference=simultaneous_confidence_bands(D,primary);
assert(inference.inferenceValid);
primary.config.resampling.independentSubjectsVerified=false;
blocked=false;
try
    simultaneous_confidence_bands(D,primary);
catch ME
    blocked=strcmp(ME.identifier,'simultaneous_confidence_bands:IndependenceUnresolved');
end
assert(blocked);

% Default whole-subject bootstrap-t: exact reference with globally shared
% subject draws and recomputed SE, including partially overlapping masks.
cfgT=cfg; cfgT.scb=rmfield(cfgT.scb,'method');
rt=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfgT);
bt=rt.simultaneousBands;
assert(strcmp(bt.method,'subject-bootstrap-t'));
stream=RandStream('mt19937ar','Seed',bt.config.seed);
draws=randi(stream,32,32,bt.config.nBootstrap);
reference=zeros(bt.config.nBootstrap,1);
for k=1:bt.config.nBootstrap
    for j=1:3
        selected=draws(:,k);
        selected=selected(rt.subjectMasks(selected,j));
        x=reshape(D(selected,j,:),numel(selected),7)-rt.mean_D(j,:);
        z=mean(x,1)./(std(x,0,1)/sqrt(numel(selected)));
        reference(k)=max(reference(k),max(abs(z)));
    end
end
assert(max(abs(bt.bootstrapMaxima-reference))<1e-11);
cfgT.scb.batchSize=73; cfgT.scb.testBatchSize=1;
rt2=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',cfgT);
assert(max(abs(bt.bootstrapMaxima-rt2.simultaneousBands.bootstrapMaxima))<1e-11);
small=temporal_omnibus_signflip(reshape([1;3],2,1,1),(1:2)',[],[],0,cfgT);
assert(small.simultaneousBands.degenerateBootstrapReplicates>0);
assert(isinf(small.simultaneousBands.criticalValue));
assert(~any(small.simultaneousBands.positive,'all'));

% Adding more cells cannot reduce the jointly calibrated maximum.
one=temporal_omnibus_signflip(D(:,1,:),(1:32)',[],[],(0:6)',cfgT);
assert(all(bt.bootstrapMaxima>=one.simultaneousBands.bootstrapMaxima-1e-11));

% Generic persistence and band-only overwrite collision before any writes.
scratch=tempname; mkdir(scratch); cleanup=onCleanup(@() rmdir(scratch,'s'));
saveCfg=cfg; saveCfg.saveOutputs=true; saveCfg.outputDir=scratch;
saved=temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',saveCfg);
csv=readtable(saved.outputFiles.simultaneousBandsCSV);
assert(height(csv)==21 && max(abs(csv.lower_simultaneous-b.table.lower_simultaneous))<1e-10);
payload=load(saved.outputFiles.mat,'results');
assert(isequal(payload.results.simultaneousBands.lower,b.lower));
collision=fullfile(scratch,'collision'); mkdir(collision);
writetable(b.table,fullfile(collision,'simultaneous_confidence_bands.csv'));
saveCfg.outputDir=collision; blocked=false;
try
    temporal_omnibus_signflip(D,(1:32)',[],[],(0:6)',saveCfg);
catch ME
    blocked=contains(ME.message,'overwrite');
end
assert(blocked && ~isfile(fullfile(collision,'temporal_omnibus_results.mat')));
% Four-category names and explicit sampled-lag plotting helper smoke check.
plot_simultaneous_confidence_bands(rt,1:3,fullfile(scratch,'bands.png'));
assert(isfile(fullfile(scratch,'bands.png')));
fprintf('All simultaneous confidence band tests passed.\n');
end

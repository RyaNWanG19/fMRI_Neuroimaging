function plot_simultaneous_confidence_bands(results, testIndices, outFile)
% Plot selected tests from bands calibrated over the FULL planned family.
% For run_wm_load_omnibus, index = (category_index-1)*489 + atlas_roi_index.
% Display selection does not recalibrate bands; no between-lag coverage claim.
assert(isfield(results,'simultaneousBands'),'Result has no simultaneous bands.');
b=results.simultaneousBands;
validateattributes(testIndices,{'numeric'},{'vector','integer','positive','<=',size(b.mean,1)});
assert(~isempty(testIndices) && numel(testIndices)<=12,'Choose 1-12 tests per figure.');
assert(~isfile(outFile),'Refusing to overwrite an existing figure.');
outDir=fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
ncol=min(3,ceil(sqrt(numel(testIndices)))); nrow=ceil(numel(testIndices)/ncol);
fig=figure('Visible','off','Color','w','Position',[100 100 1200 320*nrow]);
cleanup=onCleanup(@() close(fig));
tl=tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
for r=testIndices(:)'
    ax=nexttile(tl); hold(ax,'on');
    lo=b.lower(r,:); hi=b.upper(r,:); ts=b.timeSec;
    edges=diff([false isfinite(lo)&isfinite(hi) false]);
    starts=find(edges==1); stops=find(edges==-1)-1;
    for k=1:numel(starts)
        ii=starts(k):stops(k);
        fill(ax,[ts(ii);flipud(ts(ii))],[lo(ii)';flipud(hi(ii)')], ...
            [0.78 0.85 0.95],'EdgeColor','none','FaceAlpha',0.6);
    end
    plot(ax,ts,b.mean(r,:),'k.-','LineWidth',1.2,'MarkerSize',7);
    yline(ax,0,'k:');
    pos=b.positive(r,:); neg=b.negative(r,:);
    scatter(ax,ts(pos),b.mean(r,pos),28,[0.8 0.15 0.15],'filled');
    scatter(ax,ts(neg),b.mean(r,neg),28,[0.1 0.35 0.8],'filled');
    name=results.roiNames(r);
    if isfield(results,'testMapping')
        name=results.testMapping.category(r)+" | ROI "+results.testMapping.atlas_roi_index(r);
    end
    title(ax,name,'Interpreter','none');
    xlabel(ax,'Lag since onset (s)'); ylabel(ax,'Mean residual difference'); grid(ax,'on');
end
title(tl,sprintf(['%.1f%% simultaneous bands | %s | %s\n' ...
    'Sampled lags; red/blue dots: positive/negative bands excluding zero; unbounded bands omitted'], ...
    100*(1-b.alpha),results.config.hrfModelName,char(b.inferenceStatus)), ...
    'Interpreter','none','FontSize',10);
exportgraphics(fig,outFile,'Resolution',180);
end

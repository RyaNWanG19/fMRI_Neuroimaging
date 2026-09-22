function report = report_paired_tensor(resultFile, cfg)
%REPORT_PAIRED_TENSOR Local descriptive interpretation, without tensor refitting.
% report_paired_tensor(resultFile, struct('ranks',[1 2]))
% Optional cfg.roiTable: table or CSV with roiIndex, roiName, and optionally
% network. Indices must exactly match the input Results ROI order. No atlas
% ordering is inferred from the number of parcels alone.
if nargin<2, cfg=struct(); end
if ~isfield(cfg,'ranks'), cfg.ranks=[1 2]; end
if ~isfield(cfg,'topN'), cfg.topN=15; end
if ~isfield(cfg,'makePlots'), cfg.makePlots=true; end
if ~isfield(cfg,'outputDir'), cfg.outputDir=fullfile(fileparts(resultFile),'interpretation'); end
validateattributes(cfg.topN,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.ranks,{'numeric'},{'vector','integer','positive'});
assert(numel(unique(cfg.ranks))==numel(cfg.ranks),'Ranks must be unique.');
loaded=load(resultFile,'results'); r=loaded.results;
available=cellfun(@(f) f.rank,r.fits);
assert(all(ismember(cfg.ranks,available)),'Requested rank is absent from this result file.');
nr=size(r.meanContrast,1); t=r.config.timeSec(:);
roi=table((1:nr)',string(r.config.roiNames(:)), ...
    'VariableNames',{'roiIndex','roiName'});
labelSource='Stored labels (may be generic parcel identifiers)';
if isfield(cfg,'roiTable')
    roi=cfg.roiTable;
    if ~istable(roi), roi=readtable(roi,'TextType','string'); end
    assert(all(ismember({'roiIndex','roiName'},roi.Properties.VariableNames)), ...
        'roiTable requires roiIndex and roiName.');
    roi=sortrows(roi,'roiIndex');
    assert(isequal(roi.roiIndex,(1:nr)'), 'roiTable must cover each input ROI index exactly once.');
    roi.roiName=string(roi.roiName);
    labelSource='Explicit cfg.roiTable mapping; caller must verify extraction order';
end
assert(all(~ismissing(roi.roiName)),'Missing ROI labels.');
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
diagnostics=r.rankSummary;
n=height(diagnostics);
diagnostics.convergedStarts=zeros(n,1);
diagnostics.totalStarts=zeros(n,1);
diagnostics.cancellationRatio=nan(n,1);
diagnostics.convergedSplitFits=zeros(n,1);
diagnostics.totalSplitFits=zeros(n,1);
diagnostics.flaggedSplitFits=zeros(n,1);
diagnostics.groupMeanEnergyFit=nan(n,1);
meanEnergy=sum(r.meanContrast(:).^2);
for j=1:n
    idx=find(available==diagnostics.rank(j),1); f=r.fits{idx}; s=r.stability{idx};
    diagnostics.convergedStarts(j)=sum(f.startConverged);
    diagnostics.totalStarts(j)=numel(f.startConverged);
    diagnostics.cancellationRatio(j)=f.cancellationRatio;
    diagnostics.convergedSplitFits(j)=sum(s.converged(:));
    diagnostics.totalSplitFits(j)=numel(s.converged);
    diagnostics.flaggedSplitFits(j)=sum(s.possibleDegeneracy(:));
    e=r.meanContrast-f.A*diag(mean(f.C,1))*f.B';
    if meanEnergy>0, diagnostics.groupMeanEnergyFit(j)=1-sum(e(:).^2)/meanEnergy; end
end
writetable(diagnostics,fullfile(cfg.outputDir,'rank_interpretation.csv'));
writetable(roi,fullfile(cfg.outputDir,'roi_mapping.csv'));
components=table(); topParcels=table(); extremeScores=table(); reconstructions=cell(numel(cfg.ranks),1);
for j=1:numel(cfg.ranks)
    k=cfg.ranks(j); f=r.fits{find(available==k,1)};
    mu=mean(f.C,1); sd=std(f.C,0,1);
    reconstructed=f.A*diag(mu)*f.B'; reconstructions{j}=reconstructed;
    for c=1:k
        [~,peak]=max(abs(f.B(:,c)));
        [~,order]=sort(abs(f.A(:,c)),'descend'); take=order(1:min(cfg.topN,nr));
        ct=table(k,c,mu(c),sd(c),mean(f.C(:,c)>0),t(peak), ...
            f.B(peak,c),f.componentNorm(c), ...
            'VariableNames',{'rank','component','meanScore','sdScore','fractionPositive', ...
            'peakAbsTimeSec','signedTemporalWeightAtPeak','componentNorm'});
        components=[components;ct]; %#ok<AGROW>
        pt=roi(take,:); pt.rank=repmat(k,numel(take),1); pt.component=repmat(c,numel(take),1);
        pt.spatialWeight=f.A(take,c);
        % This slice has original residual units, unlike normalized A/B weights.
        pt.meanContributionAtComponentPeak=f.A(take,c)*mu(c)*f.B(peak,c);
        topParcels=[topParcels;pt]; %#ok<AGROW>
        [~,si]=sort(abs(f.C(:,c)),'descend'); si=si(1:min(cfg.topN,numel(si)));
        scoreEnergy=sum(f.C(:,c).^2);
        share=nan(numel(si),1);
        if scoreEnergy>0, share=f.C(si,c).^2/scoreEnergy; end
        st=table(repmat(k,numel(si),1),repmat(c,numel(si),1), ...
            string(r.retainedSubjectIDs(si)),f.C(si,c),share, ...
            'VariableNames',{'rank','component','subjectID','score','fractionComponentScoreEnergy'});
        extremeScores=[extremeScores;st]; %#ok<AGROW>
    end
    % Correlations are descriptive; near-collinear factors plus cancellation
    % can make individual components harder to interpret than their sum.
    saveCorrelations=struct('spatialCosine',f.A'*f.A,'temporalCosine',f.B'*f.B, ...
        'subjectScoreCorrelation',corrcoef(f.C));
    save(fullfile(cfg.outputDir,sprintf('rank_%02d_factor_similarity.mat',k)),'saveCorrelations');
    if cfg.makePlots
        plot_mean(r.meanContrast,reconstructed,t,cfg.outputDir,k);
        for c=1:k
            plot_component(f,c,roi,t,cfg.topN,cfg.outputDir);
        end
    end
end
writetable(components,fullfile(cfg.outputDir,'component_summary.csv'));
writetable(topParcels,fullfile(cfg.outputDir,'top_parcels.csv'));
writetable(extremeScores,fullfile(cfg.outputDir,'extreme_subject_scores.csv'));
report=struct('sourceFile',resultFile,'config',cfg,'labelSource',labelSource, ...
    'rankDiagnostics',diagnostics,'componentSummary',components, ...
    'topParcels',topParcels,'extremeSubjectScores',extremeScores, ...
    'reconstructions',{reconstructions},'roiMapping',roi);
save(fullfile(cfg.outputDir,'interpretation.mat'),'report');
fid=fopen(fullfile(cfg.outputDir,'README.txt'),'w');
assert(fid>=0,'Cannot open report README.'); closeFile=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Source: %s\nROI labels: %s\n\n',resultFile,labelSource);
fprintf(fid,['Descriptive in-sample report; no significance tests or final rank selection.\n' ...
    'Group-mean fit = 1 - squared reconstruction error / observed mean energy.\n' ...
    'A/B weights are normalized; group-mean component = A(:,k)*B(:,k)''*mean(C(:,k)).\n' ...
    'Individual component energies do not add because components can overlap/cancel.\n' ...
    'Top parcels are ranked by absolute spatial weight, not statistical significance.\n' ...
    'Extreme subject scores are review candidates, not grounds for exclusion.\n' ...
    'Score-energy fractions are not refit-based influence measures.\n' ...
    'Split-half components are not matched to full-data component indices.\n' ...
    'Inspect convergence and cancellation before interpreting any selected rank.\n' ...
    'Generic Parcel labels require a verified ROI mapping before anatomical interpretation.\n']);
disp(diagnostics); fprintf('Interpretation report: %s\n',cfg.outputDir);
end

function plot_mean(observed,reconstructed,t,out,k)
f=figure('Visible','off','Color','w','Position',[100 100 1250 650]);
cleanup=onCleanup(@() close(f)); %#ok<NASGU>
tiledlayout(2,3,'TileSpacing','compact');
limit=max(abs([observed(:);reconstructed(:);observed(:)-reconstructed(:)]));
if limit==0, limit=1; end
arrays={observed,reconstructed,observed-reconstructed};
titles={'Observed group mean','CP group-mean reconstruction','Unexplained group mean'};
for i=1:3
    nexttile; imagesc(t,1:size(observed,1),arrays{i}); axis xy; clim([-limit limit]);
    colorbar; title(titles{i}); xlabel('Time (s)'); ylabel('Parcel index');
end
% Diverging blue-white-red map, with common symmetric scale across panels.
u=linspace(0,1,128)'; colormap(f,[u u ones(128,1);ones(128,1) flipud(u) flipud(u)]);
nexttile; scatter(observed(:),reconstructed(:),3,'.'); hold on;
lim=max(abs([observed(:);reconstructed(:)])); if lim==0, lim=1; end
plot([-lim lim],[-lim lim],'k:'); xlabel('Observed mean'); ylabel('Reconstructed mean'); grid on;
nexttile; plot(t,[sqrt(mean(observed.^2,1));sqrt(mean(reconstructed.^2,1))]');
xlabel('Time (s)'); ylabel('RMS across parcels'); legend('Observed','Reconstructed','Location','best');
nexttile; plot(t,sqrt(mean((observed-reconstructed).^2,1)));
xlabel('Time (s)'); ylabel('Unexplained RMS across parcels');
sgtitle(sprintf('Rank %d: descriptive group-mean reconstruction',k));
exportgraphics(f,fullfile(out,sprintf('rank_%02d_group_mean.png',k)),'Resolution',150);
end

function plot_component(f,c,roi,t,topN,out)
[~,order]=sort(abs(f.A(:,c)),'descend'); take=order(1:min(topN,height(roi)));
fig=figure('Visible','off','Color','w','Position',[100 100 1300 650]);
cleanup=onCleanup(@() close(fig)); %#ok<NASGU>
tiledlayout(2,2,'TileSpacing','compact');
nexttile([2 1]); barh(f.A(take,c)); set(gca,'YTick',1:numel(take), ...
    'YTickLabel',roi.roiName(take),'YDir','reverse','TickLabelInterpreter','none');
xlabel('Signed spatial weight'); title('Largest absolute parcel weights'); xline(0,':');
nexttile; plot(t,f.B(:,c),'LineWidth',1.5); yline(0,':'); xlabel('Time (s)'); ylabel('Normalized temporal weight');
nexttile; histogram(f.C(:,c)); xline(0,':'); xline(mean(f.C(:,c)),'r-','Mean');
xlabel('Signed subject score'); ylabel('Count');
sgtitle(sprintf('Rank %d, component %d: exploratory weights',f.rank,c));
exportgraphics(fig,fullfile(out,sprintf('rank_%02d_component_%02d.png',f.rank,c)),'Resolution',150);
end

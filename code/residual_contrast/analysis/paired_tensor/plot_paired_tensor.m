function plot_paired_tensor(results)
%PLOT_PAIRED_TENSOR Save rank diagnostics and signed component summaries.
cfg=results.config; s=results.rankSummary;
f=figure('Visible','off','Color','w','Position',[100 100 1100 400]);
closer=onCleanup(@() close(f));
tiledlayout(1,2);
nexttile; plot(s.rank,[s.cpEnergyFit,s.svdEnergyFit],'-o');
xlabel('Rank'); ylabel('1 - SSE / uncentered energy'); legend('CP','Uncentered SVD','Location','best'); grid on;
nexttile; plot(s.rank,s.medianSplitHalfCongruence,'-o'); ylim([0 1]);
xlabel('Rank'); ylabel('Median matched split-half congruence'); grid on;
sgtitle([cfg.hrfModelName ': ' cfg.contrastName],'Interpreter','none');
exportgraphics(f,fullfile(cfg.outputDir,'rank_diagnostics.png'),'Resolution',150);
clear closer
for j=1:numel(results.fits)
    fit=results.fits{j};
    f=figure('Visible','off','Color','w','Position',[100 100 1250 250*fit.rank]);
    closer=onCleanup(@() close(f));
    tiledlayout(fit.rank,3,'TileSpacing','compact');
    for k=1:fit.rank
        nexttile; plot(fit.A(:,k)); yline(0,':'); xlabel('Parcel index'); ylabel('Spatial weight'); title(sprintf('Component %d',k));
        nexttile; plot(cfg.timeSec,fit.B(:,k)); yline(0,':'); xlabel('Time from onset (s)'); ylabel('Temporal weight');
        nexttile; histogram(fit.C(:,k)); xline(0,':'); xlabel('Signed subject score'); ylabel('Count');
    end
    sgtitle(sprintf('%s: rank %d; descriptive components',cfg.hrfModelName,fit.rank),'Interpreter','none');
    exportgraphics(f,fullfile(cfg.outputDir,sprintf('rank_%02d_components.png',fit.rank)),'Resolution',150);
    clear closer
end
end

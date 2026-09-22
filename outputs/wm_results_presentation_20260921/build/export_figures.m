build=fileparts(mfilename('fullpath'));
out=fullfile(fileparts(build),'output','figures');
if ~isfolder(out), mkdir(out); end
d=jsondecode(fileread(fullfile(build,'results.json')));
colors=[22 117 169;212 119 38;106 74 162]/255;
names=string({d.models.name});
f=figure('Visible','off','Color','w','Position',[100 100 1000 650]);
b=bar([d.models.significant]); b.FaceColor='flat'; b.CData=colors;
set(gca,'XTick',1:3,'XTickLabel',names,'FontSize',17,'Box','off');
ylim([0 489]); ylabel('ROIs passing Holm correction');
text(1:3,[d.models.significant]+12,string([d.models.significant]), ...
    'HorizontalAlignment','center','FontSize',21,'FontWeight','bold');
title('Reliable WM residual contrasts');
subtitle('Body load difference minus face load difference; Holm p <= 0.05 per model');
exportgraphics(f,fullfile(out,'discovery_counts.png'),'Resolution',250);close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1350 650]);
bar(d.patternCounts,'FaceColor',colors(1,:));
set(gca,'XTick',1:7,'XTickLabel',{'All three','cHRF + deriv','cHRF + sHRF', ...
    'deriv + sHRF','cHRF only','deriv only','sHRF only'},'FontSize',15,'Box','off');
ylim([0 220]);ylabel('Number of ROIs');
text(1:7,d.patternCounts+7,string(d.patternCounts),'HorizontalAlignment','center','FontSize',20);
title('Overlap of significant ROIs across HRF models');
subtitle('Mutually exclusive groups; 104 shared, 343 in any model, 146 in none');
exportgraphics(f,fullfile(out,'model_overlap.png'),'Resolution',250);close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1500 590]);
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for j=1:3
    ax=nexttile;hold(ax,'on');
    for k=1:3
        plot(d.timeSec,d.examples(j).mean(k,:),'Color',colors(k,:),'LineWidth',2,'DisplayName',names(k));
    end
    yline(0,':k','HandleVisibility','off');
    title(sprintf('%s (ROI %d)',d.examples(j).name,d.examples(j).id),'Interpreter','none');
    xlabel('Time (s)'); if j==1, ylabel('Mean residual contrast');end
    xlim([0 28.8]);xticks(0:7.2:28.8);set(ax,'FontSize',16,'Box','off');
    if j==2, legend('Location','southoutside','Orientation','horizontal'); end
end
title(tl,'Illustrative population-average residual contrast curves','FontSize',23);
subtitle(tl,'Selected after testing; separate y scales; no individual-lag significance implied','FontSize',16);
exportgraphics(f,fullfile(out,'example_trajectories.png'),'Resolution',250);close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1000 650]);
bar([[d.models.significant]' [d.models.bootstrap]']);
set(gca,'XTick',1:3,'XTickLabel',names,'FontSize',17,'Box','off');
ylabel('ROIs passing Holm correction');ylim([0 400]);
legend({'Sign flips','Centered bootstrap'},'Location','northoutside','Orientation','horizontal');
title('Sensitivity to the resampling method');
subtitle('Approximate bootstrap comparison; 99,999 replicates per method');
exportgraphics(f,fullfile(out,'bootstrap_comparison.png'),'Resolution',250);close(f);

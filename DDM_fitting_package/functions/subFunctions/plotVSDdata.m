function plotVSDdata(f, color)


loglog(f.stimDurs, f.dataSens,'ko', 'Markerfacecolor', 'k', 'Markersize', 1.5)
hold on
errorbar(f.stimDurs, f.dataSens, f.dataSensSE, 'linestyle','none', 'Color', 'k', 'CapSize', 0) ;
loglog(f.stimDurs, f.predSens, 'Color',color,'LineWidth', 1)
xlim([.05 1])
set(gca, 'xtick', [.05 .1 .2 .4 .8]) ;
set(gca, 'ytick', [5 10 20 40]) ;
% ylim([5 25])
formatFig(gcf, [1.4 1.4], 'jnsci');
xlabel('Stimulus duration')
ylabel('Sensitivity') 

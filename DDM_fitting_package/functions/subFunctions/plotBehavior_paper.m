function plotBehavior_paper(d, f, plotFits, color) 

% clf
C = [d.coherence] ;
RTs = nanmean([d.RTs]) ;
choice = nansum([d.choice]) ./ [d.nTrials] ;

p_sem = sqrt( (choice .* (1-choice)) ./ [d.nTrials]) ;
rt_sem = nanstd([d.RTs]) ./ sqrt([d.nTrials]) ;

% figure


% xBreak1 = [-.95 -.55] ;
% xBreak2 = [.55 .95] ;

xl = .55 ;
xl2 = .95 ;

newXl = .65 ;
dx = xl2 - newXl ;

l = f.c4p>-xl & f.c4p<xl ;
c2 = f.c4p(l) ;
pc2 = f.predChoice(l) ;
rt2 = f.predRTmean(l) ;

l2 = f.c4p<-xl2 | f.c4p>xl2 ;
c3 = f.c4p(l2) ;
pc3 = f.predChoice(l2) ;
rt3 = f.predRTmean(l2) ;
c3(c3<0) = c3(c3<0) + dx ;
c3(c3>0) = c3(c3>0) - dx ;
c3 = [c3(c3<0) 0 c3(c3>0)] ;
pc3 = [pc3(1:end/2) nan pc3(end/2+1:end)] ;
rt3 = [rt3(1:end/2) nan rt3(end/2+1:end)] ;

C2 = C ;
C2([1 end]) = [-.7 .7] ;


% l1 = f.c4p > xBreak1(1) & f.c4p < xBreak1(2) ;
% l2 = f.c4p > xBreak2(1) & f.c4p < xBreak2(2) ;
% l = logical(l1 + l2) ;
% 
% f.predChoice(l) = nan ;
% f.predRTmean(l) = nan ;





subplot(2,1,2)
hold on
if plotFits
plot(c2, pc2, 'Color', color, 'linewidth', 1)
plot(c3, pc3, 'Color', color, 'linewidth', 1)
end
hold on
plot(C2, choice, 'ko', 'Markerfacecolor', 'k', 'MarkerSize', 1.5)
% BreakPlot(choice, C, -.9, -.6, .2)
hold on
errorbar(C2, choice, p_sem, 'linestyle', 'none', 'Color', 'k', 'linewidth', .5, 'CapSize', 0) 


hold on
plot([-.6 -.6], [0 .05], 'k-')
plot([.6 .6], [0 .05], 'k-')
% axis square
ylim([-.01 1])
xlim([-.512 .512])
xlabel('Motion strength (% coh)')
ylabel('Proportion rightward choice')
set(gca, 'xtick',[-.7 -.5 0 .5 .7]) ;
set(gca, 'xticklabel', [-100 -50 0 50 100]) ;
set(gca, 'ytick',[0 .5 1]) ;

subplot(2,1,1)
hold on
if plotFits
plot(c2, rt2, 'Color', color, 'linewidth', 1)
plot(c3, rt3, 'Color', color, 'linewidth', 1)
end
hold on
plot(C2, RTs, 'ko', 'Markerfacecolor', 'k', 'Markersize', 1.5)
% BreakPlot(choice, C, -.9, -.6, .2)
hold on
errorbar(C2, RTs, rt_sem, 'linestyle', 'none', 'Color', 'k', 'linewidth',.5, 'CapSize', 0) 


hold on
% axis square
% ylim([.2 2.5])
xlim([-.512 .512])
% xlabel('Motion strength (% coh)')
ylabel('Reaction time (s)')
set(gca, 'xtick',[-.7 -.5 0 .5 .7]) ;
set(gca, 'xticklabel', []) ;
yl= get(gca, 'ylim') ;
y = yl(1) ;
plot([-.6 -.6], [y y+.05], 'k-')
plot([.6 .6], [y y+.05], 'k-')
ylim([y yl(2)])
set(gca, 'ytick', .2:.2:yl(2)+.3)
% dMat = dStruct2dMatrix(d) ;



% Get percent correct and mean rt for each coherence

% uCoh = unique(dMat(:,4)) ;
% 
% for i = 1:length(uCoh)
%     
%     l = dMat(:,4) == uCoh(i) ;
%     n(i) = sum(l) ;
%     
%     pCorrect(i) = sum(dMat(l,5)) / n(i) ;
%     meanRT_pc(i) = mean(dMat(l,3)) ;
%     stdRT_pc(i) = std(dMat(l,3)) ;
%     
% end

% % uCoh(uCoh==0.999) = 1 ;
% uCoh = uCoh*100 ;
% uCoh(1) = uCoh(2) * .3125 ;
% % f.c4p_pCorrect(1) = eps ;
% f.predChoice_pCorrect(f.c4p_pCorrect*100 <uCoh(2)*.8) = nan ;
% f.predRTmean_pCorrect(f.c4p_pCorrect*100 <uCoh(2)*.8 & f.c4p_pCorrect*100 > uCoh(1)*1.5) = nan ;
% 
% pCorrect(uCoh==min(uCoh)) = .5 ;
% 
% pC_sem = sqrt( (pCorrect.*(1-pCorrect)) ./ n) ;
% pC_sem(uCoh==min(uCoh)) = 0 ;
% rtPC_sem = stdRT_pc ./ sqrt(n) ;
% 
% subplot(2,1,2)
% hold on
% plot(f.c4p_pCorrect*100, f.predChoice_pCorrect, 'Color', color, 'linewidth',2)
% plot(uCoh, pCorrect, 'ko','markerfacecolor','k')
% errorbar(uCoh, pCorrect, pC_sem, 'linestyle','none', 'Color', 'k')
%     
% set(gca,'xscale','log', 'xtick',uCoh,'TickDir','out','FontSize',10)
% xlim([uCoh(1)*.8 100])
% ylim([.45 1])
% ylabel('Percent correct')
% xlabel('Stimulus strength')
% 
% 
% subplot(2,1,1)
% hold on
% plot(f.c4p_pCorrect*100, f.predRTmean_pCorrect, 'Color', color,'linewidth',2)
% plot(uCoh, meanRT_pc, 'ko','markerfacecolor', 'k')
% errorbar(uCoh, meanRT_pc, rtPC_sem, 'linestyle','none', 'Color', 'k')
%     
% set(gca,'xscale','log', 'xtick',uCoh,'TickDir','out','FontSize',10)
% xlim([uCoh(1)*.8 100])
% ylabel('Reaction time (s)')
% 
formatFig(gcf, [1.4 2.5], 'jnsci');

% subplot(2,2,3)
% hold on
% plot(f.predRT_t, f.bU)
% plot(f.predRT_t, f.bL)
% xlim([0 7])
% axis square



teste = results_distribution.RDS.duration.ripple;
% 
% histogram(teste{1,1},'Normalization','probability')
% hold on
% histogram(teste{2,1},'Normalization','probability')


% Total events
subplot(1,3,[1 2])
% Calculate medians
d1 = teste{1,1};
d2 = teste{2,1};
d3 = teste{3,1};

m1 = median(d1);
m2 = median(d2);
m3 = median(d3);

x = (0.04:0.01:0.2);
[N,edges] = histcounts(d1,x, 'Normalization', 'probability');
semilogx(edges(1:end-1),N,'Color',[.6 .6 .6],'LineWidth',3)
hold on
[N,edges] = histcounts(d2,x, 'Normalization', 'probability');
semilogx(edges(1:end-1),N,'Color',[255 204 102]/255,'LineWidth',3)


% Plot medians
xline(m1,'Color',[.6 .6 .6],'LineStyle','--','LineWidth',2)
xline(m2,'Color',[255 204 102]/255,'LineStyle','--','LineWidth',2)
ylabel('Fraction of SWR')
xlabel('Duration (s)')
legend('BL','D1')
title('SWR Duration BL vs D1')
legend boxoff
box off
set(gca,'fontsize',18)
set(gca,'Tickdir','out')
set(gca,'Linewidth',1.5)
set(gcf,'color',[1 1 1]);


hold off


% subplot(1,3,2)
% % Only events higher than 100 ms
% d1 = teste{1,1}(teste{1,1} >= 0.1);
% d2 = teste{2,1}(teste{2,1} >= 0.1);
% 
% % Calculate medians
% m1 = median(d1);
% m2 = median(d2);
% 
% x = (0.1:0.01:0.3);
% [N,edges] = histcounts(d1,x, 'Normalization', 'probability');
% semilogx(edges(2:end),N,'b')
% hold on
% [N,edges] = histcounts(d2,x, 'Normalization', 'probability');
% semilogx(edges(2:end),N,'r')
% 
% % Plot medians
% xline(m1,'Color','b')
% xline(m2,'Color','r')
% hold off


subplot(1,3,3)
% Join days
d = [d1; d2];

g1 = repmat({'BL'},length(d1),1);
g2 = repmat({'Second'},length(d2),1);
g = [g1; g2];

boxplot(d,g)
title('SWR events >= 100 ms')
xlabel('Day')
ylabel('Seconds')

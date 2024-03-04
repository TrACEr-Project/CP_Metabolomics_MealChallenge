function plot_allscores(allscores)


boxplot(allscores); ylabel('FMS'); xlabel('Number of Components')
set(gca,'XTick', 1:1:4, 'XTickLabel',{'$R=1$', '$R=2$','$R=3$', '$R=4$'},'TickLabelInterpreter','latex')
index= round(450*0.95);
for r=1:4
    aaa = sort(allscores(:,r),'descend');
    hold on;
    plot(r-0.25:0.005:r+0.25, ones(101,1)*aaa(index),'g-');
end

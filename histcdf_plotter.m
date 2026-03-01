function histcdf_plotter(all,tagged,BinEdges)
figure
hold on
h1 = histogram(all,'BinEdges',BinEdges,'Normalization','count','FaceColor', 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
histogram(tagged,'BinEdges',h1.BinEdges,'Normalization','count', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
ylabel('Count')
yyaxis('right')
h2 = histogram(all,'BinEdges',BinEdges,'Normalization','cdf', 'DisplayStyle','stairs', 'EdgeColor', 'k');
histogram(tagged,'BinEdges',h2.BinEdges,'Normalization','cdf', 'DisplayStyle','stairs', 'EdgeColor', 'b');
ylim([0,1])
 setmyplot_balazs
end
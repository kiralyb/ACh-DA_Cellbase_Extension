function combHhist(H_d,H_a)
figure
hold on
histogram(H_d,'BinEdges',0:0.01:1, 'FaceColor', 'm', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
histogram(H_a,'BinEdges',0:0.01:1, 'FaceColor', 'c', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
setmyplot_balazs
xlabel('H-index')
ylabel('Count')
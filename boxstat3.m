function boxstat3(a,b,c,labels,ispaired)
if nargin<5
    ispaired = 0;
end
%figure
hold on
boxplot(a,'position',0)
boxplot(b,'position',1)
boxplot(c,'position',2)
if ispaired
    p(1) = signrank(a,b)
    p(2) = signrank(a,c)
    p(3) = signrank(b,c)
else
    p(1) = ranksum(a,b)
    p(2) = ranksum(a,c)
    p(3) = ranksum(b,c)
end
sigstar({[0,1], [0,2], [1,2]},p)
xticks([0,1,2])
xticklabels(labels)
setmyplot_balazs
n = length(a);
if ispaired
for i = 1:n
    plot([0, 1, 2], [a(i), b(i), c(i)], 'Color', [0.6 0.6 0.6])
end
end

allLines = findall(gca, 'Type', 'Line');

for i = 1:length(allLines)
    allLines(i).LineStyle = '-'; % change from dashed '--' to solid '-'
end
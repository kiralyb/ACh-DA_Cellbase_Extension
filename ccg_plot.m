function ccg_plot(CCG,mlim)
if nargin < 2
    mlim = 1;
end
minx = (size(CCG,2)-1)/2 + 1;
lags = -(minx - 1) : (minx - 1);
[~,sorter] = sort(sum(CCG(:,minx : minx + 100),2));% - sum(CCG(:,minx - 100 : minx),2));
imagesc(lags,1:size(CCG,1),CCG(sorter,:))
%mlim = 0.7;
caxis([-mlim,mlim])
xlabel('lag (ms)')
%xlim([-500,500])
yyaxis right
plot(lags(1:20:end),nanmean(CCG(:,1:20:end)),'r','LineWidth',2)
hold on
line([0,0],[-mlim,mlim],'Color','w','LineWidth',1);
ylim([-mlim,mlim])
%line([0,0],[1,length(pairstypes{i})],'Color','w','LineWidth',1)
setmyplot_balazs
[~,minI]=min(CCG(sorter,:),[],2);
signrank(minI-lags(end)-1)
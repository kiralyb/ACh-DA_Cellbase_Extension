%ACCG = function ccg_plotter(file,ie_pairs)
%load(file)
CCR_norm2 = bsxfun(@rdivide,CCR,sum(CCR,2));
CCR_Poi_norm = bsxfun(@rdivide,CCR_Poi,sum(CCR_Poi,2));
CCR_norm = (CCR_norm2 - CCR_Poi_norm);
L = (length(CCR)-1)/2;
lags = -(length(CCR)-1)/2:(length(CCR)-1)/2;
%CCR_c = CCR(:,lags+1+(size(CCR,2)-1)/2);
%CCR_norm = bsxfun(@rdivide,CCR_c,sum(CCR_c,2));
%CCR_smooth = nan(size(CCR_norm));
clear CCR_smooth
for ind = 1:size(CCR_norm,1)
    %if sum(CCR_c(ind,:),2) > 0
    CCR_smooth(ind,:) = smoothdata(CCR_norm(ind,:),2,'gaussian',(length(CCR)/sqrt(SegmentLength(ind))));
    %end
end
figure
for i = 1:2
subplot(1,2,i)
CCR_smoothi = CCR_smooth(ie_pairs == i-1,:)*500;
[~,sorter] = sort(sum(CCR_smoothi(:,L:L+L/5),2));%-sum(CCR_smooth(pairstypes{i},900:1000)'));
imagesc(lags,1:size(CCR_smoothi,1),CCR_smoothi(sorter,:))
mlim = 0.0005;%max(CCR_smoothi,[],'all');
%caxis([-mlim,mlim]);
caxis([-mlim,mlim])
xlabel('lag (ms)')
xlim([-500,500])
yyaxis right
plot(lags,nanmean(CCR_smoothi),'r','LineWidth',2)
hold on
line([0,0],[-mlim,mlim],'Color','w','LineWidth',1);
ylim([-mlim,mlim])
%line([0,0],[1,length(pairstypes{i})],'Color','w','LineWidth',1)
setmyplot_balazs
[~,minI]=min(CCR_smoothi(sorter,:),[],2);
signrank(minI-lags(end))
end
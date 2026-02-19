function AVG_CCG_plot(path,ie_pairs)
figure
Colors = {'k','g','r'};
for col = 1:3
    switch col
        case 1
            [CCG,CCG_r] = CCG_normalizer([path,'Cue\CCG_matrices.mat']);
        case 2
            [CCG,CCG_r] = CCG_normalizer([path,'Hit\CCG_matrices.mat']);
        case 3
            [CCG,CCG_r] = CCG_normalizer([path,'FA\CCG_matrices.mat']);
    end
for i = 1:2

subplot(2,2,i)
hold on
lags = -(length(CCG)-1)/2:(length(CCG)-1)/2;
errorshade(lags,nanmean(CCG(ie_pairs == i-1,:)),nanstd(CCG(ie_pairs == i-1,:))/sqrt(sum(ie_pairs == i-1)),'LineColor',Colors{col},'ShadeColor',Colors{col})
line([0,0],[-1,2])
line([lags(1),lags(end)],[0,0])
setmyplot_balazs

subplot(2,2,i+2)
hold on
[~,mind] = max(abs(nanmean(CCG(ie_pairs == i-1,(length(CCG)+2)/2 - 20:end))),[],2);
mind = mind+(length(CCG)+1)/2 - 20;
boxplot(mean(CCG_r(ie_pairs == i-1,mind-20:mind+20),2),'position',col)

line([0,3],[0,0])
sigstar([col,col],signrank(mean(CCG_r(ie_pairs == i-1,mind-20:mind+20),2)))
ylim([-1.2,1.7])
setmyplot_balazs
end
end
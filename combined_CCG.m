function [p,meanlag] = combined_CCG
%choosecb('')
load([getpref('cellbase','datapath'),'\CCG.mat'])
load([getpref('cellbase','datapath'),'\CCG_DREADD.mat'])
AVCCG = [Avg_CCG;Avg_CCG_DREADD];
[MX,MI]= max(AVCCG,[],2)
figure
imagesc(lags/sr,1:size(AVCCG,1),AVCCG./MX)
xlim([-6,6])
caxis([-0.8,0.8]);
yyaxis right
plot(lags/sr,mean(AVCCG./MX),'r')
xlim([-6,6])
ylim([-1,1])
setmyplot_balazs
meanlag = (mean(lags(MI)))/sr;
p = signrank(lags(MI));
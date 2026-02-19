function [BurstIndex,Refractory,ThetaIndex] = ACG_plotter(i_neuron)

load('L:\_CellBases\ACh_DA_Cellbase\ACG\hdb_gaba\stim_filter\ACG_matrices_.mat')
inx = [];
for i = 1: length(cellids)
   if ismember(findcellpos(cellids(i)),i_neuron)
       inx = [inx, i];
   end
end
swin = 10;
CCR_norm = bsxfun(@rdivide,CCR(inx,:),sum(CCR(inx,:),2));
CCR_smooth = movmean(CCR_norm,swin,2);
BurstIndex = BurstIndex(inx);
Refractory = Refractory(inx);
ThetaIndex = ThetaIndex(inx);
figure
subplot(2,2,1)
errorshade(-1000:0.5:1000,mean(CCR_smooth(:,:)),std(CCR_smooth(:,:))/sqrt(length(inx)))
hold on
%errorshade(-1000:0.5:1000,mean(CCR_smooth(i_ach_burst,:)),std(CCR_smooth(i_ach_burst,:))/sqrt(length(i_ach_burst)))
line([0,0],[0.0000,0.002])
ylim([0,0.001])
setmyplot_balazs;
subplot(2,2,3)
[~,sorter] = sort(BurstIndex);
imagesc(-1000:0.5:1000,1:length(cellids(inx)),CCR_smooth(sorter,:))
caxis([0,0.0005]);
xlabel('lag (ms)')
setmyplot_balazs;
subplot(2,2,2)
scatter(Refractory,BurstIndex)
xlim([0,101])
ylim([-1,1])
xlabel('Refractory')
ylabel('BurstIndex')
setmyplot_balazs;
subplot(2,2,4)
scatter(Refractory,ThetaIndex)
xlabel('Refractory')
ylabel('ThetaIndex')
xlim([0,100])
ylim([-0.2,0.7])
setmyplot_balazs;

%i_ach_reg = i_neuron(inx(Refractory' >= 25));
%i_ach_burst = i_neuron(inx(Refractory' < 25));
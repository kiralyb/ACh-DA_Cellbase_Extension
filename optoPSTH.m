function p = optoPSTH(IDS,stimtype,color,x)
j = 1;
TrigEvent = 'BurstOn';
win = [-3 5];
dt = 0.02;
sigma = 0.2;

for iCell = 1:length(IDS)
        [~,spsth_BO] = ultimate_psth(IDS{iCell},stimtype,TrigEvent,win,'dt',dt,'sigma',sigma);
        Zspsth(j,:) = (spsth_BO - mean(spsth_BO((sigma*3/dt)+1:-win(1)/dt))) ./ std(spsth_BO((sigma*3/dt)+1:-win(1)/dt)); 
        j = j + 1;
end
subplot(1,2,1)
hold on
errorshade(win(1):dt:win(2),nanmean(Zspsth),nanstd(Zspsth)./sqrt(length(IDS)),'LineColor',color,'ShadeColor',color)
xlabel('Time from cholinergic stimulation start (s)')
ylabel('Average normalized firing rate')
setmyplot_balazs
xlim([-2,4])
subplot(1,2,2)
violin(mean(Zspsth(:,-win(1)/dt:((-win(1)+2)/dt)),2),'x',x)
hold on
line([-0.5,0.5]+x',[median(mean(Zspsth(:,-win(1)/dt:((-win(1)+2)/dt)),2)),median(mean(Zspsth(:,-win(1)/dt:((-win(1)+2)/dt)),2))])
p = signrank(mean(Zspsth(:,-win(1)/dt:((-win(1)+2)/dt)),2));
sigstar([x,x],p)
xlim([-0.5,2.5])


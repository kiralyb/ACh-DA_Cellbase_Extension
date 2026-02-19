function avg_viewlick(RATs)

figure
hold on
Colors = [0,1,0,0.3;1,0,0,0.3;0,0.5,0,0.3;0.5,0,0,0.3]
twind = [-2,3];
tbins = (twind(2) - twind(1)) * 1000 + 1;

LIST = listtag('session');
for ratindex = 1:length(RATs)
    RATSES = LIST(ismember(LIST(:,1), RATs(ratindex)), :);
    RatPSTH = nan(length(RATSES),4,tbins);
    for sesindex = 1:length(RATSES)
        [~,spsth,~,TAGS] = ultimate_psth(RATSES(sesindex,:),'lick','StimulusOn',twind,'parts','#Type5');
        RatPSTH(sesindex,cellfun(@(x) str2double(x(end)), TAGS),:)=spsth; 
    end
    MPSTHs(ratindex,:,:) = nanmean(RatPSTH);           
end

figure
for part = 1:4
        errorshade(linspace(-2,3,501),squeeze(mean(MPSTHs(:,part,1:10:end)))',squeeze(std(MPSTHs(:,part,1:10:end)))'/sqrt(length(RATs)),'LineColor',Colors(part,1:3),'ShadeColor',Colors(part,1:3)) 
end
setmyplot_balazs
xlabel('Time from cue (s)')
ylabel('Lick rate (Lick/s)')
legend({'Fixed go','Fixed no-go','New go','New no-go'})
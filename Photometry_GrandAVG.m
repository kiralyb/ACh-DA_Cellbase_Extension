function Photometry_GrandAVG(partition,trigger,signal,isstat)
cellbasepath = getpref('cellbase','datapath');
sr = 12048;
win = [-2,3];
Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
labels = {'Constant go','New go','Constant no-go','New no-go'};
%if %exist()
load([cellbasepath,partition,'_',trigger,'_',signal,'FM.mat'])
%else
%Photometry_AVG_runner
%end
figure
subplot(1,3,1:2)
FM_P = FM(:,:,1:50:end);
parts = 1:size(FM,2);
for k = parts
    errorshade(linspace(-2,3,length(FM_P)),squeeze(nanmean(FM_P(:,k,:))),squeeze(nanstd(FM_P(:,k,:)))/sqrt(size(FM_P,1)),'LineColor',Color{k},'ShadeColor',Color{k})
    hold on
end
xlabel('Time (s)')
ylabel('Average dF/F')
legend(labels(parts))
setmyplot_balazs

subplot(1,3,3)
for parti = parts
    if isequal(trigger,'DeliverFeedback')
        W = 1;
        if signal == 'D' & mod(parti,2)==0
            inv = -1;
        else
            inv = 1;
        end
    else
        inv = 1;
        W = 2;
    end
    [~,resps(parti,:)] = Photoresponse(squeeze(FM(:,parti,-win(1)*sr : -win(1)*sr + sr*W)),5,0.1,inv,sr);
end

if isstat
    resps([2 3],:) = resps([3 2],:);
    boxplot(resps')
    hold on
    plot(resps)
    setmyplot_balazs
    xticks(parts)
    xticklabels(labels(parts))
    
    p(1) = signrank(resps(1,:),resps(2,:))
    if length(parts) < 4
        p(2) = signrank(resps(1,:),resps(3,:))
        sigstar({[1,2],[1,3]},p)
    else
        p(2) = signrank(resps(3,:),resps(4,:))
        sigstar({[1,2],[3,4]},p)
    end
end

function p = C21_COMP(M,MAX,type_animal,parti)
Color = {'g',[0.85,1,0];'r',[1,0.8,0];[0,0.70,0],[0.6,0.70,0];[0.70,0,0],[0.70,0.6,0]};
n = [sum(type_animal==1), sum(type_animal==0)];
downsamplingrate = 100;
figure
    subplot(1,3,1:2)
    hold on
    for at = 1:2
        errorshade(linspace(-2,3,ceil(60241/downsamplingrate)),squeeze(mean(M(type_animal==at-1,parti,1:downsamplingrate:end))),squeeze(std(M(type_animal==at-1,parti,1:downsamplingrate:end)))/sqrt(n(at)),'LineColor',Color{parti,at},'ShadeColor',Color{parti,at},'LineStyle','-')
    end
    ylabel('Average dF/F')
    xlabel(['Time from ','Feedback',' (s)'])
    setmyplot_balazs
    xlim([-0,2])
    subplot(1,3,3)
    scatter([zeros(1,n(1)),ones(1,n(2))],[MAX(type_animal==0,parti);MAX(type_animal==1,parti)],'MarkerEdgeColor','k','MarkerFaceColor','k')
    p = ranksum(MAX(type_animal==0,parti),MAX(type_animal==1,parti));
    sigstar([0,1],p)
    xlim([-0.5,1.5])
    setmyplot_balazs
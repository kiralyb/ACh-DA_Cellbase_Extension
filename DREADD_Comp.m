function p = DREADD_Comp(M,MAX,type_animal,parti)
Color = {'g',[0.85,1,0];'r',[1,0.8,0];[0,0.70,0],[0.6,0.70,0];[0.70,0,0],[0.70,0.6,0]};
n_d = sum(type_animal==1);
downsamplingrate = 100;
figure
    subplot(1,3,1:2)
    hold on
    for dreadd = 1:2
    errorshade(linspace(-2,3,ceil(60241/downsamplingrate)),squeeze(mean(M(type_animal==1,dreadd,parti,1:downsamplingrate:end))),squeeze(std(M(type_animal==1,dreadd,parti,1:downsamplingrate:end)))/sqrt(n_d),'LineColor',Color{parti,dreadd},'ShadeColor',Color{parti,dreadd},'LineStyle','-')
    end
    ylabel('Average dF/F')
    xlabel(['Time from ','Stimulus',' (s)'])
    setmyplot_balazs
    xlim([-0,2])
    
    subplot(1,3,3)
    line([zeros(1,n_d);ones(1,n_d)],[MAX(type_animal==1,:,parti)]','Color','k');
    hold on
    scatter([zeros(1,n_d),ones(1,n_d)],[MAX(type_animal==1,1,parti);MAX(type_animal==1,2,parti)],'MarkerEdgeColor','k','MarkerFaceColor','k')
    p = signrank(MAX(type_animal==1,1,parti),MAX(type_animal==1,2,parti));
    sigstar([0,1],p)
    xlim([-0.5,1.5])
    setmyplot_balazs
end


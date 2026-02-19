function p = firstnew_photometry(M_SO,MAX_SO,type_animal)
Color = {'g',[0.85,1,0];'r',[1,0.8,0];[0,0.70,0],[0.6,0.70,0];[0.70,0,0],[0.70,0.6,0]};
n = [sum(type_animal==0), sum(type_animal==1)];
downsamplingrate = 100;
figure
for k = 1:3
    subplot(3,3,(k-1)*3+1:(k-1)*3+2)
    hold on
    for at = 1:2
        errorshade(linspace(-2,3,ceil(60241/downsamplingrate)),squeeze(mean(M_SO(type_animal==at-1,k,1:downsamplingrate:end))),squeeze(std(M_SO(type_animal==at-1,k,1:downsamplingrate:end)))/sqrt(n(at)),'LineColor',Color{k,at},'ShadeColor',Color{k,at},'LineStyle','-')
    end
    ylabel('Average dF/F')
    xlabel(['Time from ','Stimulus',' (s)'])
    setmyplot_balazs
    xlim([-1,1])
    subplot(3,3,(k-1)*3+3)
    scatter([zeros(1,n(1)),ones(1,n(2))],[MAX_SO(type_animal==0,k);MAX_SO(type_animal==1,k)],'MarkerEdgeColor','k','MarkerFaceColor','k')
    p(k) = ranksum(MAX_SO(type_animal==0,k),MAX_SO(type_animal==1,k));
    sigstar([0,1],p(k))
    xlim([-0.5,1.5])
    setmyplot_balazs
end

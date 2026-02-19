function p = fix_DREADD_comp(MRT,type_animal)

n_d = sum(type_animal==1);
n_c = sum(type_animal==0);
figure
hold on
scatter([zeros(1,n_d),ones(1,n_d)],reshape(MRT(type_animal==1,:),[],1)','MarkerEdgeColor','k','MarkerFaceColor','k')
scatter([zeros(1,n_c),ones(1,n_c)],reshape(MRT(type_animal==0,:),[],1)','MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8])
line([zeros(1,n_d);ones(1,n_d)],[MRT(type_animal==1,:)]','Color','k');
line([zeros(1,n_c);ones(1,n_c)],[MRT(type_animal==0,:)]','Color',[0.8,0.8,0.8]);
xlim([-0.5,1.5])
xticks([0,1])
xticklabels({'Control','C21'})
p = signrank(MRT(type_animal==1,1),MRT(type_animal==1,2));
sigstar([0,1],p)
setmyplot_balazs
legend({'DREADD','Control'})
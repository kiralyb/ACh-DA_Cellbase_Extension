%animalID = {'DHBAVD-04','DHBAVD-07','DHBAVD-10','DHBAVD-12','DHBAVD-17'};
%animalID = {'DHBAVD-01','DHBAVD-02','DHBAVD-13','DHBAVD-15','DHBAVD-16'};
[Hitrate(1),HitrateF(1),jF(1)] = psychometric_7hzbehav_analysis('DHBAVD-04','230317d')
[Hitrate(2),HitrateF(2),FarateF(2)] = psychometric_7hzbehav_analysis('DHBAVD-07','230609d')
[Hitrate(3),HitrateF(3),FarateF(3)] = psychometric_7hzbehav_analysis('DHBAVD-10','230816d')
[Hitrate(4),HitrateF(4),FarateF(4)] = psychometric_7hzbehav_analysis('DHBAVD-12','231002d')
[Hitrate(5),HitrateF(5),FarateF(5)] = psychometric_7hzbehav_analysis('DHBAVD-17','240301d')
[Hitrate(6),HitrateF(6),FarateF(6)] = psychometric_7hzbehav_analysis('DHBAVD-18','240426d')
[Hitrate(7),HitrateF(7),FarateF(7)] = psychometric_7hzbehav_analysis('DHBAVD-19','240418d')

[Hitrate2(1),HitrateF2(1),FarateF2(1)] = psychometric_7hzbehav_analysis('DHBAVD-04','230320d')
[Hitrate2(2),HitrateF2(2),FarateF2(2)] = psychometric_7hzbehav_analysis('DHBAVD-07','230612d')
[Hitrate2(3),HitrateF2(3),FarateF2(3)] = psychometric_7hzbehav_analysis('DHBAVD-10','230817d')
[Hitrate2(4),HitrateF2(4),FarateF2(4)] = psychometric_7hzbehav_analysis('DHBAVD-12','231003d')
[Hitrate2(5),HitrateF2(5),FarateF2(5)] = psychometric_7hzbehav_analysis('DHBAVD-17','240304d')
[Hitrate2(6),HitrateF2(6),FarateF2(6)] = psychometric_7hzbehav_analysis('DHBAVD-18','240429d')
Hitrate2(7) = nan;
HitrateF2(7) = nan;
FarateF2(7) = nan;

[Hitrate3(1),HitrateF3(1),FarateF3(1)] = psychometric_7hzbehav_analysis('DHBAVD-04','230322a')
[Hitrate3(2),HitrateF3(2),FarateF3(2)] = psychometric_7hzbehav_analysis('DHBAVD-07','230613a')
[Hitrate3(3),HitrateF3(3),FarateF3(3)] = psychometric_7hzbehav_analysis('DHBAVD-10','230821a')
[Hitrate3(4),HitrateF3(4),FarateF3(4)] = psychometric_7hzbehav_analysis('DHBAVD-12','231004a')
[Hitrate3(5),HitrateF3(5),FarateF3(5)] = psychometric_7hzbehav_analysis('DHBAVD-17','240305a')
[Hitrate3(6),HitrateF3(6),FarateF3(6)] = psychometric_7hzbehav_analysis('DHBAVD-18','240430a')
Hitrate3(7) = nan;
HitrateF3(7) = nan;
FarateF3(7) = nan;

[Hitrate_c(1),HitrateF_c(1),FarateF_c(1)] = psychometric_7hzbehav_analysis('DHBAVD-01','230119d')
[Hitrate_c(2),HitrateF_c(2),FarateF_c(2)] = psychometric_7hzbehav_analysis('DHBAVD-02','230118d')
[Hitrate_c(3),HitrateF_c(3),FarateF_c(3)] = psychometric_7hzbehav_analysis('DHBAVD-13','230925d')
[Hitrate_c(4),HitrateF_c(4),FarateF_c(4)] = psychometric_7hzbehav_analysis('DHBAVD-15','240318d')
[Hitrate_c(5),HitrateF_c(5),FarateF_c(5)] = psychometric_7hzbehav_analysis('DHBAVD-16','240320d')
[Hitrate_c(6),HitrateF_c(6),FarateF_c(6)] = psychometric_7hzbehav_analysis('DHBAVD-20','240508d')
[Hitrate_c(7),HitrateF_c(7),FarateF_c(7)] = psychometric_7hzbehav_analysis('DHBAVD-21','240424d')

h = figure
%boxplot([Hitrate_c,Hitrate,Hitrate2,Hitrate3],[zeros(1,length(Hitrate_c)),ones(1,length(Hitrate)),ones(1,length(Hitrate))*2,ones(1,length(Hitrate))*3])
%set(findobj(h, 'Tag', 'Whisker'), 'LineStyle', '-');
hold on
scatter([ones(1,length(Hitrate_c)),ones(1,length(Hitrate))*2,ones(1,length(Hitrate))*3,ones(1,length(Hitrate))*4],[Hitrate_c,Hitrate,Hitrate2,Hitrate3])
sigstar([1,2],ranksum(Hitrate_c,Hitrate))
sigstar([1,3],ranksum(Hitrate_c,Hitrate2))
sigstar([1,4],ranksum(Hitrate_c,Hitrate3))
sigstar([2,4],signrank(Hitrate,Hitrate3))
sigstar([2,3],signrank(Hitrate,Hitrate2))
sigstar([3,4],signrank(Hitrate2,Hitrate3))

line([ones(length(Hitrate),1)*2,ones(length(Hitrate),1)*3,ones(length(Hitrate),1)*4]',[Hitrate',Hitrate2',Hitrate3']','Color',[0.5,0.5,0.5,0.5])
xticks([1:4])
xticklabels({'Control_C21','DREADD_C21_1','DREADD_C21_2','DREADD_h2o'})
xlim([0,5])
setmyplot_balazs
ylabel('First new tone (7 Hz - go) hitrate')

h = figure
%boxplot([Hitrate_c,Hitrate,Hitrate2,Hitrate3],[zeros(1,length(Hitrate_c)),ones(1,length(Hitrate)),ones(1,length(Hitrate))*2,ones(1,length(Hitrate))*3])
%set(findobj(h, 'Tag', 'Whisker'), 'LineStyle', '-');
hold on
scatter([ones(1,length(HitrateF_c)),ones(1,length(HitrateF))*2,ones(1,length(HitrateF))*3,ones(1,length(HitrateF))*4],[HitrateF_c,HitrateF,HitrateF2,HitrateF3])
sigstar([1,2],ranksum(HitrateF_c,HitrateF))
sigstar([1,3],ranksum(HitrateF_c,HitrateF2))
sigstar([1,4],ranksum(HitrateF_c,HitrateF3))
sigstar([2,4],signrank(HitrateF,HitrateF3))
sigstar([2,3],signrank(HitrateF2,HitrateF3))
sigstar([3,4],signrank(HitrateF2,HitrateF3))

line([ones(length(HitrateF),1)*2,ones(length(HitrateF),1)*3,ones(length(HitrateF),1)*4]',[HitrateF',HitrateF2',HitrateF3']','Color',[0.5,0.5,0.5,0.5])
xticks([1:4])
xticklabels({'Control_C21','DREADD_C21_1','DREADD_C21_2','DREADD_h2o'})
xlim([0,5])
setmyplot_balazs
ylabel('4 Hz - go hitrate')

h = figure
hold on
scatter([ones(1,length(FarateF_c)),ones(1,length(FarateF))*2,ones(1,length(FarateF))*3,ones(1,length(FarateF))*4],[FarateF_c,FarateF,FarateF2,FarateF3])
sigstar([1,2],ranksum(FarateF_c,FarateF))
sigstar([1,3],ranksum(FarateF_c,FarateF2))
sigstar([1,4],ranksum(FarateF_c,FarateF3))
sigstar([2,4],signrank(FarateF,FarateF3))
sigstar([2,3],signrank(FarateF2,FarateF3))
sigstar([3,4],signrank(FarateF2,FarateF3))

line([ones(length(FarateF),1)*2,ones(length(FarateF),1)*3,ones(length(FarateF),1)*4]',[FarateF',FarateF2',FarateF3']','Color',[0.5,0.5,0.5,0.5])
xticks([1:4])
xticklabels({'Control_C21','DREADD_C21_1','DREADD_C21_2','DREADD_h2o'})
xlim([0,5])
setmyplot_balazs
ylabel('10 Hz - no-go false alarm rate')

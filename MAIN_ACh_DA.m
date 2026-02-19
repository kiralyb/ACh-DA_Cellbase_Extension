% initialization
choosecb('ACh_DA_FinalCellbase')
loadcb
%getvalue()
i_vta = find(~isnan([TheMatrix{:,25}]));
i_hdb = find(~isnan([TheMatrix{:,24}]));
i_da=find(abs([TheMatrix{:,25}])==4 | abs([TheMatrix{:,25}])==5);
i_da_i=find(abs([TheMatrix{:,25}])==5);
i_da_e=find(abs([TheMatrix{:,25}])==4);
i_ach=find(abs([TheMatrix{:,24}])==5);
i_da_tagged=i_vta([TheMatrix{i_vta,9}] < 0.01 & [TheMatrix{i_vta,12}] > 0.85);
i_ach_tagged=i_hdb([TheMatrix{i_hdb,7}] < 0.01 & [TheMatrix{i_hdb,11}] > 0.85);

stimmedAs_HDB = unique([TheMatrix{i_ach_tagged,3}])
stimmedAs_VTA = unique([TheMatrix{i_da_tagged,3}])

HDB_Clust = getvalue('HDB_Cluster_num');
VTA_Clust = getvalue('VTA_Cluster_num');


%% Behav
choosecb('ACh_DA_FinalCellbase')
MICEs = listtag('animal');MICEs([3])=[]; % !!!!

% Fig1c
psychometric_curves(MICEs)

% Sup1
[Hit, FA, RT, Mstate,trialnum,newcuenum,newtrialnum] = psy_behav_runner(MICEs,0)

% S1b
[RatHit1,RatFA1,RATRT1,RATstate1,RATtrialnum1] = psy_behav_runner(MICEs,1);
[RatHit2,RatFA2,RATRT2,RATstate2,RATtrialnum2] = psy_behav_runner(MICEs,2);

boxstat(cell2mat(RATRT1),cell2mat(RATRT2),'Parralel with new go tones','Parralel with new no-go tones',0.05,'paired')
boxstat(cell2mat(RatHit1),cell2mat(RatHit2),'Parralel with new go tones','Parralel with new no-go tones',0.05,'paired')
boxstat(cell2mat(RatFA1),cell2mat(RatFA2),'Parralel with new go tones','Parralel with new no-go tones',0.05,'paired')

%S1c
RTL = cell2mat(newtrialnum);
RNC = cell2mat(newcuenum);
RNC = RNC(RTL>100 & ~isnan(RNC));
nanmean(RNC) + 1;
figure
histogram(RNC,'Normalization','CDF')
setmyplot_balazs

%S1d
exampleID = LIST(200,:);
veiwresponsiveness(exampleID)

%S1e
viewlick(exampleID,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents','DeliverFeedback','Partitions','#Type3','window',[-2,3])
avg_viewlick(MICEs)

%% Photometry
choosecb('ACh_DA_photometry')
% Fig 1d
figure
viewphotometry('OAD16','220404','TriggerEvent','DeliverFeedback','Partitions','#Type','Signal',['dff_A']);
figure
viewphotometry('OAD16','220404','TriggerEvent','DeliverFeedback','Partitions','#Type','Signal',['dff_D']);

% Fig 1e
Photometry_GrandAVG('Soundtype','StimulusOn','A',1)
Photometry_GrandAVG('Soundtype','StimulusOn','D',1)

% Fig 1f
Photometry_GrandAVG('Type4','DeliverFeedback','A',1)
Photometry_GrandAVG('Type4','DeliverFeedback','D',1)

% Sup 3
Photometry_GrandAVG('Soundtype','StimulusOn','s405_A')
Photometry_GrandAVG('Soundtype','StimulusOn','s405_D')
Photometry_GrandAVG('Soundtype','StimulusOn','s465_A')
Photometry_GrandAVG('Soundtype','StimulusOn','s465_D') %file shitty 

% Sup 7a
Photometry_GrandAVG('CorrectRejection','StimulusOn','dff_A',0)
Photometry_GrandAVG('CorrectRejection','StimulusOn','dff_D',0)

% Fig 2i
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
(mean(lags(MI)))/sr
signrank(lags(MI))


%% Ephys


% Fig 1g
figure
viewcell2b('VVH11_200924a_9.1','TriggerName','BurstOn','eventtype','stimb','window',[-0.2 0.5],'dt',0.001,'sigma',0.001)

% Fig 1h
exampleIDs={'VVH6_191112a_8.2','WHV1_191123a_8.1','VVH7_200320a_16.1'}
for i = 1:length(exampleIDs)
figure
viewcell2b(exampleIDs{i},'TriggerName','DeliverFeedback','SortEvent','StimulusOn','ShowEvents','StimulusOn',...
    'eventtype','behav','window',[-1 1],'Partitions','#Type4')
end

% S4
i_vta = find(~isnan([TheMatrix{:,25}]));
i_hdb = find(~isnan([TheMatrix{:,24}]));
i_da_tagged=i_vta([TheMatrix{i_vta,9}] < 0.01 & [TheMatrix{i_vta,12}] > 0.85);
i_ach_tagged=i_hdb([TheMatrix{i_hdb,7}] < 0.01 & [TheMatrix{i_hdb,11}] > 0.85);
stimmedAs_HDB = unique([TheMatrix{i_ach_tagged,3}])
stimmedAs_VTA = unique([TheMatrix{i_da_tagged,3}])

% S4b
H_hdb_a = getvalue('Hindex',CELLIDLIST(i_hdb(ismember([TheMatrix{i_hdb,3}],stimmedAs_HDB))));
H_hdb_d = getvalue('Hindex_DA',CELLIDLIST(i_hdb(ismember([TheMatrix{i_hdb,3}],stimmedAs_HDB))));
H_vta_d = getvalue('Hindex_DA',CELLIDLIST(i_vta(ismember([TheMatrix{i_vta,3}],stimmedAs_VTA))));
H_vta_a = getvalue('Hindex',CELLIDLIST(i_vta(ismember([TheMatrix{i_vta,3}],stimmedAs_VTA))));
figure
subplot(2,1,1)
hold on
h1 = histogram(H_hdb_d,'BinEdges',0:0.01:1, 'FaceColor', 'm', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
h2 = histogram(H_hdb_a,'BinEdges',0:0.01:1, 'FaceColor', 'c', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
setmyplot_balazs
xlabel('H-index')
ylabel('Count')
subplot(2,1,2)
hold on
h1 = histogram(H_vta_d,'BinEdges',0:0.01:1, 'FaceColor', 'm', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
h2 = histogram(H_vta_a,'BinEdges',0:0.01:1, 'FaceColor', 'c', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
setmyplot_balazs
xlabel('H-index')
ylabel('Count')

% S4c
histcdf_plotter(getvalue('LR_PC',CELLIDLIST([i_hdb,i_vta])),getvalue('LR_PC',CELLIDLIST([i_ach_tagged,i_da_tagged])),linspace(0,0.15,100));
xlabel('L-ratio')
histcdf_plotter(getvalue('ID_PC',CELLIDLIST([i_hdb,i_vta])),getvalue('ID_PC',CELLIDLIST([i_ach_tagged,i_da_tagged])),linspace(0,500,100));
xlabel('Isolation distance')

% S4d
ChAT_R_L_J(CELLIDLIST(i_ach_tagged),'tag') 
ChAT_R_L_J(CELLIDLIST(i_da_tagged),'tagb')

% Sup 5
auROC_analysis2
auroc_k_elbow

% Fig 1i
trigger = 'StimulusOn';
partitions = '#Soundtype';
sigma = 0.02;
wind = [-0.5,1.5];
baslinewind=[-0.5,-0.0];
respwind_a = [0,0.4];
respwind_d = [0,0.4];
time_res = 0.001;
[P_ACh_SO,Latency_ACh_SO] = GAverage_PSTH(CELLIDLIST(i_ach),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_a,1)
[P_DA_i_SO,Latency_DA_i_SO] = GAverage_PSTH(CELLIDLIST(i_da_i),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_d,1)
[P_DA_e_SO,Latency_DA_e_SO] = GAverage_PSTH(CELLIDLIST(i_da_e),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_d,1)

% Fig 1j
trigger = 'DeliverFeedback';
partitions = '#Type4';
sigma = 0.01;
wind = [-1,1];
baslinewind=[-1,-0.6];
respwind_a = [0,0.4];
respwind_d = [0,0.4];
time_res = 0.001;
[p_ACh,~] = GAverage_PSTH(CELLIDLIST(i_ach),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_a,1)
[p_DA_i,~] = GAverage_PSTH(CELLIDLIST(i_da_i),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_d,-1)
[p_DA_e,~] = GAverage_PSTH(CELLIDLIST(i_da_e),trigger,partitions,sigma,time_res,wind,baslinewind,respwind_d,1)

% Fig 4d
latency_plot(Latency_ACh_SO(:,1),Latency_DA_i_SO(:,1),Latency_DA_e_SO(:,1))

trigger = 'DeliverFeedback';
[~,Latency_ACh_FA] = GAverage_PSTH(CELLIDLIST(i_ach),trigger,'FalseAlarm',sigma,time_res,wind,baslinewind,respwind_a,1,0)
[~,Latency_DA_i_FA] = GAverage_PSTH(CELLIDLIST(i_da_i),trigger,'FalseAlarm',sigma,time_res,wind,baslinewind,respwind_d,-1,0)
latency_plot(Latency_ACh_FA(:,1),Latency_DA_i_FA(:,1));

% S10h
[~,Latency_ACh_Hit] = GAverage_PSTH(CELLIDLIST(i_ach),trigger,'Hit',sigma,time_res,wind,baslinewind,respwind_a,1,0)
[~,Latency_DA_i_Hit] = GAverage_PSTH(CELLIDLIST(i_da_i),trigger,'Hit',sigma,time_res,wind,baslinewind,respwind_d,1,0)
[~,Latency_DA_e_Hit] = GAverage_PSTH(CELLIDLIST(i_da_e),trigger,'Hit',sigma,time_res,wind,baslinewind,respwind_d,1,0)
latency_plot(Latency_ACh_Hit(:,1),Latency_DA_i_Hit(:,1),Latency_DA_e_Hit(:,1))

% S7b
GAverage_PSTH(CELLIDLIST(i_ach),trigger,'#CorrectRejection',sigma,time_res,wind,baslinewind,respwind_a,1,0)
GAverage_PSTH(CELLIDLIST(i_da_i),trigger,'#CorrectRejection',sigma,time_res,wind,baslinewind,respwind_d,1,0)
GAverage_PSTH(CELLIDLIST(i_da_e),trigger,'#CorrectRejection',sigma,time_res,wind,baslinewind,respwind_d,1,0)

%% Signal Correlation

% Fig 2a, S6a
drifting_ach = signal_correlation_Hz(i_ach);
drifting_da = signal_correlation_Hz(i_da);
drifting = [drifting_ach,drifting_da];

% S6b
signal_correlation_Hz(i_da_i);
signal_correlation_Hz(i_da_e);

% Fig 2b
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'cue',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'cue',2)

% S6c
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'cue',3)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'cue',4)

% S6d
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'rewp',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'rewp',2)
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'rewp',3)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'rewp',4)

%S6e
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'pun',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'pun',2)
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'pun',3,-1)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'pun',4)

%S6f
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'rewr',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'rewr',2)
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'rewr',3)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'rewr',4)

%S6g
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'rewrn',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'rewrn',2)
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'rewrn',3)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'rewrn',4)

%S6h
figure
signal_correlation_pf(i_ach(~ismember(i_ach,drifting)),'punrn',1)
signal_correlation_pf(i_da(~ismember(i_da,drifting)),'punrn',2)
signal_correlation_pf(i_da_i(~ismember(i_da_i,drifting)),'punrn',3,-1)
signal_correlation_pf(i_da_e(~ismember(i_da_e,drifting)),'punrn',4)

%% Change Dynamics
[crosspairs, ie_pairs] = pairIDfinder(i_ach,i_da,i_da_i);
% Fig 2c
changedyn(crosspairs(45,:),1)

% Fig S8a
changedyn({'VVH6_191113a_1.2'},1)
changedyn({'VVH11_200923a_9.2'},1)
changedyn({'VVH11_201006a_9.1'},1)

% Fig 2d
for j =1:length(crosspairs)
    [R_a(j), P_a(j),R_a_b(j), P_a_b(j),R_da(j), P_da(j), R_d_b(j), P_d_b(j),...
        linparams_a(j,:),R2_a(j),p_a(j),linparams_da(j,:),R2_da(j),p_da(j),...
        R_a_d(j), P_a_d(j),R_da_d(j), P_da_d(j),R_par(j,:,:),P_par(j,:,:),...
        VAR(j),B(:,j),B_p(:,j),B_pc(:,j),Bd(:,j),Bd_p(:,j),Bc(:,j),Bc_p(:,j),...
        trialnumber(j),linparams_cross(j,:),linparams_cross_reverse(j,:),R2_cross(j)...
        ,R2_cross_reverse(j)] ...
    = changedyn(crosspairs(j,:),0);
end

figure
for ie = 0:1
    ax(1+3*ie)=subplot(2,3,1+3*ie)
    boxstat(B(ie_pairs==ie),Bd(ie_pairs==ie),'ach','da',0.05,'paired',ax(1))
    ylim([0.5,1])
    ax(2+3*ie)=subplot(2,3,2+3*ie)
    boxstat(Bd(ie_pairs==ie),Bc(ie_pairs==ie),'da','combined',0.05,'paired',ax(2))
    ax(3+3*ie)=subplot(2,3,3+3*ie)
    boxstat(B(ie_pairs==ie),Bc(ie_pairs==ie),'ach','combined',0.05,'paired',ax(3))
    ylim([0.5,1])
    linkaxes([ax], 'y')
end

figure
    ax(1)=subplot(1,3,1)
    boxstat(B,Bd,'ach','da',0.05,'paired',ax(1))
    ax(2)=subplot(1,3,2)
    boxstat(Bd,Bc,'da','combined',0.05,'paired',ax(2))
    ax(3)=subplot(1,3,3)
    boxstat(B,Bc,'ach','combined',0.05,'paired',ax(3))
    linkaxes([ax], 'y')

[dapairs, die_pairs] = pairIDfinder(i_da_e,i_da,i_da_i);
for j =1:length(dapairs)
    [~, ~,~, ~,~, ~,~, ~,...
        ~,~,~,~,~,~,...
        ~,~,~,~,~,~,...
        ~,~,~,~,~,~,Bc_DANS(:,j),~,...
        ~,~,~,~...
        ,~] ...
    = changedyn(dapairs(j,:),0);
end
boxstat(Bc,Bc_DANS,'ACh+DA','T1-DAN + T2-DAN');

% Fig 2f right
for i = 1:2
    boxstat(linparams_cross(ie_pairs==i-1,2),linparams_cross_reverse(ie_pairs==i-1,2),'ACh->DA','DA->ACh',0.05,'paired')
end

% S8c
H=figure
maximize_figure(H);
H = subplot(1,6,1)
boxstat(R_par(ie_pairs==0,1,3),R_par(ie_pairs==0,3,1),'ach','controlled_e',0.05,'paired',H)
H = subplot(1,6,2)
boxstat(R_par(ie_pairs==1,1,3),R_par(ie_pairs==1,3,1),'ach','controlled_i',0.05,'paired',H)
H = subplot(1,6,3)
boxstat(R_par(ie_pairs==1,2,3),R_par(ie_pairs==1,3,2),'da_i','controlled',0.05,'paired',H)
H = subplot(1,6,4)
boxstat(R_par(ie_pairs==0,2,3),R_par(ie_pairs==0,3,2),'da_e','controlled',0.05,'paired',H)
H = subplot(1,6,5)
boxstat(R_par(ie_pairs==0,1,2),R_par(ie_pairs==0,2,1),'ach_da_e','controlled',0.05,'paired',H)
H = subplot(1,6,6)
boxstat(R_par(ie_pairs==1,1,2),R_par(ie_pairs==1,2,1),'ach_da_i','controlled',0.05,'paired',H)

% S8b right
boxstat(R_par(ie_pairs==1,1,2),R_par(ie_pairs==0,1,2),'da_i','da_e')

% Fig 2e
inx = i_ach
numcells = length(inx);
j = 1;
for k =1:numcells
    try % ex k = 55, k = 39 81  (old 45) 
    [R_ach(j), P_ach(j),R_ach_b(j),P_ach_b(j),~, ~,~, ~,linparams_ach(j,:),R2_ach(j),p_ach(j),~,~,~,R_ach_d(j),P_ach_d(j),~,~,~,~,var_ach(j),B_ach(j),B_p_ach(j),~,~,~,~,~,trialnum_ach(j)] = changedyn({CELLIDLIST{inx(k)}},0);
    j = j + 1;
    catch
    k
    end
end

inx = i_da_i;
numcells = length(inx);
j = 1;
for k =1:numcells
    try % k = 136
    [R_da_i(j), P_da_i(j),R_da_i_b(j),P_da_i_b(j),~, ~,~, ~,linparams_da_i(j,:),R2_da_i(j),p_da_i(j),~,~,~,R_da_i_d(j),P_da_i_d(j),~,~,~,~,var_da_i(j),B_da_i(j),B_p_da_i(j),~,~,~,~,~,trialnum_da_i(j)] = changedyn({CELLIDLIST{inx(k)}},0);
    j = j + 1;
    end
end

inx = i_da_e;
numcells = length(inx);
j = 1;
for k =1:numcells
    try % k = 61k = 61
    [R_da_e(j), P_da_e(j),R_da_e_b(j),P_da_e_b(j),~, ~,~, ~,linparams_da_e(j,:),R2_da_e(j),p_da_e(j),~,~,~,R_da_e_d(j),P_da_e_d(j),~,~,~,~,var_da_e(j),B_da_e(j),B_p_da_e(j),~,~,~,~,~,trialnum_da_e(j)] = changedyn({CELLIDLIST{inx(k)}},1);
    j = j + 1;
    end
end

% S8b left
labels = {'ACh', 'DA_i' ,'DA_e'}
boxstat3(R_ach,R_da_i,R_da_e,labels)
ylabel('trend R')
% Fig 2e
boxstat3(R_ach_d,R_da_i_d,R_da_e_d,labels)
ylabel('Diff R')

% Fig 2f
labels = {'ACh', 'DA_i' ,'DA_e'}
boxstat3(linparams_ach(:,2),linparams_da_i(:,2),linparams_da_e(:,2),labels)
signrank(linparams_ach(:,2),0.1)
signrank(linparams_da_i(:,2),0.0)
ylabel('dt')
% S8d
boxstat3(R2_ach,R2_da_i,R2_da_e,labels)
ylabel('R2')

%% Interaction
% CCGs
% F2g
path = [getpref('cellbase','datapath'),'CCG_final\'];
mkdir(path);
ccg(crosspairs,0.5,'issave',false,'resdir',[path,'All'],'segfilter','@stim_excl_dual','filterinput',{'light_activation_duration',[0 1],'margins',[0,0]})
CCG = CCG_normalizer([path,'ALL\CCG_matrices.mat']);
figure
subplot(1,2,1)
ccg_plot(CCG(ie_pairs == 0,:))
subplot(1,2,2)
ccg_plot(CCG(ie_pairs == 1,:))

% F2h
ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Cue'],'segfilter','@cue_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Hit'],'segfilter','@Hit_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
ccg(crosspairs,0.25,'issave',false,'resdir',[path,'FA'],'segfilter','@FA_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)

AVG_CCG_plot(path,ie_pairs)

% F3a & % S9a-b

if ~exist([path,'ISI\CCG_matrices.mat'])
ccg(CELLIDLIST(~isnan(HDB_Clust)|~isnan(VTA_Clust)),...
    0.5,'issave',false,'resdir',[path,'ISI'],'segfilter','@stimfb_excl_dual','filterinput',{'light_activation_duration',[0 1],'feedback_duration',[-0.6 0.6],'margins',[0 0]})
end
CCG_matrix_plotter(path,HDB_Clust,VTA_Clust);
CCG_matrix_plotter(path,HDB_Clust,HDB_Clust);
CCG_matrix_plotter(path,VTA_Clust,VTA_Clust);

% F2h


% JPSTH
% S8e
jpsth(crosspairs,ie_pairs,'StimulusOn','all',[0,0.6],25)
% S8f
jpsth(crosspairs,ie_pairs,'DeliverFeedback','#Hit',[0,0.6],25)
% S8g
jpsth(crosspairs,ie_pairs,'DeliverFeedback','#FalseAlarm',[0,0.6],25)

% Noise correlation
% Fig 2j
noise_correlation(crosspairs(25),1,1) %ex25
noise_correlation(crosspairs(29),1,1) %ex29

% Fig 2k
[PP,RR] = noise_correlation(crosspairs(ie_pairs==1,:),1,0) %ex25
x=histogram(PP,0:0.05:1);
PH(1,:) = x.Values;
xx=histogram(RR(PP<0.05),-0.5:0.1:1);
RH(1,:) = xx.Values;
x=histogram(RR,-0.5:0.1:1);
RH(2,:) = x.Values-RH(1,:);
[PP,RR] = noise_correlation(crosspairs(ie_pairs==0,:),1,0)  % ex 29
x=histogram(PP,0:0.05:1);
PH(2,:) = x.Values;
xx=histogram(RR(PP<0.05),-0.5:0.1:1);
RH(3,:) = xx.Values;
x=histogram(RR,-0.5:0.1:1);
RH(4,:) = x.Values-RH(3,:);

figure
subplot(1,3,1)
bar(0.025:0.05:(1-0.025),PH(1:2,:)','stacked')
legend({'a - d_i','a - d_e'})
setmyplot_balazs
subplot(1,3,2)
bar(-0.45:0.1:0.95,RH([1,3],:)','stacked')
setmyplot_balazs
subplot(1,3,3)
bar(-0.45:0.1:0.95,RH([1:4],:)','stacked')
xlim([-1,1])
legend({'a - d_i','a - d_e'})
setmyplot_balazs


%% Optogenetic activation
miceID = getvalue('RatID')
% Fig 3c and S9c
colors={'c','y','m','r'}
nums = [5,3,5,4]
area = [1,1,2,2]
figure
hold on
for i = 1:4
    if area(i) == 1
        IDS = CELLIDLIST(abs(HDB_Clust)==nums(i));
        IDS = IDS(ismember(miceID(abs(HDB_Clust)==nums(i)),stimmedAs_HDB)); 
    else
        IDS = CELLIDLIST(abs(VTA_Clust)==nums(i));
        IDS = IDS(ismember(miceID(abs(VTA_Clust)==nums(i)),stimmedAs_HDB));
    end
p(i) = optoPSTH(IDS,'stim',colors{i},i-1)
end
legend({'ACh','pGABA','DA_i','DA_e'})

% S9c
IDS = CELLIDLIST(abs(VTA_Clust) == 4 | abs(VTA_Clust) == 5);
IDS = IDS(ismember(miceID(abs(VTA_Clust) == 4 | abs(VTA_Clust) == 5),stimmedAs_VTA));
figure
pvta(1) = optoPSTH(IDS,'stimb',{'m'},1)
legend({'DA'})

% S9d
colors={'c','y','k'}
nums = [5,3,1]
area = [1,1,1]
figure
hold on
for i = 1:3
    if area(i) == 1
        IDS = CELLIDLIST(abs(HDB_Clust)==nums(i));
        IDS = IDS(ismember(miceID(abs(HDB_Clust)==nums(i)),stimmedAs_VTA)); 
    else
        IDS = CELLIDLIST(abs(VTA_Clust)==nums(i));
        IDS = IDS(ismember(miceID(abs(VTA_Clust)==nums(i)),stimmedAs_VTA));
    end
pvta(i + 1)=optoPSTH(IDS,'stimb',colors{i},i-1)
end
legend({'ACh','pHGABA','DA_i'})

%% Chemogenetic supression
DREADD

 function MAIN_ACh_DA
%MAIN_ACH_DA Main wrapper for the figures presented in
%   'Cholinergic–dopaminergic interplay underlies prediction error
%   broadcasting' by Király et al. (2026)
%
%   The code performs:
%       - Behavioral analyses in a psychometric task
%       - Photometry analyses of ACh and DA release
%       - Single unit analyses of HDB and VTA neurons
%       - Signal correlation analyses of BFCNs and DANS in broadcastign RPE 
%       - Change dynamics analyses of BFCNs and DANs during learning
%       - Interaction analyses of BFCNs and DANs 
%       - Optogenetic activation analyses of cross-BFCN-DAN effects 
%       - Chemogenetic suppression analyses of BFCNs effect on learning & DA
%
%   Required datasets:
%       - 'ACh_DA_FinalCellbase' (Single unit recordings from the HDB and VTA with optogenetical tagging of BFCNs and DANs)
%       - 'ACh_DA_photometry' (ACh and DA release measured from the BLA and the VS)
%       - 'ACh_DA_psychometric' (Chemogenetics supression of BFCNs with ACh (BLA) and DA (VS) release measurements)

%   Bálint Király
%   Division of Neurophysiology
%   Medical University of Vienna
%   balint.kiraly@meduniwien.ac.at
%   23-Feb-2026


%% Behavioral analysis in the psychometric learning task.
% Figures 1c & S1
%--------------------------------------------------------------------------

% Initialization
choosecb('ACh_DA_FinalCellbase')
MICEs = listtag('animal');
LIST = listtag('session');

% Psychometric learning curves
% Fig 1c
psychometric_curves(MICEs)

% Fixed association performance
% Fig S1a
[~,~,~,~,~, newcuenum, newtrialnum] = psy_behav_runner(MICEs,0);

% Fig S1b
[MiceHit1, MiceFA1] = psy_behav_runner(MICEs,1);
[MiceHit2, MiceFA2] = psy_behav_runner(MICEs,2);

boxstat(cell2mat(MiceHit1),cell2mat(MiceHit2), ...
    'Parralel with new go tones','Parralel with new no-go tones', ...
    0.05,'paired')

boxstat(cell2mat(MiceFA1),cell2mat(MiceFA2), ...
    'Parralel with new go tones','Parralel with new no-go tones', ...
    0.05,'paired')

% New assocation formation rate
% Fig S1c
RTL = cell2mat(newtrialnum);
RNC = cell2mat(newcuenum);
RNC = RNC(RTL > 100 & ~isnan(RNC)); % only consider sufficiantly long sessions (>100 trials)
figure
histogram(RNC,'Normalization','CDF')

% Responsivness in example session and on average 
exampleID = {'VVH2','190608a'};
% Fig S1d
veiwresponsiveness(exampleID)
 
% Fig S1e
viewlick(exampleID,'TriggerName','StimulusOn','SortEvent','TrialStart',...
    'eventtype','behav','ShowEvents','DeliverFeedback',...
    'Partitions','#Type3','window',[-2,3])
avg_viewlick(MICEs)

%% Fiber photometry analyses of ACh and DA release during task completition.
% Figures 1d–f, S3, S7a & 2i
%--------------------------------------------------------------------------

%Initialization
choosecb('ACh_DA_photometry')

% Example session with parallel ACh and DA release measurments 
% Fig 1d
figure
viewphotometry('OAD16','220404','TriggerEvent','DeliverFeedback','Partitions','#Type4','Signal','dff_A');
figure
viewphotometry('OAD16','220404','TriggerEvent','DeliverFeedback','Partitions','#Type4','Signal','dff_D');

% Average ACh and DA release across animals and events 
% Fig 1e - Cue aligned
Photometry_GrandAVG('Soundtype','StimulusOn','A',1)
Photometry_GrandAVG('Soundtype','StimulusOn','D',1)

% Fig 1f - Feedback aligned
Photometry_GrandAVG('Type4','DeliverFeedback','A',1)
Photometry_GrandAVG('Type4','DeliverFeedback','D',1)

% Fig S3 - Ligand dependent and isosbestic signals
Photometry_GrandAVG('Soundtype','StimulusOn','s405_A',0)
Photometry_GrandAVG('Soundtype','StimulusOn','s405_D',0)
Photometry_GrandAVG('Soundtype','StimulusOn','s465_A',0)
Photometry_GrandAVG('Soundtype','StimulusOn','s465_D',0)

% Fig S7a - Correct rejection
Photometry_GrandAVG('CorrectRejection','StimulusOn','dff_A',0)
Photometry_GrandAVG('CorrectRejection','StimulusOn','dff_D',0)

% Cross-correlograms of ACh and DA release during the task
% Fig 2i
combined_CCG

%% Analysis of HDB and VTA single unit activities.
% Figure 1g-j, S4, S5, S7b, 4d & S10h
%-------------------------------------------------
% Initialization
choosecb('ACh_DA_FinalCellbase')
loadcb

% Example BFCN and DA neurons
% Fig 1g - optogenetical tagging
figure
viewcell2b('VVH11_200924a_9.1','TriggerName','BurstOn','eventtype','stimb','window',[-0.2 0.5],'dt',0.001,'sigma',0.001)

% Fig 1h - feedback aligned
exampleIDs={'VVH6_191112a_8.2','WHV1_191123a_8.1','VVH7_200320a_16.1'};
for i = 1:length(exampleIDs)
    figure
    viewcell2b(exampleIDs{i},'TriggerName','DeliverFeedback','SortEvent','StimulusOn','ShowEvents','StimulusOn',...
        'eventtype','behav','window',[-1 1],'Partitions','#Type4')
end

% Optogenetical tagging properties of BFCNs and DANs
% Fig S4
i_hdb = find(~isnan(getvalue('HDB_Cluster_num')));
i_vta = find(~isnan(getvalue('VTA_Cluster_num')));
i_da_tagged=i_vta(getvalue('Hindex_DA',CELLIDLIST(i_vta)) < 0.01 & getvalue('R_DA',CELLIDLIST(i_vta)) > 0.85);
i_ach_tagged=i_hdb(getvalue('Hindex',CELLIDLIST(i_hdb)) < 0.01 & getvalue('R',CELLIDLIST(i_hdb)) > 0.85);
stimmedAs_HDB = unique(getvalue('RatID',CELLIDLIST(i_ach_tagged)));
stimmedAs_VTA = unique(getvalue('RatID',CELLIDLIST(i_da_tagged)));

% S4b - H index
H_hdb_a = getvalue('Hindex',CELLIDLIST(i_hdb(ismember(getvalue('RatID',CELLIDLIST(i_hdb)),stimmedAs_HDB))));
H_hdb_d = getvalue('Hindex_DA',CELLIDLIST(i_hdb(ismember(getvalue('RatID',CELLIDLIST(i_hdb)),stimmedAs_HDB))));
H_vta_d = getvalue('Hindex_DA',CELLIDLIST(i_vta(ismember(getvalue('RatID',CELLIDLIST(i_vta)),stimmedAs_VTA))));
H_vta_a = getvalue('Hindex',CELLIDLIST(i_vta(ismember(getvalue('RatID',CELLIDLIST(i_vta)),stimmedAs_VTA))));
combHhist(H_hdb_d,H_hdb_a)
combHhist(H_vta_d,H_vta_a)

% S4c - L-ratio and ID
histcdf_plotter(getvalue('LR_PC',CELLIDLIST([i_hdb;i_vta])),getvalue('LR_PC',CELLIDLIST([i_ach_tagged;i_da_tagged])),linspace(0,0.15,100));
xlabel('L-ratio')
histcdf_plotter(getvalue('ID_PC',CELLIDLIST([i_hdb;i_vta])),getvalue('ID_PC',CELLIDLIST([i_ach_tagged;i_da_tagged])),linspace(0,500,100));
xlabel('Isolation distance')

% S4d - Latency and jitter
ChAT_R_L_J(CELLIDLIST(i_ach_tagged),'tag')
ChAT_R_L_J(CELLIDLIST(i_da_tagged),'tagb')

% Cluster analysis of PSTH characteristics during task execution
% Fig S5
populationclust_elbow('HDB')
populationclust_elbow('VTA')
i_da_i = find(abs(getvalue('VTA_Cluster_num')) == 5);
i_da_e = find(abs(getvalue('VTA_Cluster_num')) == 4);
i_da = [i_da_i;i_da_e];
i_ach = find(abs(getvalue('HDB_Cluster_num')) == 5);

% Average PSTHs and response latencies across neurons types and events
% Fig 1i - Cue aligned
cellgroups = {i_ach, i_da_i, i_da_e};
window = [-0.5,0.5];
baseline_window = [-0.5,-0.0];
sigma = 0.02;
Latency_SO = cell(1,length(cellgroups));
for i_group = 1:length(cellgroups)
    [Latency_SO{i_group}] = GAverage_PSTH(CELLIDLIST(cellgroups{i_group}),'StimulusOn','#Soundtype',sigma,window,baseline_window);
end

% Fig 1j - Feedback aligned
window = [-0.5,0.5];
baseline_window = [-0.5,-0.0];
sigma = 0.02;
for i_group = 1:length(cellgroups)
    GAverage_PSTH(CELLIDLIST(cellgroups{i_group}),'DeliverFeedback','#Type4',sigma,window,baseline_window);
end

% Fig 4d - Punishment response latency
boxstat3(Latency_SO{1}(:,1),Latency_SO{2}(:,1),Latency_SO{3}(:,1),{'ACh', 'DA_i' ,'DA_e'})
Latency_ACh_FA = GAverage_PSTH(CELLIDLIST(i_ach),'DeliverFeedback','#FalseAlarm',sigma,window,baseline_window,1,0);
Latency_DA_i_FA = GAverage_PSTH(CELLIDLIST(i_da_i),'DeliverFeedback','#FalseAlarm',sigma,window,baseline_window,-1,0);
boxstat(Latency_ACh_FA(:,1),Latency_DA_i_FA(:,1),'ACh','DA_i');

% S10h - Reward response latency
Latency_Hit = cell(1,length(cellgroups));
for i_group = 1:length(cellgroups)
    [Latency_Hit{i_group}] = GAverage_PSTH(CELLIDLIST(cellgroups{i_group}),'DeliverFeedback','#Hit',sigma,window,baseline_window);
end
boxstat3(Latency_Hit{1}(:,1),Latency_Hit{2}(:,1),Latency_Hit{3}(:,1),{'ACh', 'DA_i' ,'DA_e'})

% S7b - Correct rejections
for i_group = 1:length(cellgroups)
    GAverage_PSTH(CELLIDLIST(cellgroups{i_group}),'DeliverFeedback','#CorrectRejection',0.01,[-1,0.5],[-1,-0.6]);
end

%% Signal correlations of BFCNs and DANs in broadcasting RPE
% Figures 2a-b & S6
%-----------------------------------------------------------

% Initialization
cellgroups = {i_ach, i_da, i_da_i, i_da_e};
peakorthrough = [1,1,-1,1]; % (-1) indicate that Type 1 DANs are supressed by punishments 
drifting = [];

% Task difficulty correlations
% Fig 2a & S6a-b
for i_group = 1:length(cellgroups)
    d = signal_correlation_Hz(cellgroups{i_group});
    drifting = [drifting,d];
end

% Learning progress correlations
% Fig 2b & S6c - Cue
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'cue',i_group)
end

% Fig S6d - Reward
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'rewp',i_group)
end

% Fig S6e - Punishment
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'pun',i_group,peakorthrough(i_group))
end

% Temporal aspects of RPE
% Fig S6f - Reaction time
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'rewr',i_group)
end

% Fig S6g - Reward delay
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'rewrn',i_group)
end

% Fig S6h - Punishment delay 
figure
for i_group = 1:length(cellgroups)
    signal_correlation_pf(setdiff(cellgroups{i_group}, drifting),'punrn',i_group,peakorthrough(i_group))
end

%% Change Dynamics of BFCN and DAN activity during learning
% Figures 2c-f & S8a-e
%-----------------------------------------------------------
% Initialize list of simulatneously recorded BFCN-DAN pairs
[crosspairs, ie_pairs] = pairIDfinder(i_ach,i_da_i,i_da_e); 

% Example neurons' new cue responsivenes trends
% Fig 2c
changedyn2({'VVH6_191127a_8.1' , 'VVH6_191127a_16.1' },1)

% Fig S8a
changedyn2({'VVH6_191113a_1.2'},1)
changedyn2({'VVH11_200923a_9.2'},1)
changedyn2({'VVH11_201006a_9.1'},1)

% Logistic regression of new cue responses
% Fig 2d
for j =1:length(crosspairs)
    [B_ach(:,j),B_da(:,j),Bc(:,j)] = lin_reg_pair(crosspairs(j,:));
end
figure
ax(1)=subplot(1,3,1);
boxstat(B_ach,B_da,'ach','da',0.05,'paired',ax(1))
ax(2)=subplot(1,3,2);
boxstat(B_da,Bc,'da','combined',0.05,'paired',ax(2))
ax(3)=subplot(1,3,3);
boxstat(B_ach,Bc,'ach','combined',0.05,'paired',ax(3))
linkaxes(ax, 'y')

% Fig S8b
figure
for ie = 0:1
    ax(1+3*ie)=subplot(2,3,1+3*ie);
    boxstat(B_ach(ie_pairs==ie),B_da(ie_pairs==ie),'ach','da',0.05,'paired',ax(1))
    ax(2+3*ie)=subplot(2,3,2+3*ie);
    boxstat(B_da(ie_pairs==ie),Bc(ie_pairs==ie),'da','combined',0.05,'paired',ax(2))
    ax(3+3*ie)=subplot(2,3,3+3*ie);
    boxstat(B_ach(ie_pairs==ie),Bc(ie_pairs==ie),'ach','combined',0.05,'paired',ax(3))
end
linkaxes(ax, 'y')

% Fig 2d right
dapairs = pairIDfinder(i_da_e,i_da_i); % find pairs across DA subtypes
for j =1:length(dapairs)
    [~,~,Bc_DANS(:,j)] = lin_reg_pair(dapairs(j,:));
end
boxstat(Bc,Bc_DANS,'ACh+DA','T1-DAN + T2-DAN');

% Correlations with behavioral trends
for k = 1:length(i_ach)
    [R_ach(k), R_ach_d(k), linparams_ach(k,:), R2_ach(k)] = changedyn2(CELLIDLIST(i_ach(k)),0);
end
for k = 1:length(i_da_i)
    [R_da_i(k), R_da_i_d(k), linparams_da_i(k,:), R2_da_i(k)] = changedyn2(CELLIDLIST(i_da_i(k)),0);
end
for k = 1:length(i_da_e)
    [R_da_e(k), R_da_e_d(k), linparams_da_e(k,:), R2_da_e(k)] = changedyn2(CELLIDLIST(i_da_e(k)),0);
end

% Fig S8c left - correlation of trends
labels = {'ACh', 'DA_i' ,'DA_e'};
boxstat3(R_ach,R_da_i,R_da_e,labels)
ylabel('trend R')
% Fig 2e - correlation of updates
boxstat3(R_ach_d,R_da_i_d,R_da_e_d,labels)
ylabel('Diff R')

% Linear regression of behavioral responses with temporal shifts allowed
% Fig 2f - temporal shift 
boxstat3(linparams_ach(:,2),linparams_da_i(:,2),linparams_da_e(:,2),labels)
ylabel('dt')
signrank(linparams_ach(:,2),0.1)
% Fig S8e - explanatory power
boxstat3(R2_ach,R2_da_i,R2_da_e,labels)
ylabel('R2')

% Cross-system trend correlations
% Fig 2f
for j =1:length(crosspairs)
    [~,~,~,~,R_par(j,:,:),linparams_cross(j,:),linparams_cross_reverse(j,:)] ...
        = changedyn2(crosspairs(j,:),0);
end

% Fig 2f right - linear regression
for i = 1:2
    boxstat(linparams_cross(ie_pairs==i-1,2),linparams_cross_reverse(ie_pairs==i-1,2),'ACh->DA','DA->ACh',0.05,'paired')
end

% Fig S8c right - correlation of new cue responses
boxstat(R_par(ie_pairs==1,1,2),R_par(ie_pairs==0,1,2),'da_i','da_e')

% Fig S8d - partial correlations
boxstat(R_par(ie_pairs==0,1,3),R_par(ie_pairs==0,3,1),'ACh','DA_e controlled',0.05,'paired')
boxstat(R_par(ie_pairs==1,1,3),R_par(ie_pairs==1,3,1),'ACh','DA_i controlled',0.05,'paired')
boxstat(R_par(ie_pairs==1,2,3),R_par(ie_pairs==1,3,2),'DA_i','ACh controlled',0.05,'paired')
boxstat(R_par(ie_pairs==0,2,3),R_par(ie_pairs==0,3,2),'DA_e','ACh controlled',0.05,'paired')
boxstat(R_par(ie_pairs==0,1,2),R_par(ie_pairs==0,2,1),'AChxda_e','behav. controlled',0.05,'paired')
boxstat(R_par(ie_pairs==1,1,2),R_par(ie_pairs==1,2,1),'AChxda_i','behav. controlled',0.05,'paired')
%% Interactions between BFCNs and DANs
% Figures 2g,j,k, 3a, S8e-g & S9-ab, 
%-------------------------------------
% Initialization
path = [getpref('cellbase','datapath'),'CCG_final\'];

% Crosscorrelograms(CCG)
% Fig 2g BFCN-DA CCGs
mkdir(path);
if ~exist([path,'ALL\CCG_matrices.mat'],'file')
    ccg(crosspairs,0.5,'issave',false,'resdir',[path,'All'],'segfilter','@stim_excl_dual','filterinput',{'light_activation_duration',[0 1],'margins',[0,0]})
end
CCG = CCG_normalizer([path,'ALL\CCG_matrices.mat']);
figure
subplot(1,2,1)
ccg_plot(CCG(ie_pairs == 0,:))
subplot(1,2,2)
ccg_plot(CCG(ie_pairs == 1,:))

% Fig 2h - event restricted BFCN-DA CCGs 
if ~exist([path,'Cue\CCG_matrices.mat'],'file')
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Cue'],'segfilter','@cue_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Hit'],'segfilter','@Hit_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'FA'],'segfilter','@FA_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
end
AVG_CCG_plot(path,ie_pairs)

% Fig 3a & S9a-b - CCG matrix of cell clusters
HDB_Clust = abs(getvalue('HDB_Cluster_num'));
VTA_Clust = abs(getvalue('VTA_Cluster_num'));
if ~exist([path,'ISI\CCG_matrices.mat'],'file')
    ccg(CELLIDLIST(~isnan(HDB_Clust)|~isnan(VTA_Clust)),0.5,'issave',false,'resdir',[path,'ISI'],...
        'segfilter','@stimfb_excl_dual','filterinput',{'light_activation_duration',[0 1],'feedback_duration',[-0.6 0.6],'margins',[0 0]})
end
CCG_matrix_plotter(path,HDB_Clust,VTA_Clust);
CCG_matrix_plotter(path,HDB_Clust,HDB_Clust);
CCG_matrix_plotter(path,VTA_Clust,VTA_Clust);

% Joint PSTHs
% Fig S8e - cue
jpsth(crosspairs,ie_pairs,'StimulusOn','all')
% Fig S8f - Reward
jpsth(crosspairs,ie_pairs,'DeliverFeedback','#Hit')
% Fig S8g - Punishment
jpsth(crosspairs,ie_pairs,'DeliverFeedback','#FalseAlarm')

% Noise correlation
% Fig 2j - examples
noise_correlation({'VVH9_200622b_8.1','VVH9_200622b_9.1'})
noise_correlation({'VVH6_191114a_3.1','VVH6_191114a_16.2'})

% Fig 2k - distribution
noise_correlation_hist(crosspairs,ie_pairs);

%% Cross-system effects of optogenetic activation
% Figures 3c & S9c-d
%------------------------------------------------

% Fig 3c and S9c - BFCN stimulation
optoPSTH_plot([5,3,5,4],[1,1,2,2],stimmedAs_HDB,{'c','y','m','r'},{'ACh','pBGABA','DA_i','DA_e'})

% Fig S9d - DAN stimulation
optoPSTH_plot([6,5,3,1],[2,1,1,1],stimmedAs_VTA,{'m','c','y','k'},{'DA','ACh','pBGABA','BF1'})

%% Chemogenetic supression of BFCNs
% Figures 4 & S10
%----------------------------------

% Initialization
choosecb('Cellbase_psychometric'); loadcb;
AnimalIDs = CELLIDLIST;
animal_group = getvalue('animal_group');
ach_mask = getvalue('isBLA');
da_mask = getvalue('isVS');

% Behavioral effects
% Fig 4b - First new association
firstNovelBehav_DREADD(animal_group)

%Psychometric learning curves
% Fig 4c - Treated mice
[DIFF_D] = psychometric_DREADD_COMP(AnimalIDs(animal_group==1),ones(1,sum(animal_group==1))*4);
% S10c - Control mice
[DIFF_C] = psychometric_DREADD_COMP(AnimalIDs(animal_group==0),ones(1,sum(animal_group==0))*4);
% S10d - Difference
boxstat(DIFF_C,DIFF_D,'control','dreadd')

% Fixed associations
% S10e 
for i = 1:length(CELLIDLIST)
    load([getpref('cellbase','datapath'),'\',AnimalIDs{i},'\',AnimalIDs{i},'.mat'],'Hit','FalseAlarm','RT','type_SO');    
    for type = 1:2 % extract C21 (1) and Control (0) days
    MHit(i,type) = mean(Hit(type_SO{1}==type-1));
    MFA(i,type) = mean(FalseAlarm(type_SO{1}==type-1));
    MRT(i,type) = mean(RT(type_SO{1}==type-1));
    end
end
fix_DREADD_comp(MRT,animal_group)
ylabel('Reaction time (s)')
fix_DREADD_comp(MHit,animal_group)
ylabel('Hit %')
fix_DREADD_comp(MFA,animal_group)
ylabel('FA %')

% Effect of BFCN supression on neuromodulator release
% Fig 4e, 4g, S10f - first new association
firstnew_photometry_run(AnimalIDs,ach_mask,da_mask);

% 4f,h,i,j, S10g,i - long term
dreadd_photometry_comp(AnimalIDs,ach_mask,da_mask);


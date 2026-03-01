function preproc_ACh_DA

%% FiberPhotometry PSTHs
choosecb('ACh_DA_photometry')

animalID = {'OAD2','OAD3','OAD4','OAD5','OAD6','OAD7','OAD15','OAD16','OAD17'...
    ;'OAD3','OAD4','OAD5','OAD6','OAD7','OAD15','OAD16','OAD17','OAS1'};
Signal = {'dff','s465','s405'};
Neuromodulator = {'D','A'};

Trigger = 'StimulusOn';
Partition = 'Soundtype';
minpart = 2;
Color = {[0,1,0],[1,0,0],[0.5,0.5,0]};
Partitions = {'Constant go tone','Constant no-go tone','New tone'};
numParts = (1:3);
win = [-1,1];

for S = 1:length(Signal)
    for NM = 1:length(Neuromodulator)
        animalID_ = animalID(NM,:);
        Photometry_AVG_Run(animalID_,Trigger,[Signal{S},'_',Neuromodulator{N}],Partition,Partitions,minpart,numParts,Color,win);
    end
end

Trigger = 'DeliverFeedback';
Partition = 'Type4';
minpart = 2;
Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
Partitions = {'Constant go tone','Constant no-go tone','New go tone','New no-go tone'};
numParts = (1:4);
S = 1;
for NM = 1:length(Neuromodulator)
    animalID_ = animalID(NM,:);
    Photometry_AVG_Run(animalID_,Trigger,[Signal{S},'_',Neuromodulator{N}],Partition,Partitions,minpart,numParts,Color,win)
end

Trigger = 'DeliverFeedback';
Partition = 'CorrectRejection';
minpart = 1;
Color = {[0.75,0.5,0]};
Partitions = {'Correct Rejection'};
numParts = 1;
S = 1;
for NM = 1:length(Neuromodulator)
    animalID_ = animalID(NM,:);
    Photometry_AVG_Run(animalID_,Trigger,[Signal{S},'_',Neuromodulator{N}],Partition,Partitions,minpart,numParts,Color,win)
end

%% PSTH clustering
choosecb('ACh_DA_FinalCellbase'); 
loadcb;
popPSTH_Matrix ('HDB')
popPSTH_Matrix ('VTA')

%% CCG
path = [getpref('cellbase','datapath'),'CCG_final\'];
mkdir(path);

mkdir(path);
if ~exist([path,'ALL\CCG_matrices.mat'],'file')
    ccg(crosspairs,0.5,'issave',false,'resdir',[path,'All'],'segfilter','@stim_excl_dual','filterinput',{'light_activation_duration',[0 1],'margins',[0,0]})
end

if ~exist([path,'Cue\CCG_matrices.mat'],'file')
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Cue'],'segfilter','@cue_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'Hit'],'segfilter','@Hit_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
    ccg(crosspairs,0.25,'issave',false,'resdir',[path,'FA'],'segfilter','@FA_incl_nb','filterinput',{'margins',[0 0],'min_int',0},'minspikeno',5)
end

if ~exist([path,'ISI\CCG_matrices.mat'],'file')
    ccg(CELLIDLIST(~isnan(HDB_Clust)|~isnan(VTA_Clust)),0.5,'issave',false,'resdir',[path,'ISI'],...
        'segfilter','@stimfb_excl_dual','filterinput',{'light_activation_duration',[0 1],'feedback_duration',[-0.6 0.6],'margins',[0 0]})
end

%% Adding analysis to CEllBASE
addanalysis(@LRatio2,'property_names',{'ID_PC','Lr_PC'},'arglist',{'feature_names' {'WavePC1' 'Energy'}})
addanalysis(@cellid2vals,'property_names',{'RatId','DateNum','Tetrode','Unit'});
addanalysis(@nbisstim,'property_names',{'Hindex','D_KL'},'arglist',{'a'})
addanalysis(@nbisstim,'property_names',{'Hindex_DA','D_KL_DA'},'arglist',{'b'})
addanalysis(@spikeshapecorr,'property_names',{'R'},'arglist',{'a'})
addanalysis(@spikeshapecorr,'property_names',{'R_DA'},'arglist',{'b'})
addanalysis(@spikeshapeanalysis,'property_names',{'firing_rate','SpikeWidth','Peak_Integral_norm','PostValleyAmplNormToPeak','PeakToPostValleyTime'})

%% DREAD-FiberPhotometry PSTHs 

choosecb('Cellbase_psychometric')
loadcb
for i=1:length(CELLIDLIST)
    Partition = 'Type5';
    minpart = 3;
    Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
    Partitions ={'Constant go tone','Constant no-go tone','New go tone','New no-go tone'};%{'Constant go tone','Constant no-go tone','New tone'}
    numParts = [1,2,3,4];
    win = [-2,3];
    
    [FM_SO_D,F_SO_D,~,~] = Photometry_AVG(CELLIDLIST{i},win,'StimulusOn','dff_D',Partition,numParts,Color,Partitions,minpart);
    [FM_SO_A,F_SO_A,type_SO,validsessions] = Photometry_AVG(CELLIDLIST{i},win,'StimulusOn','dff_A',Partition,numParts,Color,Partitions,minpart);
    
    Partition = 'Type4e';
    Partitions ={'Constant hit','False alarm','New hit'};
    numParts = [1,2,3];
    minpart = 1;
    
    [FM_DF_D,F_DF_D,~,~] = Photometry_AVG(CELLIDLIST{i},win,'DeliverFeedback','dff_D',Partition,numParts,Color,Partitions,minpart,validsessions);
    [FM_DF_A,F_DF_A,type_DF,~] = Photometry_AVG(CELLIDLIST{i},win,'DeliverFeedback','dff_A',Partition,numParts,Color,Partitions,minpart,validsessions);
    for ses_inx = 1:length(validsessions)
        [Hit(ses_inx),FalseAlarm(ses_inx),RT(ses_inx)] = Punish(CELLIDLIST{i},validsessions(ses_inx ).name);
    end
    cellbasepath = [getpref('cellbase','datapath')];
    save([cellbasepath,'\',CELLIDLIST{i},'\',CELLIDLIST{i},'.mat'], 'F_SO_D','F_SO_A','F_DF_D','F_DF_A','FM_SO_D','FM_SO_A','FM_DF_D','FM_DF_A','type_SO','type_DF','RT','FalseAlarm','Hit','validsessions');
end



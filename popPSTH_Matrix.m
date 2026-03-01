function popPSTH_Matrix (AREA_)  %(cellids,varargin)

% Preprocesses data to identify groups of similary behaving cells using the
% first three principal components of the Z-scored curves 
% KB 2018/10 2019/07 2026/03

dbstop if error
% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    Area = getvalue('Area1');
    Mousenum=getvalue('RatID');
    Tetrodenum=getvalue('Tetrode');
    isArea=strcmp(Area,AREA_);
    
    isAnimal=(Mousenum~=0);
    ptinx = ID > 20 & Lratio < 0.15 & isArea & isAnimal ;   
    I = find(ptinx);
    IDS = CELLIDLIST(I);
    
% Find putative tagged cells
    Hindex = getvalue('Hindex');
    R = getvalue('R');
    Hindex_DA = getvalue('Hindex_DA');
    R_DA = getvalue('R_DA');
    ptinx = ID > 20 & Lratio < 0.15 & ((Hindex < 0.01 & R > 0.85) | (Hindex_DA < 0.01 & R_DA > 0.85) ) & isArea;  %% R 0.85 alá biztos nem mehet, H 0.01, ID-Lratio  
    I = find(ptinx);
    taggedIDS = CELLIDLIST(I);
    
skipinx=[];


default={...
    'window',               [-3 3];... % time window around the 1st trigger event: calculated  // % !!! do not set it too narrow: longer increased or decsreased block can pull away the mean value of Z-score !!!
    'pwindow'               [-0.5,1];...% time window around the 1st trigger: plotted
    'cwindow',              [0.0,0.6;0.4,0.6];... % time window around the 1st trigger: considerd for the clustering 
    'dt',                   0.005;... % time resolution (s)  
    'sigma',                0.02;...  % smoothing kernel
    'Trigger1Name',         'StimulusOn';... % 1st Trigger --> time(0:Trigger2)  //e.g. cue 
    'Trigger2Name'          'DeliverFeedback';... % 2nd Trigger --> time(Trigger2:end) //e.g. feedback
    'Partitions',           '#Hit';...  % trial partitions
    'ClusterPartitions'     [1,2];... % partison # used for clustering e.g. [1,2]; the first one wil be used for sorting
    'PlotPartitions'        [1,2];... % partition # to be plotted 
    'baselinelength'        1;... % length of the baseline period (s) from the begining of the 'window' //only for auROC
    'ROCWindow'             10;... % length of the moving window (in dt points), used for the auROC histogram comparisson
    'avgtimediff'           0.4;... % ~avg time difference between triggers
    'integrallength'        0.6;... % length of the window (from the first trigger event) used for sorting clusters and cells
    'maxtimediff'           1;...  % maximum time difference between triggers  
    'ClusterNum'            5;... % required number of clusters
    'PCA_dim'               [2,1];...
    'clustering_method'     'complete';...
    'clustering_metric'     'euclidean';...
    'LastEvents',           '';...
    'BurstPSTH'             'off';...
    'showtagged'            '*';... % '*':mark tagged cells with * / 'IDS': write tagged IDS  / 0: don't mark tagged neurons / 'showall' write all IDS
    'norm_method'           'Z-score';... % 'auROC' / 'Z-score'
    'isfilter'              0;... %1: exlude neurons from tetrodes without neurons from the tagged group
    'usedoubletrigger'      1;... % 1:yes; 0:no
    'savefig'               0;... %1:yes, 0:no 
    'savematrix'            0;... %1:yes, 0:no 
    'avgplotcolors'         ['g','r','y','b'];...
    'savepath'              'D:\CellBase_DATA\Population_maps\ALL_20220816\HDB\',... 
    };
%--------------------------------------------------------------------------
    % Calculate the Z-score and AUROC psth of the cells 
%--------------------------------------------------------------------------

%[g,error] = parse_args(default,varargin{:});
[g,error] = parse_args(default);
mkdir(g.savepath)

% defining time windows
g.baselinelength=g.baselinelength/g.dt;
margin=g.sigma*3;
time=g.window(1):g.dt:g.window(2);
time_m=g.window(1)-margin:g.dt:g.window(2)+margin; %plotted time window + margin for smoothing
time_plot=g.pwindow(1):g.dt:g.pwindow(2); % plotted time window
time_calc=g.window(1)-g.maxtimediff-margin:g.dt:g.window(2)+g.maxtimediff+margin; % extended time window for calculation
for partnum=1:length(g.ClusterPartitions)
ctime=g.cwindow(partnum,1):g.dt:g.cwindow(partnum,2); % time window considered for clustering
[~,c1]=(min(abs(time-ctime(1))));
[~,c2]=(min(abs(time-ctime(end))));
cinx{partnum,:}=c1:c2;
end

i=0;
for cellinx=1:length(IDS)
    
    cellid=IDS(cellinx);
    
    % load event files
    TE = loadcb(cellid,'TrialEvents');
    SP = loadcb(cellid,'EVENTSPIKES');
    
    % redifine partitions if needed; 
    delinx=~ismember(TE.(g.Partitions(2:end)),[g.ClusterPartitions,g.PlotPartitions]);
    TE.(g.Partitions(2:end))(delinx)=NaN;

    TE.Hit(TE.FalseAlarm==1)=2;
    TE.Hit(TE.CorrectRejection==1)=3;
    TE.Hit(TE.Miss==1)=4;
    
    % Partition trials, and skip cells without at least 1 trial from all required partitions

    
    try
        [COMPTRIALS2, ~] = partition_trials(TE,g.Partitions);
        COMPTRIALS = COMPTRIALS2(g.ClusterPartitions);
        if min(cellfun('length',COMPTRIALS))<5
            sprintf('%s skipped because of a missing partiton',cellid{1})
            skipinx=[skipinx cellinx];
            continue
        end
    catch
        sprintf('%s skipped because of a missing partiton',cellid{1})
        skipinx=[skipinx cellinx];
        continue
    end
    
    try
       COMPTRIALS = COMPTRIALS2(g.PlotPartitions);
    end
    
      
    % find valid trials in which Trigger2 (feedback) happened
    trigev1 = TE.(g.Trigger1Name);
    trigev2 = TE.(g.Trigger2Name);
    if ~iscell(trigev2)
        valid_trials = find(~isnan(trigev1) & ~isnan(trigev2));
    else
        valid_trials = find(cellfun(@(s)~isempty(s),trigev2));
    end
    
    % find Trigger2 events
    trigger_pos = findcellstr(SP.events(:,1),g.Trigger2Name);
    
    % TriggerName mismatch
    if trigger_pos == 0
        error('Trigger name not found');
    else
        TriggerEvent = SP.events{trigger_pos,2};
    end
    if ~isfield(TE,TriggerEvent)
        error('TriggerEvent mismatch: supply correct Events structure')
    end
    
    % Spike times
    alltrials = 1:size(SP.event_stimes{1},2);
    stimes  = SP.event_stimes{trigger_pos}(alltrials);
    if strcmp(g.BurstPSTH,'on'),
        stimes = detect_bursts(stimes);
    end
    
    % Event windows
    if ~iscellstr(g.LastEvents) && (strcmpi(g.LastEvents,'none') || isempty(g.LastEvents))
        ev_windows = SP.event_windows{trigger_pos};
    else
        ev_windows = get_last_evtime(TE,TriggerEvent,g.LastEvents);
    end
    
    %--------------------------------------------------------------------------
    % Make the main raster
    %--------------------------------------------------------------------------
    
    % Calculate binraster
    NUMtrials = length(TE.(g.Trigger1Name));
    if iscell(stimes{1})   % deal with lick-aligned raster
        stimes2 = [stimes{1:end}];
        binraster0 = stimes2binraster(stimes2,time_calc,g.dt);
        binraster = nan(NUMtrials,size(binraster0,2));
        %     binraster2 = nan(NUMtrials,size(binraster0,2));
        for k = 1:NUMtrials   % calculate sum of rows for each trial, which will be used for the PSTH
            sind = sum(cellfun(@length,stimes(1:k-1))) + 1;
            eind = sind + length(stimes{k}) - 1;
            %         disp([sind eind])
            binraster(k,:) = mean(binraster0(sind:eind,:),1);
            %         binraster2(k,:) = sum(stimes2binraster(stimes{k},time,g.dt),1);
        end
    else
        binraster = stimes2binraster(stimes,time_calc,g.dt);
    end
    
    % For variable windows, change padding to NaN to ensure correct averaging - BH
    if ~isempty(g.LastEvents)
        for iT = 1:NUMtrials    % loop through trials
            inx = time_calc > ev_windows(iT,2);
            binraster(iT,inx) = NaN;
        end
    end

    % Calculate double triggered rasters
    if g.usedoubletrigger==1
        trigtime=find(time_calc==0);
        EventTimes = trialevents2relativetime(TE,TriggerEvent,g.Trigger1Name); % time of trgi1 relative to trig2
        Braster=NaN(size(binraster,1),length(time_m)); % inicialize double triggered raster
        trigtimediff=floor((TE.(TriggerEvent)-TE.(g.Trigger1Name))/g.dt); % timedifference between triggers 
        trigtimediff(trigtimediff>g.avgtimediff/g.dt)=g.avgtimediff/g.dt; % cut trials with longer time diff than avg
        for j=valid_trials
            if abs(EventTimes(j))<g.maxtimediff  % check if the time diff is too long because of a recordign mistake and skip problematic trials
                B1=NaN(1,floor(g.avgtimediff/g.dt)); % inicialize raster between triggers
                % calculate raster between triggers; in case of shorter
                % then avg time diff between triggers bins reamin filled with nans 
                B1(1:trigtimediff(j))=binraster(j,trigtime+floor(EventTimes(j)/g.dt):trigtime+floor(EventTimes(j)/g.dt)+trigtimediff(j)-1); 
                startinx=trigtime+floor((EventTimes(j)/g.dt))+floor(time_m(1)/g.dt):trigtime+floor(EventTimes(j)/g.dt)-1; %time before trig1
                %fill rester [before trig1; between triggers; after trig2 ]
                Braster(j,:)=[binraster(j,startinx),B1,binraster(j,trigtime:trigtime+ceil(time_m(end)/g.dt-g.avgtimediff/g.dt))];
            else
                sprintf('Omitted trial in %s (time difference between the trigger events bigger than maxtimediff) \n', cellid{1})
            end
        end

    end
    
    
         % Calculate PSTH
    if g.usedoubletrigger==1     
    [psth, spsth, ~] = binraster2psth(Braster,g.dt,g.sigma,COMPTRIALS,valid_trials);
    psth=psth(:,find(abs((time_m-g.window(1)))<g.dt*0.99):find(abs((time_m-g.window(2)))<g.dt*0.99));
    spsth=spsth(:,find(abs((time_m-g.window(1)))<g.dt*0.99):find(abs((time_m-g.window(2)))<g.dt*0.99));
    else % temmporary solutiuon for single triggers
    
        g.Trigger1Name=g.Trigger2Name;
        [psth, spsth, ~] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);
        psth=psth(:,floor(g.maxtimediff/g.dt):end-floor(g.maxtimediff/g.dt));
        spsth=spsth(:,floor(g.maxtimediff/g.dt):end-floor(g.maxtimediff/g.dt));
        g.avgtimediff=0;
    end
    
    %--------------------------------------------------------------------------
    % Calculate the Z-score and the auROC
    %--------------------------------------------------------------------------
    test=sum(isnan(zscore(squeeze(spsth(:,:)),0,2))')>1;
    if sum(test(g.ClusterPartitions))>0
        sprintf('%s skipped because there was not enough data from a partiton',cellid{1})
        skipinx=[skipinx cellinx];
        continue
    end
    i=i+1;
    
    trigtime=find(time_m==0);
    for k=1:length(g.PlotPartitions)
        if size(spsth,1)<length(g.PlotPartitions) & ismember(k,setdiff(g.PlotPartitions,g.ClusterPartitions))
           continue 
        end
        % Z-score map
        Zpsth(i,k,:)=zscore(squeeze(spsth(k,:)),0,2); % Z-score spsth
        integralrespons(i,k)=sum(Zpsth(i,k,trigtime:trigtime+floor(g.integrallength/g.dt))); %integrated respons from trigger1 for sorting
       
        % auROC map
        maximum=floor(max(psth(k,:)))+1; %maximal length of the histogram
        b=histcounts(psth(k,1:g.baselinelength),'BinWidth',1,'BinLimits',[0 , maximum]); % computing baselinehistorgram
        baselinehist=b/g.baselinelength; % baseline histogram normalization
        z=0;
        for q=g.ROCWindow:g.ROCWindow:length(psth)
            z=z+1;
            a=histcounts(psth(k,q-(g.ROCWindow-1):q),'BinWidth',1,'BinLimits',[0 , maximum]); %computing moving window histogram
            ahist=a/g.ROCWindow;         %moving window histogram normalization
            for qq=1:1:maximum                  % moving the criteria
                p_a(qq)=sum(ahist(qq:end));     % probability that the activity during the window is bigger than the criteria
                p_b(qq)=sum(baselinehist(qq:end)); % probability that the activity during the baseline is bigger than the criteria
                
            end
            auROC(i,k,z)=trapz(p_b,p_a)*-1; % the area under the p_a-p_b curve (quantifies the degree of overlap between the two spike count distirbuitons)
            clear p_b
            clear p_a
        end
             
    end
    
end

IDS_orig=IDS;
IDS(skipinx)=[];

if g.savematrix==1
   save([g.savepath 'DATA.mat'],'Zpsth','g','IDS_orig','inx','IDS','skipinx')
end


        


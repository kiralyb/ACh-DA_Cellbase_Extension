function viewphotometry(animalID,sessionID,varargin)
%VIEWPHOTOMETRY   Colour map and PSTH of fiber photometry data.
%   VIEWPHOTOMETRY(ANIMALID,SESSIONID,'TRIGGERNAME',TRIGEVENT,'SORTEVENT',SEVENT,'EVENTTYPE
%   ',EVTYPE) pl

% Default arguments
default_args={...
    'window',               [-2 3];...
    'dt',                   0.01;...
    'sigma',                1;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'Signal'                'dff' %s465_D, s405, dff
    'TriggerEvent',         'DeliverFeedback';...
    'SortEvent',            'TrialStart';...
    'ShowEvents',           {{'StimulusOn'}};...
    'ShowEventsColors',     {{[0 0.8 0] [0.8 0.8 0] [0 0.8 0.8]}};...
    'Num2Plot',             'all';...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',      ':';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    'BurstPSTH'             'off';......
    'isephys'                0;....
    };
[g,error] = parse_args(default_args,varargin{:});

if (animalID(1)=='O' || animalID(1)=='D')
    %choosecb('Cellbase_psychometric')
    animalID_old =['VFPO'] ;
elseif (animalID(1)=='S')
    choosecb('Cellbase_sequential')
    animalID_old =['SEQ'] ;
elseif (animalID(1)=='B')
    choosecb('Cellbase_bandit')
    animalID_old =['VFP'] ;
end

% Input argument check
if nargin < 1
    error('Not enough input arguments. Please provide, at least, an animal and a session.');
end

if strcmp(g.Signal,'dff_1') 
    isnorm = 1;
elseif strcmp(g.Signal,'dff_2')
    isnorm = 1;
else
    isnorm = 0;
end

% Load data and FiberEvents
fullpath = getpref('cellbase','datapath');
path = [fullpath filesep animalID filesep sessionID filesep];

if matches(g.Signal,'ephys')    
    [signal, tss] = load_open_ephys_data([path filesep '100_CH1.continuous']);
    TE = load([path filesep 'TrialEvents.mat']);
    sr = 1000;
    signal = signal(1:30:end);
    tss = tss(1:30:end);

else
    DATA = load([path filesep 'proF.mat']);
    TE = load([path filesep 'FiberEvents.mat']);
    sr = DATA.sr;
    tss = DATA.tss-DATA.tss(1);
    signal = DATA.(g.Signal);

end
TE.Blocknum(TE.Hit~=1)=NaN;
% Extracting valid trials
[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
COMPTRIALS = COMPTRIALS(vinx);
TAGS = TAGS(vinx);
trigev = TE.(g.TriggerEvent);
if ~iscell(trigev)
    valid_trials = find(~isnan(trigev));
else
    valid_trials = find(cellfun(@(s)~isempty(s),trigev));
end

% Creating time vector
time = g.window(1):(1/sr):g.window(2);

% Sort trials
NUMevents = length(g.SortEvent);
if iscellstr(g.SortEvent)
    sort_var = nan(NUMevents,NUMtrials);
    for iS = 1:NUMevents
        sort_var(iS,:) = TE.(g.SortEvent{iS}) - TE.(g.TriggerEvent);
    end
    sort_var = min(sort_var);
elseif ~isempty(g.SortEvent)
    if ~iscell(TE.(g.TriggerEvent))
        sort_var = TE.(g.SortEvent) - TE.(g.TriggerEvent);
    else
        gte = nan(1,NUMtrials);
        inx = ~cellfun(@isempty,TE.(g.TriggerEvent));
        gte(inx) = cell2mat(cellfun(@(s)s(1),TE.(g.TriggerEvent)(inx),'UniformOutput',false));
        sort_var = TE.(g.SortEvent) - gte;
    end
else
    sort_var = NaN;
end

[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_Balint,TAGS);


% Peri-event matrix for the given contingencies
for iPAR = 1:length(TAGS)
    valT = intersect(valid_trials,COMPTRIALS{1,iPAR});
    if isempty(valT)
        continue
    end
    sortT = TE.(g.TriggerEvent)(valT);
    vdisc = [];
    if strcmp(g.TriggerEvent,'TrialStart')
        [~,vdisc] = min(abs((sortT) - (tss))); 
    else
        startT = TE.TrialStart(valT);
        %vdisc = round((startT + sortT) * DATA.sr);
        for i= 1:length(startT)
        [~,vdisc(i)] = min(abs((startT(i) + sortT(i)) - (tss)));
        end
    end
    

    
    [fibmean,fibcv,spmat] = perievent(vdisc,signal,sr,g.window,isnorm,g.dt,g.sigma); % peri-event
    if size(spmat,1) < length(valT)
        valT(size(spmat,1)+1:end)=[];
    end
    %zspmat=zscore(spmat,0,2);
    
    nT = [1:1:length(valT)];
      
    subplot(length(TAGS)+1,1,iPAR);
    [~,inx]=sort(sort_var(valT));
    N=imagesc(time(1:100:end),nT,spmat(inx,1:100:end));
    %imagesc(time,nT,spmat);
    %colormap(gray);
    if iPAR==1
    title([g.Signal ' aligned to ' g.TriggerEvent],'fontsize',12);
    c_limit = quantile(N.CData(:),[0.99,0.01]);
    end
    try
    caxis([c_limit(2) c_limit(1)]);
    end
    ylabel('Trial Number','fontsize',12);
    if iPAR == length(TAGS)+1
        xlabel(['Seconds from ' g.TriggerEvent ' onset'],'fontsize',12);
    end
%     try
%     blocknum = TE.TrialsinBlock(TE.(g.Partitions(2:end))==iPAR);
%     line(g.window,[find(diff(blocknum)<0);find(diff(blocknum)<0)],'Color','c','LineWidth',1)
%     end
    subplot(length(TAGS)+1,1,length(TAGS)+1);
    hold on
    errorshade(time(1:100:end),fibmean(:,1:100:end),fibcv(:,1:100:end),'LineColor',mycolors{iPAR},'ShadeColor',mycolors{iPAR});
    ylabel(g.Signal,'fontsize',12);
    xlabel(['Seconds from ' g.TriggerEvent ' onset'],'fontsize',12);
    
end
legend(mylabels)

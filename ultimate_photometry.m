function [fibmean,spmat,TAGS] = ultimate_photometry(animalID,sessionID,varargin)
%VIEWPHOTOMETRY   Colour map and PSTH of fiber photometry data.
%   VIEWPHOTOMETRY(ANIMALID,SESSIONID,'TRIGGERNAME',TRIGEVENT,'SORTEVENT',SEVENT,'EVENTTYPE
%   ',EVTYPE) pl

% Default arguments
default_args={...
    'window',               [-2 3];...
    'dt',                   0.01;...
    'sigma',                0.05;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'Signal'                'dff_D' %s465, s405, dff
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
    'BurstPSTH'             'off';...
    'COMPTRIALS'            '[]';...
    'FilterPartition'       '';...
    'Filter'                [];...
    };
[g,error] = parse_args(default_args,varargin{:});

% if (animalID(1)=='O')
%     %choosecb('Cellbase_psychometric')
%     animalID_old =['VFPO'] ;
% elseif (animalID(1)=='S')
%     choosecb('Cellbase_sequential')
%     animalID_old =['SEQ'] ;
% elseif (animalID(1)=='B')
%     choosecb('Cellbase_bandit')
%     animalID_old =['VFP'] ;
% end

% Input argument check
if nargin < 1
    error('Not enough input arguments. Please provide, at least, an animal and a session.');
end

% Load data and FiberEvents
fullpath = getpref('cellbase','datapath');
path = [fullpath filesep animalID filesep sessionID filesep];
TE = load([path filesep 'FiberEvents.mat']);
varsInFile = who('-file', fullfile(path, 'proF.mat'));
if ~ismember(g.Signal, varsInFile)
    g.Signal = g.Signal(1:find(g.Signal == '_', 1, 'last')-1);
end        

DATA = load([path filesep 'proF.mat'],g.Signal,'tss','sr');

%TE.Blocknum(TE.Hit~=1)=NaN;

% Extracting valid trialsN
if isempty(g.COMPTRIALS)
    TAGS = 1;
    COMPTRIALS = g.COMPTRIALS;
    valid_trials = COMPTRIALS;
else
    [COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
    vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
    COMPTRIALS = COMPTRIALS(vinx);
    TAGS = TAGS(vinx);
    trigev = TE.(g.TriggerEvent);
    if ~isempty(g.Filter)
        filtered_trials = find(TE.(g.FilterPartition(2:end)) == g.Filter);
        for iPAR = length(TAGS):-1:1
            COMPTRIALS{1,iPAR} = intersect(filtered_trials,COMPTRIALS{1,iPAR});
            if isempty(COMPTRIALS{1,iPAR})
                COMPTRIALS(iPAR) = [];
                TAGS(iPAR) = [];
            end
        end
    end
    if ~iscell(trigev) 
        valid_trials = find(~isnan(trigev));
    else
        valid_trials = find(cellfun(@(s)~isempty(s),trigev));
    end
end

% Peri-event matrix for the given contingencies
for iPAR = 1:length(TAGS)
    valT = intersect(valid_trials,COMPTRIALS{1,iPAR});
    sortT = TE.(g.TriggerEvent)(valT);
    vdisc = [];
    if strcmp(g.TriggerEvent,'TrialStart')
        [~,vdisc] = min(abs((sortT) - (DATA.tss-DATA.tss(1))));
    else
        startT = TE.TrialStart(valT);
        for i= 1:length(startT)
            [~,vdisc(i)] = min(abs((startT(i) + sortT(i)) - (DATA.tss-DATA.tss(1))));
        end
    end
    [fibmean(iPAR,:),fibcv(iPAR,:),spmat{iPAR}] = perievent(vdisc,DATA.(g.Signal),DATA.sr,g.window,0,1/DATA.sr,g.sigma); % peri-event
end


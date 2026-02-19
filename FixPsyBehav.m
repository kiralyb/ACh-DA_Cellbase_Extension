function [hit,FA,RT_m,state,L,nc,LN] = FixPsyBehav(animalID,sessionID, filter) % Sums up the water obtained
dbstop if error
fullpth = [getpref('cellbase','datapath') filesep animalID filesep sessionID filesep];
filename2 = [animalID,sessionID,'.mat'];
load(fullfile(fullpth,filename2))

khz4_index = find(SessionData.TrialTypes(1:end-1)==1);
khz10_index = find(SessionData.TrialTypes(1:end-1)==2);
khzN_index = find(SessionData.TrialTypes(1:end-1)>2);

if filter ~= 0
    M = vertcat(SessionData.TrialSettings.TrialType);
    if size(M,2)==3
        index = find(M(:,3)==filter);
        %if filter == 1
        %    index = find(M(:,3)==1 & SessionData.ThrdSound > 7)
        %else
        %    index = find(M(:,3)==2 & SessionData.ThrdSound < 7)
        %end
    else
        index = [];
    end
else
    index = 1:SessionData.nTrials-1;
end

L = length(khz4_index);
LN = length(khzN_index);

if isempty(khz10_index)
    state = 1;
elseif isempty(khzN_index)
    state = 2;
else
    state = 3;
end
if SessionData.TrialSettings(SessionData.nTrials).GUI.Phase==5 && isfield(SessionData,'thrdSound')
    nc = sum(diff(SessionData.thrdSound)~=0);
else
    nc = nan;
end
IND = intersect(khz4_index,index);
hit = mean(SessionData.Outcomes(IND) == 1);
FA = mean(SessionData.Outcomes(intersect(khz10_index,index)) == 0);
for i = 1:length(IND)
    RT(i)= SessionData.RawEvents.Trial{1, IND(i)}.States.StartStimulus(2)-SessionData.RawEvents.Trial{1, IND(i)}.States.StartStimulus(1);
end
if exist('RT')
    RT_m = mean(RT);
else
    RT_m = nan;
end
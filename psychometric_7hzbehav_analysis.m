function [Hitrate,Hitrate4,FalseAlarmrate10] = psychometric_7hzbehav_analysis(animalID,sessionID)

fullpath = getpref('cellbase','datapath');
path = [fullpath '\' animalID filesep sessionID filesep];
try
    TE = load([path '\' 'TrialEvents.mat']);
catch
    TE = load([path '\' 'TE.mat']);
end

Hitrate = nansum(TE.Hit(TE.Type3 == 3))/ nansum((TE.Type3 == 3));
Hitrate4 = nansum(TE.Hit(TE.Type3 == 1))/ nansum((TE.Type3 == 1));
FalseAlarmrate10 = nansum(TE.FalseAlarm(TE.Type3 == 2))/ nansum((TE.Type3 == 2));
%Hitrate = nanmean(TE.StimulusDuration(TE.Type == 1));


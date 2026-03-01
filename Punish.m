function [hit,FalseAlarm,RT_m,state] = Punish(animalID,sessionID,issave, filter) % Sums up the water obtained
if nargin < 4
    filter = 0;
end

if issave
fullpth2 = [getpref('cellbase','datapath') filesep animalID filesep 'behav\'];
mkdir(fullpth2);
end
fullpth = [getpref('cellbase','datapath') filesep animalID filesep sessionID filesep];
filename2 = [animalID,sessionID,'.mat'];
load(fullfile(fullpth,filename2))
wind=20;
b=1;
for a=1:wind:floor((SessionData.nTrials-1)/wind)*wind
    y=1;
    z=1;
    w=1;
    Outcomes1=NaN;
    Outcomes2=NaN;
    Outcomes3=NaN;
for x = a:a+20
    if (SessionData.TrialTypes(x)==2)
    if ~isnan(SessionData.RawEvents.Trial{x}.States.Punish(1))
        Outcomes2(y) = 1;
    else
        Outcomes2(y)=0;
        
    end
    y=y+1;
    end
    if (SessionData.TrialTypes(x)==1)
    if ~isnan(SessionData.RawEvents.Trial{x}.States.Reward(1))  
        Outcomes1(z) = 1;
    else
        Outcomes1(z)=0;
        
    end
    z=z+1;
    end
        if (SessionData.TrialTypes(x)>=3)
    if ~isnan(SessionData.RawEvents.Trial{x}.States.Reward(1)) || ~isnan(SessionData.RawEvents.Trial{x}.States.Punish(1))
        Outcomes3(w) = 1;
    else
        Outcomes3(w)=0;   
    end
    w=w+1;
    end
end
Outcomes11(b)=sum(Outcomes1)/length(Outcomes1);
Outcomes22(b)=sum(Outcomes2)/length(Outcomes2);
Outcomes33(b)=sum(Outcomes3)/length(Outcomes3);

if isnan(Outcomes22)
    state = 1;
elseif isnan(Outcomes33)
    state = 2;
else
    state = 3;
end
if SessionData.TrialSettings(SessionData.nTrials).GUI.Phase==0
    state = 0;
end


b=b+1;
end

figure;
subplot(3,1,1)
%A=[Outcomes11];
A=[Outcomes11;Outcomes22;Outcomes33];
B=bar(A');
try
B(2).FaceColor='red';
B(1).FaceColor='green';
B(3).FaceColor='y';
end
xlabel('Trial (window=20)')
ylabel('Hitratio')
if ~isfield(SessionData,'thrdSound')
    SessionData.thrdSound = zeros(1,SessionData.nTrials);
end
border=find(diff(SessionData.thrdSound)~=0)+1;
B=(border/wind);
hold on
line([B;B],[zeros(1,length(B));ones(1,length(B))],'LineWidth',3,'Color','k')

subplot(3,1,2)
[hit,Miss,FalseAlarm,CorrectReject]=Punish2(SessionData,filter);
khz4_index = find(SessionData.TrialTypes==1);
for i = 1:length(khz4_index)
RT(i)= SessionData.RawEvents.Trial{1, khz4_index(i)}.States.StartStimulus(2)-SessionData.RawEvents.Trial{1, khz4_index(i)}.States.StartStimulus(1);
end
RT_m = mean(RT);

bar([hit,Miss;FalseAlarm,CorrectReject],'stacked')
xticklabels({'4 Hz','10Hz'})
ylabel('performance')

subplot(3,1,3)
hold on
lab = {};
SessionData.TrialTypes(end)=[];
for t_index= 3:max(SessionData.TrialTypes)
hit_3 = sum(SessionData.Outcomes(SessionData.TrialTypes==t_index)==1)/sum(SessionData.TrialTypes==t_index);
miss_3 = sum(SessionData.Outcomes(SessionData.TrialTypes==t_index)==2)/sum(SessionData.TrialTypes==t_index);
falsealarm_3 = sum(SessionData.Outcomes(SessionData.TrialTypes==t_index)==-1)/sum(SessionData.TrialTypes==t_index);
correctr_3 = sum(SessionData.Outcomes(SessionData.TrialTypes==t_index)==0)/sum(SessionData.TrialTypes==t_index);
set(gca,'ColorOrderIndex',1)
bar([t_index*2-1,t_index*2],[hit_3,miss_3;falsealarm_3,correctr_3],'stacked')
xticks([3*2-1:max(SessionData.TrialTypes)*2])
lab = [lab,{'go','nogo'}];
xticklabels(lab)
end
ylabel('performance')
T = unique(SessionData.thrdSound, 'stable');
title([mat2str(T)])
if issave
saveas(gcf,[fullpth,filename2,'behav.png']);
saveas(gcf,[fullpth,filename2,'behav.fig']);
saveas(gcf,[fullpth2,filename2,'behav.png']);
saveas(gcf,[fullpth2,filename2,'behav.fig']);
end
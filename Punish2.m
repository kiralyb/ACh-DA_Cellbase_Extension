function [hit,Miss,FalseAlarm,CorretReject]=Punish2(Data,filter) % Sums up the water obtained
hit=0;
slowhit=0;
FalseAlarm=0;
CorretReject=0;
Miss=0;
if filter ~= 0
    M = vertcat(Data.TrialSettings.TrialType);
    if size(M,2)==3
        index = find(M(:,3)==filter)
    else
        index = [];
    end
else
    index = 1:Data.nTrials
end
    
if isempty(index)
    hit=nan;
    slowhit = nan;
    FalseAlarm=nan;
    CorretReject=nan;
    Miss=nan;
end
    
for i=1:1:Data.nTrials
   
    if (Data.TrialTypes(i)==2)
        if ~isnan(Data.RawEvents.Trial{i}.States.Punish(1))
            FalseAlarm=FalseAlarm+1;
        else
            CorretReject=CorretReject+1;
        end
    end
    
    if (Data.TrialTypes(i)==1)
        if ~isnan(Data.RawEvents.Trial{i}.States.Reward(1))
            
            if isnan(Data.RawEvents.Trial{i}.States.Delay(1))
            hit=hit+1;
            else
            slowhit=slowhit+1;
            end
        else
           Miss=Miss+1; 
        end
    end
    
end

hit=hit/sum(Data.TrialTypes==1)
slowhit=slowhit/sum(Data.TrialTypes==1)
Miss=Miss/sum(Data.TrialTypes==1)
FalseAlarm=FalseAlarm/sum(Data.TrialTypes==2)
CorretReject=CorretReject/sum(Data.TrialTypes==2)

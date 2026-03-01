function DREADD_prerunner(animalID)
dbstop if error
choosecb('Photobase_ACh_DA')
Partition = 'Type5';
minpart = 3;
Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
Partitions ={'Constant go tone','Constant no-go tone','New go tone','New no-go tone'};%{'Constant go tone','Constant no-go tone','New tone'}
numParts = [1,2,3,4];
win = [-2,3];

[FM_SO_D,F_SO_D,type_SO,validsessions] = Photometry_AVG(animalID,win,'StimulusOn','dff_D',Partition,numParts,Color,Partitions,minpart);
[FM_SO_A,F_SO_A,type_SO,validsessions] = Photometry_AVG(animalID,win,'StimulusOn','dff_A',Partition,numParts,Color,Partitions,minpart);
Partition = 'Type4e'
Partitions ={'Constant hit','False alarm','New hit'}
numParts = [1,2,3];
minpart = 1;
[FM_DF_D,F_DF_D,type_DF,~] = Photometry_AVG(animalID,win,'DeliverFeedback','dff_D',Partition,numParts,Color,Partitions,minpart,validsessions);
[FM_DF_A,F_DF_A,type_DF,~] = Photometry_AVG(animalID,win,'DeliverFeedback','dff_A',Partition,numParts,Color,Partitions,minpart,validsessions);
for ses_inx = 1:length(validsessions)
    [Hit(ses_inx),FalseAlarm(ses_inx),RT(ses_inx)] = Punish(animalID,validsessions(ses_inx ).name); % Sums up the water obtained
end
cellbasepath = [getpref('cellbase','datapath')];
save([cellbasepath,'\',animalID,'\',animalID,'.mat'], 'F_SO_D','F_SO_A','F_DF_D','F_DF_A','FM_SO_D','FM_SO_A','FM_DF_D','FM_DF_A','type_SO','type_DF','RT','FalseAlarm','Hit','validsessions')

%%
win = [-2,3];
Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
sr = 12048;
time = linspace(win(1),win(2),sr*(win(2)-win(1))+1)
figure
linestyles={'-',':'}
for parts=1:4
   
subplot(4,4,(parts-1)*4+1)
hold on
for isdread = 0:1
errorshade(time,nanmean(F_SO_A{parts}(type_SO{parts}==isdread,:)),std(F_SO_A{parts}(type_SO{parts}==isdread,:))./sqrt(sum(type_SO{parts}==isdread)),'LineStyle',linestyles(isdread+1),'LineColor',Color(parts))
end
subplot(4,4,(parts-1)*4+2)
hold on
for isdread = 0:1
errorshade(time,nanmean(F_SO_D{parts}(type_SO{parts}==isdread,:)),std(F_SO_D{parts}(type_SO{parts}==isdread,:))./sqrt(sum(type_SO{parts}==isdread)),'LineStyle',linestyles(isdread+1),'LineColor',Color(parts))
end

if parts <4 
subplot(4,4,(parts-1)*4+3)
hold on
for isdread = 0:1
errorshade(time,nanmean(F_DF_A{parts}(type_DF{parts}==isdread,:)),std(F_DF_A{parts}(type_DF{parts}==isdread,:))./sqrt(sum(type_DF{parts}==isdread)),'LineStyle',linestyles(isdread+1),'LineColor',Color(parts))
end
subplot(4,4,(parts-1)*4+4)
hold on
for isdread = 0:1
errorshade(time,nanmean(F_DF_D{parts}(type_DF{parts}==isdread,:)),std(F_DF_D{parts}(type_DF{parts}==isdread,:))./sqrt(sum(type_DF{parts}==isdread)),'LineStyle',linestyles(isdread+1),'LineColor',Color(parts))
end
end
end


figure
for parts=1:4
    if ~mod(parts,2)
        isinv = -1;
    else
        isinv = 1;
    end
    
   
subplot(4,4,(parts-1)*4+1)
hold on
bar([Photoresponse2(mean(F_SO_A{parts}(type_SO{parts}==0,24097:30121))),Photoresponse2(mean(F_SO_A{parts}(type_SO{parts}==1,24097:30121)))])
subplot(4,4,(parts-1)*4+2)
hold on
bar([Photoresponse2(mean(F_SO_D{parts}(type_SO{parts}==0,24097:30121))),Photoresponse2(mean(F_SO_D{parts}(type_SO{parts}==1,24097:30121)))])
if parts <4 

subplot(4,4,(parts-1)*4+3)
hold on
bar([Photoresponse2(mean(F_DF_A{parts}(type_DF{parts}==0,24097:36121))),Photoresponse2(mean(F_DF_A{parts}(type_DF{parts}==1,24097:36121)))])
subplot(4,4,(parts-1)*4+4)
hold on
bar([Photoresponse2(mean(F_DF_D{parts}(type_DF{parts}==0,24097:36121))*isinv)*isinv,Photoresponse2(mean(F_DF_D{parts}(type_DF{parts}==1,24097:36121))*isinv)*isinv])
end
end
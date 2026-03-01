function Photometry_AVG_Run(animalID,Trigger,Signal,Partition,Partitions,minpart,numParts,Color,win)

sr = 12048;
for i=1:length(animalID)
    FM(i,:,:) = Photometry_AVG(animalID{i},win,Trigger,Signal,Partition,numParts,Color,Partitions,minpart);
end
cellbasepath = [getpref('cellbase','datapath')];
save([cellbasepath,'\',Partition,'_',Trigger,'_',Signal,'FM.mat'],'FM','animalID')
function performance_plotter(animalID,sessionID)

fullpth = [getpref('cellbase','datapath') filesep animalID filesep sessionID filesep];
filename2 = [animalID,sessionID,'.mat'];
load(fullfile(fullpth,filename2))




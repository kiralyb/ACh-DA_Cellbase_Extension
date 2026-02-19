function [DIFF] = psychometric_DREADD_COMP(animalID,starti)
R = cell(length(animalID),5);
P = cell(length(animalID),5);
R_c = cell(length(animalID),5);
P_c = cell(length(animalID),5);
threshold = 200;

for j = 1:length(animalID)
cellbasepath = [getpref('cellbase','datapath')];
a = dir([cellbasepath '\' animalID{j}]);
a = a(3:end);
a = a(vertcat(a.isdir));
%a = a(3:end);     
for i = starti(j): length(a)
    sessionID = a(i).name;
    if sessionID(end)=='d' || sessionID(end)=='e' || sessionID(end)=='f'
        try
            [R(j,:), P(j,:)] = psychometric_behav_analysis(animalID{j},sessionID, R(j,:), P(j,:),threshold,1);
        catch
            animalID{j}
            sessionID
        end
    else
        try
            [R_c(j,:), P_c(j,:)] = psychometric_behav_analysis(animalID{j},sessionID, R_c(j,:), P_c(j,:),threshold,1);
        catch
            animalID{j}
            sessionID
            
        end
    end
end
end

figure
hold on 
for i = 1:5
MR(i) = mean(cell2mat(R(:,i)'));
MP(i) = mean(cell2mat(P(:,i)'));
MR_c(i) = mean(cell2mat(R_c(:,i)'));
MP_c(i) = mean(cell2mat(P_c(:,i)'));

SR(i) = std(cell2mat(R(:,i)'))./sqrt(length(cell2mat(R(:,i)')));
SP(i) = std(cell2mat(P(:,i)'))./sqrt(length(cell2mat(P(:,i)')));
SR_c(i) = std(cell2mat(R_c(:,i)'))./sqrt(length(cell2mat(R_c(:,i)')));
SP_c(i) = std(cell2mat(P_c(:,i)'))./sqrt(length(cell2mat(P_c(:,i)')));

end

errorshade(5:9,MR_c,SR_c,'LineColor','g','ShadeColor','g','LineWidth',5)
errorshade(5:9,MR,SR,'LineColor',[0.66,0.66,0],'ShadeColor',[0.66,0.66,0],'LineWidth',5)
errorshade(5:9,MP_c,SP_c,'LineColor','r','ShadeColor','r','LineWidth',5)
errorshade(5:9,MP,SP,'LineColor',[1,0.8,0],'ShadeColor',[1,0.8,0],'LineWidth',5)

setmyplot_balazs;
xlim([4.5,9.5])
xlabel('New cue frequency (kHz)')
ylabel('Required number of trials')
legend({'Control','DREADD','Control','DREADD'})

RDIFF = cellfun(@mean, R) - cellfun(@mean, R_c);
PDIFF = cellfun(@mean, P) - cellfun(@mean, P_c);
DIFF = [RDIFF(:,2:4);PDIFF(:,2:4)];
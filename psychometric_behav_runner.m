%%
threshold = 200;
animalID = listtag('animal'); animalID([3])=[];
%animalID = {'VHD2','VHD4','VHD6','VVH1','VVH2','VVH4','VVH6','VVH7','VVH8','VVH9','VVH11','WHV1','WHV2'};
%animalID = {'VTA1','VHD1','VHD2','VHD4','VHD6','VVH1','VVH2','VVH4','VVH6','VVH7','VVH8','VVH9','VVH11','WHV1','WHV2'};
R = cell(length(animalID),5);
P = cell(length(animalID),5);
for j = 1:length(animalID)
cellbasepath = [getpref('cellbase','datapath')];
a = dir([cellbasepath '\' animalID{j}]);
a = a(4:end);
a = a(vertcat(a.isdir));
%a = a(3:end);    
for i = 1: length(a)
    sessionID = a(i).name;
    try
        [R(j,:), P(j,:)] = psychometric_behav_analysis(animalID{j},sessionID, R(j,:), P(j,:),threshold);
    end
end
end

%%
for q = 1:length(R)
    inx = find(R{q,3} ~= threshold,1);
    R{q,3} = R{q,3}(inx:end);
    %P{q,3} = P{q,3}(inx:end);
end

%%
figure

plot(cellfun(@mean, R)','-+','Color',[0,1,0,0.2])
plot(cellfun(@mean, P)','-+','Color',[1,0,0,0.2])

plot(nanmean(cellfun(@mean, P))','r','LineWidth',5)
plot(nanmean(cellfun(@mean, R))','g','LineWidth',5)

%%
figure
hold on 
for i = 1:5
MR(i) = mean(cell2mat(R(:,i)'));
MP(i) = mean(cell2mat(P(:,i)'));
%bar(i-0.2,mean(cell2mat(R(:,i)')),0.35,'w','EdgeColor','g')
%bar(i+0.2,mean(cell2mat(P(:,i)')),0.35,'w','EdgeColor','r')
%errorbar(i-0.2,mean(cell2mat(R(:,i)')),0,std(cell2mat(R(:,i)')),'k')
%errorbar(i+0.2,mean(cell2mat(P(:,i)')),0,std(cell2mat(P(:,i)')),'k')

%scatter(rand(1,length(cell2mat(R(:,i)')))/5+i-0.30,cell2mat(R(:,i)'),'g')    
%scatter(rand(1,length(cell2mat(P(:,i)')))/5+i+0.10,cell2mat(P(:,i)'),'r')

%scatter(rand(1,length(cell2mat(R(:,i)')))/3-1/6+i+4,cell2mat(R(:,i)'),'g')    
%scatter(rand(1,length(cell2mat(P(:,i)')))/3-1/6+i+4,cell2mat(P(:,i)'),'r')
end

hold on
plot(5:9,cellfun(@mean, R)','Color',[0,1,0,0.2])
%plot(5:9,flipud(cellfun(@mean, P)'),'Color',[1,0,0,0.2])
plot(5:9,cellfun(@mean, P)','Color',[1,0,0,0.2])


plot(5:9,fliplr(MP),'r','LineWidth',5)
plot(5:9,MR,'g','LineWidth',5)

setmyplot_balazs;
xlim([4.5,9.5])
xlabel('New cue frequency (kHz)')
ylabel('Required number of trials')

%%

figure
hold on 
for i = 1:5
MR(i) = mean(cell2mat(R(:,i)'));
MP(i) = mean(cell2mat(P(:,i)'));
end

plot(5:9,MR,'g','LineWidth',5)
plot(5:9,MP,'r','LineWidth',5)


hold on

plot(5:9,fillmissing(cellfun(@mean, R)','linear','EndValues','none'),'g')
plot(5:9,fillmissing(cellfun(@mean, P)','linear','EndValues','none'),'r')



setmyplot_balazs;
xlim([4.5,9.5])
xlabel('New cue frequency (kHz)')
ylabel('Required number of trials')
legend({'Control','DREADD','Control','DREADD'})

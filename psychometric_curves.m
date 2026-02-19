function psychometric_curves(animalID)

% initialize
threshold = 200;
R = cell(length(animalID),5);
P = cell(length(animalID),5);

%extract data
for j = 1:length(animalID)
    cellbasepath = [getpref('cellbase','datapath')];
    a = dir([cellbasepath '\' animalID{j}]);
    a = a(3:end);
    a = a(vertcat(a.isdir));
    for i = 1: length(a)
        try
            [R(j,:), P(j,:)] = psychometric_behav_analysis(animalID{j},a(i).name, R(j,:), P(j,:),threshold);
        end
    end
end

% exclude first new tone sessions
for q = 1:length(R)
    inx = find(R{q,3} ~= threshold,1);
    R{q,3} = R{q,3}(inx:end);
end


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
%legend({'Control','DREADD','Control','DREADD'})

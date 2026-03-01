function [FM,F,type,a] = Photometry_AVG(animalID,win,Trigger,Signal,Partition,numParts,Color,Partitions,minpart,FilterPart,Filter,a)
%animalID = 'BAA25'%'BAA29';
%Partitions = {'Hit','FalseAlarm','CorrectRejection','Miss'};
%Partition = 'Lucky'%'Lucky'; % Type Type3 Type4 
%Trigger = 'StimulusOn'; % DeliverFeedback, StimulusOn, PokeIn LedOn   
%Signal = 'dff_A_PrL'%'dff_A_PrL';
isnorm = 0;
%Color = {'g','r','y','m'};
%Color = {[0,1,0],[1,0,0],[0,0.5,0],[0.5,0,0]};
%choosecb('Cellbase_psychometric')


%win = [-2,3];
cellbasepath = [getpref('cellbase','datapath')];
if nargin < 12
    a = dir([cellbasepath '\' animalID]);
    a = a(3:end);
    a = a(vertcat(a.isdir));
end

if nargin < 10
    FilterPart = '#TrialType';
    Filter = [];% [] if not
end    
%a = a(3:end);                                                                                                                                           
%clear F_4kHz F_10kHz
% F_4kHz = nan(length(a),60241);
% F_10kHz = nan(length(a),60241);
F = cell(1,length(numParts));
type = cell(1,length(numParts));
isareturned = nan(length(a),1);
for i = 1:length(a)
    %for j = 1:length(Partitions) 
        %Partition = Partitions{j};
        try
         %TE = solo2trialevents_psychometric_kb([a(i).folder,'\',a(i).name,'\SEQ',animalID(4:end),a(i).name],1);
         %TE = solo2trialevents_psychometric_kb([a(i).folder,'\',a(i).name,'\',animalID,a(i).name],1);
         %FE = load([a(i).folder,'\',a(i).name,'\FiberEvents.mat']);
         %FE.Soundtype = TE.Soundtype;
         %FE.Type4 = TE.Type4;
         %FE.Type4e = TE.Type4e
         %FE.Type5 = TE.Type5;
%         FE.TrialType = TE.TrialType;
%         FE.Type_R = TE.Type_R;
%         FE.DeliverAFeedback = TE.DeliverAFeedback;
%           FE.GTHit = TE.GTHit;
%           FE.PokeIn = TE.PokeIn;
         %save([a(i).folder,'\',a(i).name,'\FiberEvents.mat'],'-struct','FE');
         [fibmean,~,Tags] = ultimate_photometry(animalID,a(i).name,'window',win,'Signal',Signal,'TriggerEvent',Trigger,'Partitions',['#',Partition],'FilterPartition',FilterPart,'Filter',Filter);
        if length(Tags) < minpart
           isareturned(i) = 0;
           continue 
        end
        for j=1:length(numParts)
            index = find(strcmp(Tags, [Partition,'=',num2str(numParts(j))]));
            %index = find(strcmp(Tags, [Partition,'=1']));
            if ~isempty(index)
                if isnorm
                    F{j}(size(F{j},1)+1,:) = fibmean(index,:)/max(fibmean(1,:));
                else
                    F{j}(size(F{j},1)+1,:) = fibmean(index,:);
                end
            else
                F{j}(size(F{j},1)+1,:) = nan(1,length(fibmean));
            end
                if a(i).name(end)=='d' | a(i).name(end)=='e' | a(i).name(end)=='f'
                    type{j}(size(F{j},1)) = 1;
                else
                    type{j}(size(F{j},1)) = 0;
                end
            
        end
        isareturned(i) = 1 ;
        catch
        isareturned(i) = 0  ; 
      end
   
 end
a = a(logical(isareturned)); 
figure
%Color = [0,0,1];
for k = 1:length(numParts)
    errorshade(linspace(win(1),win(2),length(F{k})),nanmean(F{k}),nanstd(F{k})/sqrt(size(F{k},1)),'LineColor',Color{k},'ShadeColor',Color{k})
    %errorshade(linspace(-2,3,length(F{k})),nanmean(F{k}),nanstd(F{k})/sqrt(size(F{k},1)),'LineColor',Color*(k-4)/6,'ShadeColor',Color)
    FM(k,:) = nanmean(F{k});
    hold on
end
save([cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'F.mat'], 'F','FM')
setmyplot_balazs
xlabel(['Time from ', Trigger, ' (s)'])
ylabel('Average normalized response')
title([Signal])
if strcmp(Partition,'Type') || strcmp(Partition,'Type3') || strcmp(Partition,'Type4') || strcmp(Partition,'Type5')
legend('Constant go tone','Constant no-go tone','New go tone','New no-go tone')
elseif strcmp(Partition,'Lucky')
legend('Likely Reward','Lucky Reward','Likely Omission','Unlucky Omission')    
else
legend(Partitions)    
end
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'.png'])
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'.fig'])
try
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{1},1),F{1})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_1.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{2},1),F{2})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_2.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{3},1),F{3})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_3.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{4},1),F{4})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_4.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{5},1),F{5})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_5.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{6},1),F{6})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_6.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{1},1),F{1}-F{4})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_7.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{2},1),F{2}-F{5})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_8.fig'])
figure
imagesc(linspace(win(1),win(2),length(F{k})),1:size(F{3},1),F{3}-F{6})
saveas(gcf,[cellbasepath,'\',animalID,'\',animalID,'_',Partition,'_',Trigger,'_',Signal,'_9.fig'])
end
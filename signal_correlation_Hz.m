function drifting = signal_correlation_Hz(i_neur)
drifting = [];
global CELLIDLIST
% signal correlation analysis for one neuron tpye
    %     i_neur = i_da; %
    %     i_neur(i_neur == 3049) = []; % drifting neurons
    %     i_neur(i_neur == 2847) = [];
    %     i_neur(i_neur == 4741) = [];    
    %     i_neur = i_ach; 
    p_window = [-2,0.8];
    respwind = 0.2;
    sigma = 0.02;
    sr = 1000;
    Resps = [];
    Tags = [];
    Soundtype = [];
    Curves = [];
    for j = 1:length(i_neur)
        %try
            [psth_a,spsth_a,~,tags,X,spsth_1st, spsth_2nd, psth_1st] = ultimate_psth_kb2(CELLIDLIST{i_neur(j)},'trial','StimulusOn', p_window,'parts','#Soundfrequency','sigma',sigma);   
        %catch
        %    j
        %    continue
        %end
        %spsth_a(1,:) = spsth_1st(1,:);
        %psth_a(1,:) = psth_1st(1,:);
             
    if iscell(X)
    part_lengths = cell2mat(cellfun(@size,X,'UniformOutput',false));
    part_lengths = part_lengths(1:2:end);
    psth_a = psth_a (part_lengths>29,:);
    spsth_a = spsth_a (part_lengths>29,:);
    tags = tags (part_lengths>29);
    X = X(part_lengths>29);
    else
        continue
    end
    if isempty(tags)
       continue 
    end
    % load session data and find out the sound type
    
    VE = loadcb(CELLIDLIST{i_neur(j)},'TrialEvents');   % load events
    [COMPTRIALS, ~] = partition_trials(VE,'#Soundfrequency');

    soundtype2 = [];
    for k=1:length(COMPTRIALS)
        soundtype2 (k) = VE.TrialType(COMPTRIALS{k}(1));
    end
    soundtype2 = soundtype2(part_lengths>29);
    tags2 = [];
    delinx = [];
    for i= 1: length(tags)
        tags2(i) = str2double(tags{i}(end));
        
        if tags2(i)==1
            delinx = [delinx,i];
        end
    end
    P = [];
    base4 = sum(X{1}(:,1:p_window(1)*-1*sr),2);
    base1 = sum(X{1}(1:round(size(X{1},1)/2),1:p_window(1)*-1*sr),2);
    base2 = sum(X{1}(round(size(X{1},1)/2):end,1:p_window(1)*-1*sr),2);
    for pinx = 2:length(X)-1
        basep = sum(X{pinx}(:,1:p_window(1)*-1*sr),2);
        P(pinx-1) = ranksum(base4,basep);
        if min(sum(base1),sum(base2))==0
        P(pinx-1) = 0; 
        end
    end
    delinx = [delinx,find(P < 0.001)+1];
    if ~isempty(delinx)
        drifting = [drifting, i_neur(j)];
        CELLIDLIST(i_neur(j))
    end
    psth_a(delinx,:) = [];
    spsth_a(delinx,:) = [];
    soundtype2(delinx) = [];
    tags2(delinx) = [];
    

%     [~,I_a] = max(spsth_a(1,respwind));
%     I_a = I_a + respwind(1);
%     if     sum(psth_a(1,I_a-window:I_a+window),2) == 0
%         continue
%     end
%     
    %Resps = [Resps;sum(psth(:,window),2)./sum(psth(:,baseline),2)*factor];
    Resps = [Resps;max(spsth_a(:,(-p_window(1)+0.000)*sr:(-p_window(1))*sr+respwind*sr),[],2)./max(spsth_a(1,-p_window(1)*sr:(-p_window(1))*sr+respwind*sr))];
    Curves = [Curves;spsth_a./max(spsth_a(1,-p_window(1)*sr:(-p_window(1))*sr+respwind*sr),[],2)];
    
    

    Tags = [Tags,tags2];
    Soundtype = [Soundtype,soundtype2];
    if length(Tags)~=length(Soundtype)
        j
    end
    end

    
    
    Soundtype(Tags==1)=[];
    Resps(Tags==1)=[];
    Tags(Tags==1)=[];
    %Curves(Tags==1,:)=[];
    
    Tags(Tags==0) = 10;
    Tens=find(Tags==10)
    Fours=find(Tags==4)
    
    
    Neuronnum = [1];
    for j = 1:length(Fours)
        Neuronnum(end+1:Fours(j)-1) = j; 
    end
    Neuronnum(end+1:length(Tags)) = j;
    
    figure
    hold on
    for i_col = 1:2
    for i = 1:max(Neuronnum)
        line(Tags(Neuronnum == i & Soundtype == i_col),Resps(Neuronnum == i & Soundtype == i_col),'Color',[0.3,0.3,0.3,0.3])
    end
    Colors = [0,1,0;1,0,0];
    

    hold on
    for i = 4:10
        scatter(ones(1,sum(Tags==i & Soundtype == i_col))*i,Resps(Tags==i & Soundtype == i_col),15,Colors(i_col,:))
        if i_col == 1 & i < 10
            sigstar([i,i+1],ranksum(Resps(Tags==i),Resps(Tags==i+1)))
        end
    end
    end
    xlim([4,10])
%     delinx = find(sum(isnan(Curves'))~=0 | sum(isinf(Curves'))~=0);
%     Tags(delinx) = [];
%     Soundtype(delinx) = [];
%     Curves(delinx,:) = [];
%     Neuronnum(delinx) = [];
%     Resps(delinx) = [];
    
    figure
    hold on
    for i = 4:10
    %for i_col = 1:2
    plot(p_window(1):1/sr:p_window(2),mean(Curves(Tags==i & Soundtype == 1,:)),'Color',-1*((Colors(1,:)*(i-3)/7)-[0,1,0]),'LineWidth',3,'DisplayName',[num2str(i),'kHz go'])
    plot(p_window(1):1/sr:p_window(2),mean(Curves(Tags==i & Soundtype == 2,:)),'Color',Colors(2,:)*(i-3)/7,'LineWidth',3,'DisplayName',[num2str(i),'kHz no-go'])
    %end
    end
    xlim([-0.2,0.8])

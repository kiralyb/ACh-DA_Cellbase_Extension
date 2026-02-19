function [P,R]=noise_correlation_df(i_a,i_b,isdoublepart,isplot)

load('L:\_CellBases\ACh_DA_Cellbase\CCG_new\all\stim_light_filter2\CCG_matrices.mat')
crosspairs=[];

for i = 1: length(PairOfCells)
    a1 = ismember(findcellpos(PairOfCells(i,1)),i_a);
    a2 = ismember(findcellpos(PairOfCells(i,2)),i_a);
    d2 = ismember(findcellpos(PairOfCells(i,2)),i_b);
    d1 = ismember(findcellpos(PairOfCells(i,1)),i_b);
    if a1 & d2
       crosspairs = [crosspairs,i]; 
    end
    
end
partition = 1;
window_a = 10; 
window_d = 40;
sigma_a = 0.01;
sigma_d = 0.02;
ach_wind = [-0.0 0.2];
da_wind = [0.0 0.4];
    
    pairs = crosspairs;
    for i = 1:length(pairs)
        [psth_a,spsth_a,~,tags,raszt_a] =  ultimate_psth(PairOfCells{pairs(i),1},'trial','DeliverFeedback',ach_wind,'parts','#Type','sigma',sigma_a);
        [psth_d,spsth_d,~,~,raszt_d] =  ultimate_psth(PairOfCells{pairs(i),2},'trial','DeliverFeedback',da_wind,'parts','#Type','sigma',sigma_d);
        for tag_inx= 1: length(tags)
            tags2(tag_inx) = str2double(tags{tag_inx}(end));
        end
        partition4Hz = find(tags2==2,1,'first');
        partition = find(tags2==2,1,'first');
        %[~,I_a] = max(spsth_a(partition4Hz,:));
        [peaks,locs,width,prominance] = findpeaks(spsth_a(partition4Hz,:));
        [~,ind] = max(prominance);
        I_a = locs(ind);
        window_a = ceil(width(ind));
        response_a = (sum(raszt_a{partition}(:,max([1,I_a-window_a]):min([length(spsth_a),I_a+window_a])),2));
        [peaks,locs,width,prominance] = findpeaks(spsth_d(partition4Hz,:));
        [~,ind] = max(prominance);
        I_d = locs(ind);
        window_d = ceil(width(ind));
        response_d = (sum(raszt_d{partition}(:,max([1,I_d-window_d]):min([length(spsth_d),I_d+window_d])),2));
        meanr_a = mean(response_a);
        meanr_d =  mean(response_d);
        if isdoublepart
            [psth_a,spsth_a,~,tags,raszt_a] =  ultimate_psth(PairOfCells{pairs(i),1},'trial','StimulusOn',ach_wind,'parts','#Type3','sigma',sigma_a);
            [psth_d,spsth_d,~,~,raszt_d] =  ultimate_psth(PairOfCells{pairs(i),2},'trial','StimulusOn',da_wind,'parts','#Type3','sigma',sigma_d);
            for tag_inx= 1: length(tags)
                tags2(tag_inx) = str2double(tags{tag_inx}(end));
            end
            partition4Hz = find(tags2==1,1,'first');
            partition = find(tags2==2,1,'first');
            VE = loadcb(PairOfCells{pairs(i),1},'TrialEvents');   % load events
            trialnum1 = find(VE.Type3==2);
            trialnum2 = find(VE.CorrectRejection==1);
            trialinx = ismember(trialnum1,trialnum2);
            %response_a2 = (sum(raszt_a{partition}(trialinx,I_a-window_a:I_a+window_a),2));
            response_a2 = (sum(raszt_a{partition}(trialinx,max([1,I_a-window_a]):min([length(spsth_a),I_a+window_a])),2));        
            %response_d2 = (sum(raszt_d{partition}(trialinx,I_d-window_d:I_d+window_d),2));
            response_d2 = (sum(raszt_d{partition}(trialinx,max([1,I_d-window_d]):min([length(spsth_d),I_d+window_d])),2));
        
            meanr_a2 = mean(response_a2);
            meanr_d2 =  mean(response_d2);
            
            x = [response_a-round(meanr_a);response_a2-round(meanr_a2)];
            y = [response_d-round(meanr_d);response_d2-round(meanr_d2)];
            
        else
            x = response_a-meanr_a;
            y = response_d-meanr_d;
        end
        [R(i),P(i)] = corr(x,y,'Type','Spearman');
        if isplot
        xx = x + (rand(length(x),1)-0.5)*.5;
        yy = y + (rand(length(y),1)-0.5)*.5;
        figure
        scatter(xx(1:length(response_a)),yy(1:length(response_a)),'x','k','MarkerEdgeAlpha',0.5)
        hold on
        if isdoublepart
        scatter(xx(length(response_a):end),yy(length(response_a):end),'x','r','MarkerEdgeAlpha',0.5)
        end
        coef = polyfit(x, y, 1);
        refline(coef(1), coef(2));
        setmyplot_balazs;
        ylim = get(gca,'ylim');
        xlim = get(gca,'xlim');
        text(xlim(1),ylim(2) * 0.95,sprintf('r=%f',R(i)))
        text(xlim(1),ylim(2) * 0.85,sprintf('p=%f',P(i)))
        end
    end
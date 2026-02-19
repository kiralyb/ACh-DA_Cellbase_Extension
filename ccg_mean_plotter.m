figure
for col = 1:3
    switch col
        case 1
            load('L:\_CellBases\ACh_DA_Cellbase\CCG_final\Cue\CCG_matrices.mat')
            %L = length(CCR);
            lags = -250:250;
            CCR = CCR(:,lags+1+(size(CCR,2)-1)/2);
            CCR_Poi = CCR_Poi(:,lags+1+(size(CCR_Poi,2)-1)/2);
            %mlag{1} = 300:350;
            %mlag{2} = 225:275;
        case 2
            load('L:\_CellBases\ACh_DA_Cellbase\CCG_final\Hit\CCG_matrices.mat')
            %L = length(CCR);
            lags = -250:250;
            CCR = CCR(:,lags+1+(size(CCR,2)-1)/2);
            CCR_Poi = CCR_Poi(:,lags+1+(size(CCR_Poi,2)-1)/2);
            %mlag{1} = 300:350;
            %mlag{2} = 225:275;
        case 3
            load('L:\_CellBases\ACh_DA_Cellbase\CCG_final\FA\CCG_matrices.mat')
            %L = length(CCR);
            lags = -250:250;
            CCR = CCR(:,lags+1+(size(CCR,2)-1)/2);
            CCR_Poi = CCR_Poi(:,lags+1+(size(CCR_Poi,2)-1)/2);
            %mlag{1} = 350:400;
            %mlag{2} = 225:275;
        case 4
            load('L:\_CellBases\ACh_DA_Cellbase\CCG_new\all\behav6\CCG_matrices.mat')
            %lags = -1000:1000;
     end
i_pairs=[];
e_pairs = [];

flipi = [];
swin = 20;
for i = 1: length(PairOfCells)
    a1 = ismember(findcellpos(PairOfCells(i,1)),i_ach);
    de = ismember(findcellpos(PairOfCells(i,2)),i_da_e);
    di = ismember(findcellpos(PairOfCells(i,2)),i_da_i);
    
    if (a1 & di) 
       i_pairs = [i_pairs,i]; 
    end
    if (a1 & de)
        e_pairs = [e_pairs, i];
    end
end
if col == 3
    %i_pairs([13,14,15]) = [];
    %i_pairs([10,13,14,15]) = [];
end
%i_pairs = i_pairs([10,11,13,14,15,23]);
pairstypes = {i_pairs,e_pairs};
Colors = {'b','g','r','k'};
for i = 1:2
CCR_norm2 = bsxfun(@rdivide,CCR,sum(CCR,2));
CCR_Poi_norm = bsxfun(@rdivide,CCR_Poi,sum(CCR_Poi,2));
CCR_norm = (CCR_norm2 - CCR_Poi_norm);
%CCR_norm = (CCR - CCR_Poi)./sqrt(CCR_Poi + 1) ./ SegmentLength;
%CCR_norm = (CCR_norm2 - CCR_Poi_norm)./sqrt(CCR_Poi+1);
%CCR_norm = bsxfun(@rdivide,CCR_norm,)
%CCR_norm = bsxfun(@rdivide,CCR-CCR_Poi,1)%(CCR-CCR_Poi);%bsxfun(@rdivide,CCR-CCR_Poi,sum(CCR,2));%CCR-CCR_Poi;%1%sqrt(CCR_Poi)%sum(CCR_Poi,2))
%CCR_norm(flipi,:)=flip(CCR_norm(flipi,:),2);
%CCR_smooth = movmean(CCR_norm,50,2)
%CCR_smooth = smoothdata(CCR_norm(:,1:1:end),2,'gaussian',70);

 CCR_smooth = nan(size(CCR_norm));
  for ind = 1:size(CCR_norm,1)
      if sum(CCR(ind,:),2) > 0
          CCR_smooth(ind,:) = smoothdata(CCR_norm(ind,:),2,'gaussian',length(lags)/sqrt(SegmentLength(ind)));
      end
  end
subplot(2,2,i)
hold on
%plot(lags(1:1:end),movmean(nanmedian(CCR_smooth(pairstypes{i},:)),20),Colors{col},'LineWidth',2)
%errorshade(lags(1:1:end),movmean(nanmedian(CCR_smooth(pairstypes{i},:)),25),iqr(CCR_smooth(pairstypes{i},:))/sqrt(length(pairstypes{i}))/1.35,'LineColor',Colors{col},'ShadeColor',Colors{col})
%errorshade(lags(1:1:end),movmean(nanmedian(CCR_smooth(pairstypes{i},:)),30),std(CCR_smooth(pairstypes{i},:))/sqrt(length(pairstypes{i})),'LineColor',Colors{col},'ShadeColor',Colors{col})
errorshade(lags(1:1:end),nanmean(CCR_smooth(pairstypes{i},:)),nanstd(CCR_smooth(pairstypes{i},:))/sqrt(length(pairstypes{i})),'LineColor',Colors{col},'ShadeColor',Colors{col})


% %xlim([-300,300])
line([0,0],[-0.001,0.002])
line([lags(1),lags(end)],[0,0])
%ylim([-0.0015,0.0015])
%xlim([-225,225])
setmyplot_balazs

subplot(2,2,i+2)
hold on
% [~,mind] = max(abs(CCR_smooth(pairstypes{i},250:450)),[],2);
% mind = mind + 249;
%B = [];
%for b = 1:length(mind)
    %B = nanmean(CCR_norm(pairstypes{i},mind-50:mind),2)
%end
[~,mind] = max(abs(nanmean(CCR_smooth(pairstypes{i},:))),[],2);
mind = mind;
boxplot(mean(CCR_norm(pairstypes{i},mind-25:mind+25),2),'position',col)
ylim([-.002,0.002])
line([0,3],[0,0])

signrank(mean(CCR_norm(pairstypes{i},mind-25:mind+25),2))
% %[~,sorter] = sort(-sum(CCR_smooth(pairstypes{i},800:900)') + sum(CCR_smooth(pairstypes{i},1001:1100)')-sum(CCR_smooth(pairstypes{i},1150:1300)'));
% [~,sorter] = sort(sum(CCR_smooth(pairstypes{i},300:400)'));
% subplot(1,3,i)
% imagesc(lags,1:length(pairstypes{i}),CCR_smooth(pairstypes{i}(sorter),:))
% %caxis([0,0.0010]);
% xlabel('lag (ms)')
% yyaxis right
% plot(lags,nanmean(CCR_smooth(pairstypes{i},1:end)),'r','LineWidth',2)
% hold on
% %xlim([-300,300])
% %line([0,0],[0.0000,0.001],'Color','w','LineWidth',1);
% %ylim([0,0.001])
% %line([0,0],[1,length(pairstypes{i})],'Color','w','LineWidth',1)
end

end
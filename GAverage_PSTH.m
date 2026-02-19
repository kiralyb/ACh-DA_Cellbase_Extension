function [p,RESPEAKL] = GAverage_PSTH(cellids,trigger,partitions,sigma,dt,wind,baselinewind,respwind,invfact,isplot)
dbstop if error  
    if nargin < 10
      isplot = 1;
  end
  spsth_T2 = nan(length(cellids),4,((wind(2)-wind(1))/dt)+1);
  time = wind(1):dt:wind(2);
  baselinet = ((baselinewind-wind(1))/dt)+1;
  baselinet = baselinet(1):baselinet(2);
  respt = ((respwind-wind(1))/dt)+1;
  respt = respt(1):respt(2);
  sizes = nan(length(cellids),4);
  for  j = 1:length(cellids)

      try % CorrectReject miatt
      [~,spsth_T2_temp,~,tags,spt] =  ultimate_psth(cellids{j},'trial',trigger,wind,'parts',partitions,'sigma',sigma,'dt', dt);
      if strcmp(partitions,'#FalseAlarm') | strcmp(partitions,'#Hit') | strcmp(partitions,'#CorrectRejection') | strcmp(partitions,'Miss')  
      if min(size(spt))>1
      spsth_T2(j,1,:) = (spsth_T2_temp(1,:) - mean(spsth_T2_temp(1,baselinet),2)) ./ std(spsth_T2_temp(1,baselinet),[],2);    
      end
      elseif length(tags) == 1
          continue
      else
      index4Hz = find(strcmp(tags, [partitions(2:end), '=', num2str(1)]));
      for part = 1:4
        index = find(strcmp(tags, [partitions(2:end), '=', num2str(part)]));
        if ~isempty(index) 
            sizes(j,part)=size(spt{index},1);
            if size(spt{index},1)>1
            spsth_T2(j,part,:) = (spsth_T2_temp(index,:) - mean(spsth_T2_temp(index4Hz,baselinet),2)) ./ std(spsth_T2_temp(index4Hz,baselinet),[],2);   
            end
        end
      end
      end
      %k = k + 1
      catch
      %j
      end
  end
spsth_T2(spsth_T2==inf) = nan; 
if isplot
figure
subplot(1,2,1)
hold on
colors = [0,1,0;1,0,0;0,0.5,0;0.5,0,0];
for plo=1:4
  errorshade(time,squeeze(nanmean(spsth_T2(:,plo,:))),squeeze(nanstd(spsth_T2(:,plo,:)))/sqrt(size(spsth_T2,1)),'LineColor',colors(plo,:),'ShadeColor',colors(plo,:)) 
end
setmyplot_balazs 
subplot(1,2,2)
%plot(max(spsth_T2(:,:,1000:1200),[],3)','o')
hold on
%if invfact == -1
    spsth_T2(:,[2,4],:) = spsth_T2(:,[2,4],:)*invfact;
%    [RESPAMP,RESPEAKL]=max(spsth_T2(:,:,respt)-spsth_T2(:,:,round((respt(1)-sigma/dt))),[],3,'omitnan');
%else
     [RESPAMP,RESPEAKL]=max(spsth_T2(:,:,respt),[],3);
%end
RESPEAKL(isnan(RESPAMP))=nan;
RESPAMP(:,[2,3,4])=RESPAMP(:,[3,4,2]);
RESPAMP(:,[3,4])=RESPAMP(:,[3,4]) * invfact;
sizes(:,[2,3,4])=sizes(:,[3,4,2]);


boxplot(RESPAMP,'symbol','')
h=findobj('LineStyle','--'); set(h, 'LineStyle','-');


%plot(RESPAMP','g')
p = nan;
if ~(strcmp(partitions,'FalseAlarm') | strcmp(partitions,'Hit'))
filt = sizes(:,1)>7 & sizes(:,2)>7;
p(1) = signrank(RESPAMP(filt,1),RESPAMP(filt,2));
try
filt = sizes(:,3)>7 & sizes(:,4)>7;
p(2) = signrank(RESPAMP(filt,3),RESPAMP(filt,4));
catch
filt = sizes(:,2)>7 & sizes(:,4)>7;
p(2) = signrank(RESPAMP(filt,2),RESPAMP(filt,4));
end
sigstar({[1,2],[3,4]},p)
setmyplot_balazs


%figure
%boxplot(RESPAMP(:,4)-RESPAMP(:,3))
end
end
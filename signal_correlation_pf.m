function signal_correlation_pf(i_neur,type,plotinx,inverter)
if nargin < 3
   plotinx = 1 
end
if nargin < 4
    inverter = 1;
end

  global CELLIDLIST
  cue_wind = [-0.3,0.8];
  respwind = 0.20;
  sr = 1000;
  for j = 1:length(i_neur)

  VE = loadcb(CELLIDLIST{i_neur(j)},'TrialEvents');   % load events
  newcues = VE.Soundtype==3;
  oldcues = VE.Soundtype == 1;
    resps = ~isnan(VE.Type(newcues));
  performance = smoothdata(resps,"gaussian",7);
  
  %[COMPTRIALS, ~] = partition_trials(VE,'#Type3');

  %reactiontime = VE.DeliverFeedback(newcues) - VE.StimulusOff(newcues) ;%VE.Delay(newcues);
  switch type
      case ('cue')
          NC = find(newcues);
          rop = ones(1,length(NC));
          trigger = 'StimulusOn';
          sigma_cue = 0.02;
          control = VE.Soundtype==1;
          Comptrials{1} = NC(find(performance < 0.4 & rop));
          Comptrials{2} = NC(find(performance >= 0.4 & performance <= 0.6 & rop));
          Comptrials{3} = NC(find(performance > 0.6 & rop));
          Comptrials{4} = find(control);
          LEGENDS = {'below 40%','between 40% and 60%','above 60%'};
          ycolors = [0.5,0.5,0;0.7,0.7,0;0.9,0.9,0];
          norm = [];
          Xlab = {'Time from cue (s)'};
      case ('cueo')
          NC = find(VE.Soundtype==2);
          rop = ones(1,length(NC));
          trigger = 'StimulusOn';
          sigma_cue = 0.02;
          control = VE.Soundtype==1;
          if ~isfield(VE, 'thrdSound')
              VE.thrdSound = nan(1,VE.NTrials(1));
          end
          thrdsoundorder = cumsum([true, diff(VE.thrdSound) ~= 0]);
          for q = 1:max(thrdsoundorder)
              try
             thrdsoundtype(q)=mean(VE.TrialType(thrdsoundorder==q & VE.Soundtype == 3)) ;
              catch
              thrdsoundtype(q)= 0;
              end
          end
          thrdsoundtype2 = arrayfun(@(x) thrdsoundtype(x), thrdsoundorder);
          
          Comptrials{1} = NC(find( thrdsoundtype2(NC) == 1));
          Comptrials{2} = NC(find( thrdsoundtype2(NC) == 2));
          Comptrials{3} = find(control);
          LEGENDS = {'reward','punihsment'};
          ycolors = [0.5,0.5,0;0.9,0.9,0];
          norm = [];
          Xlab = {'Time from cue (s)'};
      case ('rewp')
          NC = find(newcues);
          rop = VE.Hit(newcues)==1;
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          control = rop;
          Comptrials{1} = NC(find(performance < 0.5 & rop));
          %Comptrials{2} = NC(find(performance >= 0.4 & performance <= 0.6 & rop));
          Comptrials{2} = NC(find(performance > 0.5 & rop));

          Comptrials{3} = NC(find(control));
          %LEGENDS = {'below 40%','between 40% and 60%','above 60%'};
          %ycolors = [0,0.5,0.1;0,0.7,0.2;0,1,0.3];
          norm = [];
          Xlab = {'Time from reward (s)'};
          LEGENDS = {'below 50%','above 50%'};
          ycolors = [0.0,0.6,0;0.0,0.9,0.5];
      case('rewr')
          NC = find(oldcues);
          control = find(VE.Hit(oldcues)==1);
          rop = VE.Hit(oldcues)==1;
          reactiontime = VE.StimulusDuration(oldcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          Comptrials{1} = NC(find(reactiontime < quantile(reactiontime(rop),1/3) & rop));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),2/3) & reactiontime >= quantile(reactiontime(rop),1/3) & rop));
          Comptrials{2} = NC(find(reactiontime > quantile(reactiontime(rop),2/3) & rop));
          Comptrials{3} = NC(find(control));
          %LEGENDS = {'fast','medium','slow'};
          %ycolors = [0,1,0.9;0,0.9,0.7;0,0.8,0.5];
          norm = [];
          Xlab = {'Time from reward (s)'};
          LEGENDS = {'fast','slow'};
          ycolors = [0,1,0.9;0,0.8,0.5];
      case('rewr2')
          NC = find(newcues);
          control = find(VE.Hit(newcues)==1);
          rop = VE.Hit(newcues)==1;
          reactiontime = VE.StimulusDuration(newcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          Comptrials{1} = NC(find(reactiontime < quantile(reactiontime(rop),1/2) & rop));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),2/3) & reactiontime >= quantile(reactiontime(rop),1/3) & rop));
          Comptrials{2} = NC(find(reactiontime > quantile(reactiontime(rop),1/2) & rop));
          Comptrials{3} = NC(find(control));
          %LEGENDS = {'fast','medium','slow'};
          %ycolors = [0,1,0.9;0,0.9,0.7;0,0.8,0.5];
          norm = [];
          Xlab = {'Time from reward (s)'};
          LEGENDS = {'fast','slow'};
          ycolors = [0,1,0.9;0,0.8,0.5];
      case('rewrn')
          NC = find(oldcues);
          reactiontime = VE.StimulusDuration(oldcues);
          rop = (VE.Hit(oldcues)==1);
          
          delaytime = VE.DeliverFeedback(oldcues) - VE.StimulusOff(oldcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          control = rop;
          Comptrials{1} = NC(find(delaytime < quantile(delaytime(rop),1/6) & rop));
          %Comptrials{1} = NC(find(reactiontime <= quantile(reactiontime(rop),3/4) & reactiontime >= quantile(reactiontime(rop),1/4) & rop & delaytime < quantile(delaytime(rop),1/3)));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),3/4) & reactiontime >= quantile(reactiontime(rop),1/4) & rop & delaytime > quantile(delaytime(rop),2/3)));
          Comptrials{2} = NC(find(delaytime <= quantile(delaytime(rop),5/6) & delaytime >= quantile(delaytime(rop),1/6) & rop));
          Comptrials{3} = NC(find(delaytime > quantile(delaytime(rop),5/6) & rop));
          Comptrials{4} = NC(find(control));
          %LEGENDS = {'fast','medium','slow'};
          %ycolors = [0,1,0.9;0,0.9,0.7;0,0.8,0.5];
          norm = [];
          Xlab = {'Time from reward (s)'};
          LEGENDS = {'short','average','long'};
          ycolors = [0,1,0.9;0,0.9,0.7;0,0.8,0.5];
      case ('pun')
          NC = find(newcues);
          rop = VE.FalseAlarm(newcues)==1;
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;  
          control = rop;
          Comptrials{1} = NC(find(performance < 0.5 & rop));
          %Comptrials{2} = NC(find((performance >= 0.4 & performance <= 0.6) & rop));
          Comptrials{2} = NC(find(performance > 0.5 & rop));
          Comptrials{3} = NC(find(control));
          %LEGENDS = {'below 40%','between 40% and 60%','above 60%'}; ycolors = [1,0.5,0;0.75,0.4,0;0.5,0.3,0];            
          LEGENDS = {'below 50%','above 50%'}; ycolors = [1,0.5,0;0.5,0.3,0];
          norm = [];%getvalue('firing_rate',CELLIDLIST{i_neur(j)});
          Xlab = {'Time from punishment (s)'};
      case ('punr')
          NC = find(oldcues);
          rop = VE.FalseAlarm(oldcues)==1;

          reactiontime = VE.StimulusDuration(oldcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          control = rop;

          Comptrials{1} = NC(find(reactiontime < quantile(reactiontime(rop),1/2) & rop));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),2/3) & reactiontime >= quantile(reactiontime(rop),1/3) & rop));
          Comptrials{2} = NC(find(reactiontime > quantile(reactiontime(rop),1/2) & rop));
          Comptrials{3} = NC(find(control));

          LEGENDS = {'fast','slow'};
          ycolors = [1,0.5,0;0.5,0.3,0];
          norm = [];%getvalue('firing_rate',CELLIDLIST{i_neur(j)});
          Xlab = {'Time from punishment (s)'};
      case ('punr2')
          NC = find(newcues);
          rop = VE.FalseAlarm(newcues)==1;

          reactiontime = VE.StimulusDuration(newcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          control = rop;

          Comptrials{1} = NC(find(reactiontime < quantile(reactiontime(rop),1/2) & rop));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),2/3) & reactiontime >= quantile(reactiontime(rop),1/3) & rop));
          Comptrials{2} = NC(find(reactiontime > quantile(reactiontime(rop),1/2) & rop));
          Comptrials{3} = NC(find(control));

          LEGENDS = {'fast','slow'};
          ycolors = [1,0.5,0;0.5,0.3,0];
          norm = [];%getvalue('firing_rate',CELLIDLIST{i_neur(j)});
          Xlab = {'Time from punishment (s)'};
      case('punrn')
          %newcues = VE.FalseAlarm==1;  
          NC = find(newcues);
          reactiontime = VE.StimulusDuration(newcues);
          rop = (VE.FalseAlarm(newcues)==1);
          delaytime = VE.DeliverFeedback(newcues) - VE.StimulusOff(newcues);
          trigger = 'DeliverFeedback';
          sigma_cue = 0.01;
          control = rop;
          Comptrials{1} = NC(find(delaytime < quantile(delaytime(rop),1/5) & rop));
          %Comptrials{1} = NC(find(reactiontime <= quantile(reactiontime(rop),3/4) & reactiontime >= quantile(reactiontime(rop),1/4) & rop & delaytime < quantile(delaytime(rop),1/3)));
          %Comptrials{2} = NC(find(reactiontime <= quantile(reactiontime(rop),3/4) & reactiontime >= quantile(reactiontime(rop),1/4) & rop & delaytime > quantile(delaytime(rop),2/3)));
          Comptrials{2} = NC(find(delaytime <= quantile(delaytime(rop),4/5) & delaytime >= quantile(delaytime(rop),1/5) & rop));
          Comptrials{3} = NC(find(delaytime > quantile(delaytime(rop),4/5) & rop));
          Comptrials{4} = NC(find(control));
          norm = [];
          Xlab = {'Time from reward (s)'};
          LEGENDS = {'short','average','long'};
          ycolors = [1,0.5,0;0.75,0.4,0;0.5,0.3,0];
  end
   if min(cellfun(@length,Comptrials))<5
      spsth_cue(j,:,:) = nan(length(Comptrials),diff(cue_wind)*sr+1);
      continue
   end

  [~,spsth_cue(j,:,:)] =  ultimate_psth_spectrials(CELLIDLIST{i_neur(j)},'trial',trigger,cue_wind,Comptrials,'sigma',sigma_cue);

  %if inverter == 1
      norm= max(spsth_cue(j,end,round((-cue_wind(1)*sr):round((-cue_wind(1)+respwind)*sr))));
      spsth_cue(j,:,:) = spsth_cue(j,:,:) ./ norm; %./ max(spsth_cue(j,end,:));
  %else
  %    norm = getvalue('firing_rate',CELLIDLIST{i_neur(j)});
  %    spsth_cue(j,:,:) = spsth_cue(j,:,:) ./ norm; 
  %end
end

%figure
%cue_wind(1) = -0.10;
RESP = max((spsth_cue(:,:,round(-cue_wind(1)/0.001):round((respwind-cue_wind(1))/0.001))-spsth_cue(:,:,round((-cue_wind(1)-sigma_cue)/0.001)))*inverter,[],3)*inverter;
keepinx = sum(isnan(RESP(:,[1:end-1]))')==0;
RESP = RESP(keepinx,1:length(Comptrials)-1);
spsth_cue=spsth_cue(keepinx,:,:);
cue_wind(1) = -0.2;
subplot(2,4,plotinx)
hold on
for p =  1:length(Comptrials)-1
errorshade(linspace(cue_wind(1),cue_wind(2),length(spsth_cue)),squeeze(mean(spsth_cue(:,p,:)))',squeeze(std(spsth_cue(:,p,:)))'/sqrt(sum(keepinx)),'LineColor',ycolors(p,:),'ShadeColor',ycolors(p,:))
end
xlabel(Xlab)
ylabel('Average normalized firing rate')
legend(LEGENDS);
setmyplot_balazs
subplot(2,4,plotinx+4)
hold on
plot(RESP','Color',[1,1,1,1]/5)
ylabel('Normalized response')
xlim([1,length(Comptrials)-1])
plot(nanmean(RESP),'r')
xticks(1:length(Comptrials)-1)
xticklabels(LEGENDS)

if length(Comptrials)-1 == 3
    sigstar([1,2],signrank(RESP(:,1),RESP(:,2)))
    sigstar([2,3],signrank(RESP(:,2),RESP(:,3)))
    sigstar([1,3],signrank(RESP(:,1),RESP(:,3)))
else
    sigstar([1,length(Comptrials)-1],signrank(RESP(:,1),RESP(:,end)))
end
setmyplot_balazs

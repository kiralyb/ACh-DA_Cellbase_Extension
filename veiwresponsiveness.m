function veiwresponsiveness(ID)
TE = loadcb(ID,'TrialEvents')
smoother = 7;
meanwindow = [smoother,smoother];
smoothmet = "gaussian";
newcues = TE.TrialCode>2;

resps = ~isnan(TE.Type);
performance_bin = resps(newcues);
performance = smoothdata(performance_bin,smoothmet,meanwindow);
figure; plot(performance)
ylim([0,1])
borders = find(diff(TE.thrdSound(newcues))~=0);
hold on
line([borders,borders],[0,1])
scatter(1:length(performance_bin),performance_bin,'|')
setmyplot_balazs


function [Hit, FA, RT, Mstate,trialnum,newcuenum,newtrialnum] = psy_behav_runner(MICE,filter)
dbstop if error
LIST = listtag('session');

for mouseindex = 1:length(MICE)
    MICESES = LIST(ismember(LIST(:,1), MICE(mouseindex)), :);
    for sesindex = 1:length(MICESES)
        [Hit{mouseindex}(sesindex),FA{mouseindex}(sesindex),RT{mouseindex}(sesindex),Mstate{mouseindex}(sesindex),trialnum{mouseindex}(sesindex),newcuenum{mouseindex}(sesindex),newtrialnum{mouseindex}(sesindex)]=FixPsyBehav(MICESES{sesindex,1},MICESES{sesindex,2},filter); %0, 1 or 2        
    end
end

if filter == 0
    
    % S1a
    figure
    minlength = 75;
    nMICE = numel(MICE);
    maxL = 40;

for state = 1:3

    % Preallocate
    MR.Hit = nan(nMICE, maxL);
    MR.FA  = nan(nMICE, maxL);
    MR.RT  = nan(nMICE, maxL);

    for i = 1:nMICE
        mask = Mstate{i} == state & trialnum{i} > minlength;
        L = 1:sum(mask);

        MR.Hit(i,L) = Hit{i}(mask);
        MR.FA(i,L)  = FA{i}(mask);
        MR.RT(i,L)  = RT{i}(mask);
    end

    % Last index with at least 5 data points
    ML = find(sum(~isnan(MR.RT)) <= 4, 1, 'first') - 1;

    % ----- plot boxplot for state 3 -----
    if state == 3
        subplot(2,6,6)
        boxplot([MR.Hit(:), MR.FA(:)])
        ylim([0 1])
        setmyplot_balazs

        subplot(2,6,12)
        boxplot(MR.RT(:))
        ylim([0 0.7])
        setmyplot_balazs

        plotinx = state:5;
    else
        plotinx = state;
    end

    % ----- Hit / FA -----
    subplot(2,6,plotinx)
    hold on
    plot(MR.Hit','.-','Color',[0 1 0])
    plot(MR.FA','.-','Color',[1 0 0])
    plot(nanmean(MR.Hit(:,1:ML)),'g','LineWidth',3)
    plot(nanmean(MR.FA(:,1:ML)),'r','LineWidth',3)
    setmyplot_balazs

    % ----- RT -----
    subplot(2,6,plotinx+6)
    hold on
    plot(MR.RT','.-','Color',[0 0 1])
    plot(nanmean(MR.RT(:,1:ML)),'b','LineWidth',3)
    ylim([0 0.7])
    setmyplot_balazs
end

end
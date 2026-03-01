function populationclust_elbow(area)
load(['L:\_CellBases\ACh_DA_Cellbase\Population_maps\',area,'_Matrix.mat'])
g.baselinelength=g.baselinelength/g.dt;
margin=g.sigma*3;
time=g.window(1):g.dt:g.window(2);
time_plot=g.pwindow(1):g.dt:g.pwindow(2); % plotted time window
for partnum=1:length(g.ClusterPartitions)
    ctime=g.cwindow(partnum,1):g.dt:g.cwindow(partnum,2); % time window considered for clustering
    [~,c1]=(min(abs(time-ctime(1))));
    [~,c2]=(min(abs(time-ctime(end))));
    cinx{partnum,:}=c1:c2;
end
PCA2=[];
VAR=[];
figure; hold on;
for f=1:length(g.ClusterPartitions)

[~,PCA1,~,~,VAR1]=pca(((squeeze(Zpsth(:,f,cinx{f}))))); % calculate principal copmponents of the Z-scored psth in the clustering time window
plot(0:length(VAR1),[0;cumsum(VAR1)])
xlim([0,10])
ylim([0,100])
ylabel('Explained Variance')
xlabel('# of PCs')
setmyplot_balazs
line([0,length(VAR1)],[70,70],'Color','r')
PCA2=[PCA2,PCA1(:,1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components
VAR=[VAR;VAR1(1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components


% Total variance (around global mean)
grandMean = mean(PCA2);
TotalSS = sum(sum((PCA2 - grandMean).^2));

Kmax = 10;
explainedVar = zeros(Kmax,1);

for k = 1:Kmax
    [Clusters, ~, sumd] = kmeans(PCA2, k, 'Replicates', 100);

    WCSS = sum(sumd);   % within-cluster SS

    explainedVar(k) = (1 - WCSS/TotalSS) * 100;
end

figure;
plot(1:Kmax, explainedVar, '-o');
hold on
yline(70, '--r', '70% threshold', 'LineWidth', 1.5);
xlabel('Number of clusters (k)');
ylabel('Explained variance (%)');
title('Elbow Method (% Variance Explained)');
grid on;
ylim([0,100])
xlim([0,8])
setmyplot_balazs

k = find(explainedVar > 70,1,'first')

figure
Clusters = kmeans(PCA2, k, 'Replicates', 100);
SORT_=zeros(1,g.ClusterNum);
for j=1:g.ClusterNum
     %SORT_(j)=mean(integralrespons(Clusters==j,g.ClusterPartitions(1)));
     SORT_(j)=mean(PCA2(Clusters==j,1))   
end
[~,rule]=sort(SORT_);
TEMP=zeros(length(Clusters),1);
for j=1:g.ClusterNum
     TEMP(Clusters==rule(j))=j;
end
scatter3(PCA2(:,1),PCA2(:,2),PCA2(:,3),50,TEMP*100,'.'); % 3D PCA cluster plot

[B,inx]=sortrows([TEMP,PCA2(:,3)]); %inx: order of the sorted cells (original index)
%clusternum=histc(Clusters,1:g.ClusterNum); %number of cells in each cluster
clusterborders=find(diff(Clusters(inx))~=0);

sortedZpsth=Zpsth(inx,:,:);
sortedCellIDS=IDS(inx);


cellnum = size(sortedZpsth,1);
partnum = 2;
figure
for k = 1:2 %rows -- corresponding to different partitions
    subplot(partnum,8,[(k-1)*8+1:(k-1)*8+6])
    N=imagesc(time,1:cellnum,squeeze(sortedZpsth(:,k,:)));
    c_limit=quantile(N.CData(:),[0.9995,0.0005]);
    colormapdefiner(c_limit(1),c_limit(2),0,100,partnum,8,[(k-1)*8+1:(k-1)*8+6],g.Trigger1Name,g.Trigger2Name,'Z-score',-1)
    hold on
    additionalplots(time,time_plot,clusterborders,cellnum(end),g,inx,[],IDS);
    hold off
    subplot(partnum,8,(k-1)*8+7)
    imagesc(PCA2);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
end

for j=1:g.ClusterNum
    P = figure;
    %g.avgplotcolors=['g','r','y','b'];
    %P=subplot(1,1,1)
    for k=1:length(g.PlotPartitions)
        %subplot(g.ClusterNum,length(g.PlotPartitions),(j-1)*length(g.PlotPartitions)+k)
        %P=subplot(length(g.PlotPartitions),1,k)
        borders=[0;clusterborders;cellnum(end)];
        avgZscore(j,k,:)=nanmean(sortedZpsth(borders(j)+1:borders(j+1),k,:),1);
        stdZscore(j,k,:)=std(sortedZpsth(borders(j)+1:borders(j+1),k,:),1);
        set(gcf, 'Renderer', 'painters');
        %errorshade(time, squeeze(avgZscore(j,k,:)),squeeze(stdZscore(j,k,:)), 'LineColor',g.avgplotcolors(k), 'ShadeColor',g. avgplotcolors(k), 'LineWidth', 3);
        plot(time,squeeze(avgZscore(j,k,:)),'Color',g.avgplotcolors(k),'LineWidth', 3);
        hold on
        %min(P.YTick)
        line([0,0],[0,8],'Color','magenta','LineWidth',1.5);
        line([g.avgtimediff,g.avgtimediff],[0,8],'Color','cyan','LineWidth',1.5);
        xlim([time_plot(1),time_plot(end)]);
    end  
end


function colormapdefiner(largest,smallest,indexValue,binnum,s1,s2,f,trig1,trig2,ctype,scalcefactor) 
% logarithmic scale for the auROC map (red - green)
% Calculate where proportionally indexValue lies between minimum and maximum values
index = abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [logspace(0,scalcefactor,round(binnum*index))',...
            zeros(round(binnum*index),1),...
            zeros(round(binnum*index),1)];

customCMap2 = [zeros(round(binnum*(1-index)),1),...
            logspace(scalcefactor,0,round(binnum*(1-index)))',...
            zeros(round(binnum*(1-index)),1)];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
colormap(subplot(s1,s2,f),customCMap);
caxis ([smallest,largest]);
%c2=colorbar;
xlabel(sprintf('time(s) - %s - %s',trig1,trig2));
ylabel('cell#');
%c2.Label.String=ctype;

function additionalplots(time,time_plot,clusternum,cellnum,g,inx,taggednum,IDS)

line([time(1),time(end)],[clusternum+0.5,clusternum+0.5],'Color','yellow');
line([0,0],[0,cellnum+0.5]','Color','magenta','LineWidth',1.5);
line([g.avgtimediff,g.avgtimediff],[0,cellnum]','Color','cyan','LineWidth',1.5);
[t1,t2]=ismember(squeeze(inx),taggednum);
t1=find(t1==1);
t3=IDS(inx(t1));
%t2 = t2(t2~=0);
%[t2,t4]=sort(t2);
if isequal(g.showtagged,'*')
scatter(ones(1,length(taggednum))*g.pwindow(1)+g.dt,t1,[],[0.5843 0.8157 0.9882],'*');
elseif isequal(g.showtagged,'IDS')
ylabel('tagged cell IDs');    
set(gca,'ytick',[t1],'yticklabel',t3);
elseif isequal(g.showtagged,'showall')
ylabel('cell IDs');    
set(gca,'ytick',[1:cellnum],'yticklabel',IDS(inx));
end
xlim([time_plot(1),time_plot(end)])
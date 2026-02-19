load('L:\_CellBases\ACh_DA_Cellbase\Population_maps\ALL_20220816\HDB\DATA.mat')
%%
Zpsth = sortedZpsth;
g.baselinelength=g.baselinelength/g.dt;
margin=g.sigma*3;
time=g.window(1):g.dt:g.window(2);
time_m=g.window(1)-margin:g.dt:g.window(2)+margin; %plotted time window + margin for smoothing
time_plot=g.pwindow(1):g.dt:g.pwindow(2); % plotted time window
time_calc=g.window(1)-g.maxtimediff-margin:g.dt:g.window(2)+g.maxtimediff+margin; % extended time window for calculation
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
if g.norm_method=='Z-score'
[~,PCA1,~,~,VAR1]=pca(((squeeze(Zpsth(:,f,cinx{f}))))); % calculate principal copmponents of the Z-scored psth in the clustering time window
plot(0:length(VAR1),[0;cumsum(VAR1)])
xlim([0,10])
ylim([0,100])
ylabel('Explained Variance')
xlabel('# of PCs')
setmyplot_balazs
line([0,length(VAR1)],[70,70],'Color','r')
%PCA1=pca(((squeeze(Zpsth(:,f,cinx)))')); % calculate principal copmponents of the Z-scored psth in the clustering time window
elseif g.norm_method=='auROC'
PCA1=pca(((squeeze(auROC(:,f,floor(cinx{f}(1)/g.ROCWindow):floor(cinx{f}(end)/g.ROCWindow))))')); % calculate principal copmponents of the auROC in the clustering time window
end
%PCA2(:,((f-1)*g.PCA_dim+1):((f-1)*g.PCA_dim+1)+(g.PCA_dim-1))=PCA1(:,1:g.PCA_dim); % take the 1st x=PCA_dim principal components
PCA2=[PCA2,PCA1(:,1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components
VAR=[VAR;VAR1(1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components
end
%%

% Total variance (around global mean)
grandMean = mean(PCA2);
TotalSS = sum(sum((PCA2 - grandMean).^2));

Kmax = 10;
explainedVar = zeros(Kmax,1);

for k = 1:Kmax
    [idx, C, sumd] = kmeans(PCA2, k, 'Replicates', 100);

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

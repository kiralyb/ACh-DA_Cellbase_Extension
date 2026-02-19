function jpsth(PairOfCells,ie_pairs,trigger,partition,wind,binSize)
isplotex = 0;

%binSize = 25; % number of original bins to combine
%wind = [0,0.76];
W = round(diff(wind)*1000/binSize);
jointMatrixT = nan(73,W,W);
InteractionT = nan(73,W,W);
ccgT  =  nan(73,W*2+1);
InteractionCCGT  = nan(73,W*2+1);
DiagT = nan(73,W);
for pairs = 1:length(PairOfCells)
cellid_a = PairOfCells{pairs,1}
cellid_d = PairOfCells{pairs,2}

[psth_a,spsth_a,~,tags,raszt_a] =  ultimate_psth(cellid_a,'trial',trigger,wind,'parts',partition);
[psth_d,spsth_d,~,~,raszt_d] =  ultimate_psth(cellid_d,'trial',trigger,wind,'parts',partition);


if iscell(raszt_a)
    try
    neuron1_raster = raszt_a{find(cellfun(@(x) x(end) == '1', tags))};
    neuron2_raster = raszt_d{find(cellfun(@(x) x(end) == '1', tags))};
    catch
        continue
    end
else
    %if size(raszt_a,1)<50
    %    continue
    %end
    S(pairs) = size(raszt_a,1);
    neuron1_raster = raszt_a;
    neuron2_raster = raszt_d;
end

% Bin the rasters

nTrials = size(neuron1_raster,1);
nBins = size(neuron1_raster,2);
% Initialize joint PSTH matrix
jointMatrix = zeros(nBins, nBins);

for ti = 1:nTrials
    st1 = find(neuron1_raster(ti,:) > 0);
    st2 = find(neuron2_raster(ti,:) > 0);
    for i = 1:length(st1)
        for j = 1:length(st2)
            jointMatrix(st1(i),st2(j)) = jointMatrix(st1(i),st2(j))+1;
        end
    end
end
jointMatrix = jointMatrix / nTrials;
A = mean(neuron1_raster, 1);
B = mean(neuron2_raster, 1);
AB = A'*B;
epsilon = 1e-7;
%sigmaA = sqrt(A .* (1 - A) + epsilon);
%sigmaB = sqrt(B .* (1 - B) + epsilon);
sigmaA = std(neuron1_raster,[], 1) + epsilon;
sigmaB = std(neuron2_raster,[], 1) + epsilon;
stdAB = sigmaA' * sigmaB; 
Interaction = (jointMatrix - AB)./stdAB;
% Normalize by number of trials to get probabilities




[nBinsOld, ~] = size(jointMatrix);
nBinsNew = floor(nBinsOld / binSize);
jointMatrixBinned = zeros(nBinsNew, nBinsNew);
InteractionBinned = zeros(nBinsNew, nBinsNew);
for i = 1:nBinsNew
    for j = 1:nBinsNew
        rows = (i-1)*binSize + (1:binSize);
        cols = (j-1)*binSize + (1:binSize);
        jointMatrixBinned(i,j) = mean(mean(jointMatrix(rows,cols)));
        InteractionBinned(i,j) = mean(mean(Interaction(rows,cols))); %/ (binSize)^2;
        %A_binned(i) = mean(A(rows));
        %B_binned(j) = mean(B(cols));
    end
end
%A = sum(jointMatrixBinned,2);
%B = sum(jointMatrixBinned,1);
%AB = A*B;
%SE = sqrt(AB .* (1 - AB) / nTrials);
%epsilon = 1e-7;
%AB_binned = A_binned'*B_binned;
jointMatrixT(pairs,:,:) = jointMatrixBinned;
InteractionT(pairs,:,:) = InteractionBinned; %- AB_binned;%./ (SE + epsilon);%./(B*A);

if isplotex
figure
subplot(2,3,5)
imagesc(1:nBinsNew, 1:nBinsNew, squeeze((InteractionBinned)));
axis xy;
xlabel('Neuron 2 time (ms)');
ylabel('Neuron 1 time (ms)');
title('Joint PSTH Matrix');


ax = subplot(2,3,4)
barh(1:nBins,sum(neuron1_raster,1))
set(ax, 'YAxisLocation', 'right');  % move vertical axis to the right
set(ax, 'XDir', 'reverse'); 

subplot(2,3,2)
bar(1:nBins,sum(neuron2_raster,1))
end
%for b = 1:nBinsNew
    % get original bin indices along this diagonal
    %diagIndices = (b-1)*binSize + (1:binSize);  %% likely not good!!!!!
    %diagValsBinned(b) = mean(arrayfun(@(k) sum(diag(jointMatrix, k)), diagIndices));
%    diagValsBinned(b) = diag(jointMatrixBinned);
%end
DiagT(pairs,:) = diag(jointMatrixBinned);

if isplotex
ax = subplot(2,3,6)
bar(1:nBinsNew,DiagT(pairs,:))
end

maxLag = nBins-1; % bins
lags = -maxLag:maxLag;
%diagValsT(neurons,:) = arrayfun(@(k) mean(diag(jointMatrix, k)), lags);
fineCCG = arrayfun(@(k) mean(diag(jointMatrix, k)), lags);   % fine-resolution CCG
fineInteractionCCG = arrayfun(@(k) mean(diag(Interaction, k)), lags);   % fine-resolution CCG
nFine = numel(fineCCG);
%binSize = 10;
nBinsNew = floor((maxLag + 1) / binSize);
lagsNew = (-nBinsNew : nBinsNew) * binSize;
ccgBinned = zeros(1, numel(lagsNew));
InteractionCCGBinned = zeros(1, numel(lagsNew));

for b = 1:numel(lagsNew)
    % offset range in fine bins
    k = lagsNew(b);
    % indices in fine CCG corresponding to this coarse lag bin
    idxStart = maxLag + 1 + k - floor(binSize/2);
    idxEnd   = idxStart + binSize - 1;
    idxStart = max(1, idxStart);
    idxEnd   = min(nFine, idxEnd);
    ccgBinned(b) = mean(fineCCG(idxStart:idxEnd));
    InteractionCCGBinned(b) = mean(fineInteractionCCG(idxStart:idxEnd));
end
ccgT(pairs, :) = ccgBinned;
InteractionCCGT(pairs, :) = InteractionCCGBinned;
if isplotex
ax = subplot(2,3,3)
bar(lagsNew, ccgBinned);
xlabel('Lag (bins)');
ylabel('Mean joint probability');
title('JPSTH lag profile');
end

end

% Average
figure
for i = 1:2

ax1 = subplot(2,4,1+(i-1)*4)
timeVec = (0:nBinsNew-1)*binSize;
imagesc(timeVec, timeVec, squeeze((nanmean(jointMatrixT(ie_pairs==i-1,:,:))))*1000);
colormap(ax1,hot)

hold on
axis xy;
line([timeVec(1), timeVec(end)],[timeVec(1), timeVec(end)],'Color','w')
xlabel('Dopaminergic neurons spike time from cue (ms)');
ylabel('Cholinergic neurons spike time from cue (ms)');
setmyplot_balazs
    
subplot(2,4,2+(i-1)*4)
bar(lagsNew,nanmean(ccgT(ie_pairs==i-1,:)./(nanmean(ccgT(ie_pairs==i-1,:),2)+0.0001))*1000)

setmyplot_balazs
xlabel('Lags (ms)');
ylabel('Count');
xlim([lagsNew(round(nBinsNew/2)),lagsNew(end-round(nBinsNew/2))+1])
 
    
ax2 = subplot(2,4,3+(i-1)*4)
timeVec = (0:nBinsNew-1)*binSize;
imagesc(timeVec, timeVec, squeeze((nanmean(InteractionT(ie_pairs==i-1,:,:)))));
colormap(ax2,turbo)
hold on
axis xy;
c = max(abs(nanmean(InteractionT(ie_pairs==i-1,:,:))),[],'all');
caxis([-c c]/1.2);
line([timeVec(1), timeVec(end)],[timeVec(1), timeVec(end)],'Color','w')
xlabel('Dopaminergic neurons spike time from cue (ms)');
ylabel('Cholinergic neurons spike time from cue (ms)');
setmyplot_balazs
% ax = subplot(2,3,2+(i-1)*3)
% bar(timeVec,mean(DiagT(ie_pairs==i-1,:)))
% hold on
% xlabel('Time from cue (ms)');
% ylabel('Coincidence frequency')

%plot(1:nBinsNew,mean(DiagT(ie_pairs==i-1,:)))
%setmyplot_balazs
ax = subplot(2,4,4+(i-1)*4)
%bar(lagsNew,nanmean(InteractionCCGT(ie_pairs==i-1,:)))
middleinx  = ceil(length(lagsNew)/2);
%boxplot(InteractionCCGT(ie_pairs==i-1,middleinx-3:middleinx+3), 'Positions', lagsNew(middleinx-3:middleinx+3),'Symbol','')
INT = InteractionCCGT(ie_pairs==i-1,middleinx-4:middleinx+4);
errorshade((-4:4)*binSize,nanmean(INT),nanstd(INT)./sqrt(sum(ie_pairs==i-1)))
%ax = gca;
%ax.XAxisLocation = 'origin';
%ax.YAxisLocation = 'origin';
ylim([-max(abs(nanmean(INT))),max(abs(nanmean(INT)))])
box off
hold on
%for bin = 1:length(lagsNew)
%    sigstar({[lagsNew(bin),lagsNew(bin)]},signrank(InteractionCCGT(ie_pairs==i-1,bin)))    
%end

[~,MI] = max(abs(nanmean(INT)));
sigstar({[([MI,MI]-5)*binSize]},signrank(INT(:,MI)))
%sigstar({[MI,MI]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx-5+MI-1)))
%sigstar({[-1,-1]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx-1:middleinx-1),2)))
%sigstar({[1,1]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx+1:middleinx+1),2)))
%sigstar({[-2,-2]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx-2:middleinx-2),2)))
%sigstar({[2,2]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx+2:middleinx+2),2)))
%sigstar({[-3,-3]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx-3:middleinx-3),2)))
%sigstar({[3,3]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx+3:middleinx+3),2)))
%sigstar({[-4,-4]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx-4:middleinx-4),2)))
%sigstar({[4,4]*binSize},signrank(mean(InteractionCCGT(ie_pairs==i-1,middleinx+4:middleinx+4),2)))


setmyplot_balazs
xlabel('Lags (ms)');
ylabel('Count');
end






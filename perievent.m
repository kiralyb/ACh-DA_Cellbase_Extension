    
function [eta,sem,spmat] = perievent(discE,data,sr,wns,needn,dt,sigma)
%ASTANORMMOD    Normalized Spike Triggered Average.
%   [MEAN SEM MATRIX] = ASTANORM(DISCE,DATA,SR,WN) calculates STA of EEG
%   and discriminated unit. Input arguments:
%       DISCE: discriminative events (in data points)
%       DATA: contiuous photometry signal
%       SR: sampling rate
%       WNS: window size (in seconds)
%       NEEDN: normalization factor
%    Output arguments:
%       ETA: normalized event triggered average
%       SEM: standard error of the mean
%       MATRIX: ETA stack for ETA image
%
%   See also ASTANORM2.

% Time window in sampling points
wns=wns*sr;
wn = wns(2) - wns(1);

% Standardize EEG
if needn == 1
    eeg = (data - mean(data)) / std(data);
else
    eeg = data;
end

% Calculate ETA
discE = discE(discE+wns(1)>0&discE+wns(2)<=length(eeg));
lenv = length(discE);
eeginx = repmat(discE(:)+wns(1),1,wn+1) + repmat(0:wn,lenv,1);
matrix = eeg(eeginx);
if lenv==1
    matrix=matrix';
end
% Removing outlayers and replacing them with NaN's
%[AA TF] = rmoutliers(matrix,'ThresholdFactor',20);
%matrix(TF == 1,:) = NaN;

for i = 1:size(matrix,1)
spmat(i,:) = smoothed_psth(matrix(i,:),dt,sigma);
end

%Normalising baselines to mean 0
%newmat = [];
% for i = 1:size(spmat,1)
%     baseline = mean(spmat(i,1:(abs(wns(1)))));
%     newmat = [newmat;(spmat(i,1:end)-baseline)];
% end 


eta = nanmean(spmat,1);
sem = nanstd(spmat,[],1) ./ sqrt(size(matrix,1));

%sem = smoothed_psth(nanstd(spmat,[],1) ./ sqrt((size(matrix,1))),dt,sigma);

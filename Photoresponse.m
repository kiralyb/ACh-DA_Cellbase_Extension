function [areasUnderHalfProminence,careasUnderHalfProminence,peakLocation] = Photoresponse(PhotoMatrix,earlierpeakfactor,peakfilterfactor,PeakorTrough,sr,win)

if nargin < 2
    earlierpeakfactor = 1;    
end
if nargin < 3
    peakfilterfactor = 0;
end
if nargin < 4
    PeakorTrough = 1; % 1 or -1
end
if nargin <5
    sr = 1;
end
if nargin <6
   win = [1,length(PhotoMatrix)]; 
end

for i= 1:size(PhotoMatrix,1)

signal = PhotoMatrix(i,:)*PeakorTrough;

    
    
    
% Find peaks with prominence
[peak, peak_locs, peak_width, peak_prom] = findpeaks(signal,'MinPeakWidth',1);
peak(peak_locs<win(1) | peak_locs>win(2)) = [];
peak_width(peak_locs<win(1) | peak_locs>win(2)) = [];
peak_prom(peak_locs<win(1) | peak_locs>win(2)) = [];
peak_locs(peak_locs<win(1) | peak_locs>win(2)) = [];
[trough, trough_locs, through_width, trough_prom] = findpeaks(-signal,'MinPeakWidth',1);


% Find the peak with the highest prominence
[highestprom, highestProminenceIndex] = max(peak_prom);
earlierpeak = find(peak_prom(1:highestProminenceIndex-1)>highestprom/earlierpeakfactor,1,'first')

if ~isempty(earlierpeak)
    highestProminenceIndex = earlierpeak(1);
end
peakValue = peak(highestProminenceIndex);
if ~isempty(highestProminenceIndex) 
peakLocation(i) = peak_locs(highestProminenceIndex);
else
areasUnderHalfProminence(i) = nan;
careasUnderHalfProminence(i) = 0;
peakLocation(i) =nan;
continue
end

if (sum(trough_prom(trough_locs < peak_locs(highestProminenceIndex)).*through_width(trough_locs<peak_locs(highestProminenceIndex))>peakValue*peak_width(highestProminenceIndex))>0) & PeakorTrough == -1
[areasUnderHalfProminence,careasUnderHalfProminence,peakLocation] = Photoresponse(PhotoMatrix,earlierpeakfactor,peakfilterfactor,1,sr,win);
%areasUnderHalfProminence(i) = nan;
%careasUnderHalfProminence(i) = nan;
%peakLocation(i) =nan;
continue
end


% Compute the half-prominence value for the highest prominence peak
halfProminence = peak_prom(highestProminenceIndex) / 2;

% Find the boundaries of the half-prominence
leftBoundary = find(signal(1:peakLocation(i)) < (peakValue - halfProminence), 1, 'last');
rightBoundary = peakLocation(i) + find(signal(peakLocation(i):end) < (peakValue - halfProminence), 1, 'first') - 1;

% Check if leftBoundary is empty and set it to 1 if needed
if isempty(leftBoundary)
    leftBoundary = 1;
end
%base = min(signal(1:leftBoundary));
%base = mean(signal(1));
base = signal(leftBoundary);

% Check if rightBoundary is empty and set it to the last index if needed
if isempty(rightBoundary)
    rightBoundary = numel(signal);
end

% Calculate the area under the highest prominence half-peak
areasUnderHalfProminence(i) = sum(signal(leftBoundary:rightBoundary))*PeakorTrough/sr;
careasUnderHalfProminence(i) = sum(signal(leftBoundary:rightBoundary)-base)*PeakorTrough/sr;  
if careasUnderHalfProminence(i) < (max(signal(1:leftBoundary))-base)*leftBoundary*peakfilterfactor
    peakLocation(i) =nan;
end

end
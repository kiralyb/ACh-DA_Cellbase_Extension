function [CCR_smooth,CCR_norm,PairOfCells]  = CCG_normalizer(path)
load(path);
if size(CCR,2) > 1001
     inx = (size(CCR,2)-1)/2 - 500 : (size(CCR,2)-1)/2 + 500;
else
    inx = 1:size(CCR,2);
end

CCR_norm2 = bsxfun(@rdivide,CCR(:,inx),sum(CCR,2));
CCR_Poi_norm = bsxfun(@rdivide,CCR_Poi(:,inx),sum(CCR_Poi,2));
CCR_norm = (CCR_norm2 - CCR_Poi_norm)*size(CCR,2);
CCR_smooth = nan(size(CCR_norm,1),length(inx));
for ind = 1:size(CCR_norm,1)
    if sum(CCR(ind,inx),2) > 100
    CCR_smooth(ind,:) = smoothdata(CCR_norm(ind,:),2,'gaussian',(size(CCR,2)/sqrt(SegmentLength(ind))));   
    end
end
%CCR_smooth = CCR_smooth * size(CCR_smooth,2);
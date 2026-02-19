function [CrossPairs,ie_pairs] = getpairs(IDS1,IDS2,subID)

load('L:\_CellBases\ACh_DA_Cellbase\CCG_new\all\stim_light_filter2\CCG_matrices.mat')
crossinx=[];
ie_pairs = [];
for i = 1: length(PairOfCells)
    a1 = ismember(findcellpos(PairOfCells(i,1)),IDS1);
    a2 = ismember(findcellpos(PairOfCells(i,2)),IDS1);
    d2 = ismember(findcellpos(PairOfCells(i,2)),IDS2);
    d1 = ismember(findcellpos(PairOfCells(i,1)),IDS2);
    if ((a1 & d2) | (d1 & a2))
       crossinx = [crossinx,i];
       ie_pairs = [ie_pairs, ismember(findcellpos(PairOfCells(i,2)),subID)];
    end
end
CrossPairs = PairOfCells(crossinx,:);
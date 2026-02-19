function [crosspairs,ie_pairs] = pairIDfinder(i_ach,i_da,i_da_i)
loadcb
RATID_a = getvalue('RatId',CELLIDLIST(i_ach))
RATID_d = getvalue('RatId',CELLIDLIST(i_da))
SESSIONID_A = getvalue('DateNum',CELLIDLIST(i_ach))
SESSIONID_D = getvalue('DateNum',CELLIDLIST(i_da))
crosspairs = cell(0,2)
ie_pairs = []
for i = 1:length(RATID_a)
    inx = find(RATID_a(i) == RATID_d & SESSIONID_A(i) == SESSIONID_D)
    for j = 1:length(inx)
        if CELLIDLIST{i_ach(i)}(find(CELLIDLIST{i_ach(i)}=='_',1,'last')-1) == CELLIDLIST{i_da(inx(j))}(find(CELLIDLIST{i_da(inx(j))}=='_',1,'last')-1)
        crosspairs(end+1,:)  = {CELLIDLIST{i_ach(i)},CELLIDLIST{i_da(inx(j))}}
        if ismember(i_da(inx(j)),i_da_i)
            ie_pairs = [ie_pairs,1]; 
        else
            ie_pairs = [ie_pairs,0];
        end
        end
    end
end
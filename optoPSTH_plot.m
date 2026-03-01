function p = optoPSTH_plot(clustnums,area,stimmedAs,colors,legends)
figure
hold on
loadcb
for i = 1:length(clustnums)
    if area(i) == 1
        Clust = getvalue('HDB_Cluster_num');
    else
        Clust = getvalue('VTA_Cluster_num');
    end
    if clustnums(i) > 5
        inx = find(abs(Clust) > 4);
    else
        inx = find(abs(Clust)==clustnums(i));
    end
    mask = ismember(getvalue('RatID',CELLIDLIST(inx)),stimmedAs);
    p(i) = optoPSTH(CELLIDLIST(inx(mask)),'stim',colors{i},i-1);
end
legend(legends)
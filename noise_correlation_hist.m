function noise_correlation_hist(crosspairs,ie_pairs)

[PP,RR] = noise_correlation(crosspairs,1,0);
pp_edges = 0:0.05:1;
rr_edges = -0.5:0.1:0.5;
for row = 1:2
    PH(row,:) = histcounts(PP(ie_pairs == row-1), pp_edges);
    RH(2*row-1,:) = histcounts(RR((ie_pairs == row-1) & PP < 0.05), rr_edges);
    RH(2*row,:) = histcounts(RR(ie_pairs == row-1), rr_edges) - RH(2*row-1,:);
end
figure
subplot(1,2,1)
bar(pp_edges(1:end-1) + diff(pp_edges)/2,PH','stacked')
legend({'a - d_i','a - d_e'})
subplot(1,2,2)
bar(rr_edges(1:end-1) + diff(rr_edges)/2,RH','stacked')
legend({'a - d_e','a - d_e not sig.','a - d_i','a - d_i not sig.'})
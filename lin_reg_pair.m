function [Ba,Bd,Bc] = lin_reg_pair(pairids)
cellid1 = pairids{1};
cellid2 = pairids{2};
Ba = nan;
Bd = nan;
Bc = nan;

VE = loadcb(cellid1,'TrialEvents');   % load events
newcues = VE.Soundtype > 2;
if sum(newcues)==0
    'no new cue trials';
    return
end

resps = ~isnan(VE.Type);
performance_bin = resps(newcues);

if var(performance_bin) < 0.01
    'no variance to explain in behaviour';
    return
end

ach_b = get_norm_spikeresp(cellid1,newcues);
da_b = get_norm_spikeresp(cellid2,newcues);

mdl = fitglm(ach_b, performance_bin', ...
    'Distribution', 'binomial');
[~, ~, ~, Ba] = perfcurve(performance_bin', predict(mdl, ach_b), 1);

mdl_d = fitglm(da_b, performance_bin', ...
    'Distribution', 'binomial');
[~, ~, ~, Bd] = perfcurve(performance_bin', predict(mdl_d, da_b), 1);

X = [ach_b,da_b];
mdl_c = fitglm(X, performance_bin', ...
    'Distribution', 'binomial');
[~, ~, ~, Bc] = perfcurve(performance_bin', predict(mdl_c, X), 1);


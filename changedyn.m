function [R_ach, P_ach, R_ach_b, P_ach_b, R_da, P_da, R_da_b, P_da_b ,linparams_ach,R2_ach,p_ach,linparams_da,R2_da,p_da,R_ach_d,P_ach_d,R_da_d,P_da_d,R_par,P_par,pe_var,B,B_p,B_pC,Bd,Bd_p,Bc,Bc_p,trialnum,linparams_cross,linparams_cross_rev,R2_cross,R2_cross_rev] = changedyn(cellids,isplot)
if nargin<2
    isplot = 1;
end


if length(cellids) == 2
    cellid1 = cellids{1};
    cellid2 = cellids{2};
    numcells = 2;

else
    cellid1 = cellids{1};
    numcells = 1;
end
    linparams_ach = nan(1,3);
    R2_ach = nan;
    p_ach = nan;
    pred_ach = nan;
    R_da = nan;
    P_da = nan;
    R_da_d = nan;
    P_da_d = nan;
    R_da_b = nan;
    P_da_b = nan;
    linparams_da = nan;
    R2_da = nan;
    p_da = nan;
    R_par = nan;
    P_par = nan;
    B_pC = nan;
    Bd = nan;
    Bd_p = nan;
    Bc = nan;
    Bc_p = nan;
    linparams_cross = nan(1,3);
    linparams_cross_rev = nan(1,3);
    R2_cross = nan;
    R2_cross_rev = nan;
    
R_ach = nan;
P_ach = nan;
R_ach_b = nan;
P_ach_b = nan;
R_ach_d = nan;
P_ach_d = nan;
B = nan;
B_p = nan;
pe_var = nan;
trialnum = nan;

sigma = 0.02;
wind = [-0.0,0.4];
smoother = 7;
meanwindow = [smoother,smoother];
meanwindow_s = [smoother,smoother];
smoothmet = "gaussian";

VE = loadcb(cellid1,'TrialEvents');   % load events
newcues = VE.Soundtype > 2;
if sum(newcues)==0
    'no new cue trials'
    return
end
resps = ~isnan(VE.Type);
performance_bin = resps(newcues);
performance = smoothdata(performance_bin,smoothmet,meanwindow);
pe_var = var(performance);
trialnum = sum(newcues);
if pe_var < 0.005
    'no variance to explain in behaviour'
    return
end

[~,spsth_a,~,tags,raszt_a] =  ultimate_psth(cellid1,'trial','StimulusOn',wind,'parts','#TrialCode','sigma',sigma);
if numcells == 2
[~,spsth_d,~,~,raszt_d] =  ultimate_psth(cellid2,'trial','StimulusOn',wind,'parts','#TrialCode','sigma',sigma);
end
for i = 1:length(tags)
    inxs(i) = str2double(tags{i}(end));
end
inx_4 = find(inxs==1);
inx_n = find(inxs==3);
[~,locs,width,prominance] = findpeaks(spsth_a(inx_4(1),:));
[~,ind] = max(prominance);
I_a_4Hz = locs(ind);
window_a_4Hz = ceil(width(ind));
[~,locs,width,prominance] = findpeaks(spsth_a(inx_n(1),:));
[~,ind] = max(prominance);
I_a_n = locs(ind);
window_a_n = ceil(width(ind));
if numcells == 2
    [~,locs,width,prominance] = findpeaks(spsth_d(inx_4(1),:));
    [~,ind] = max(prominance);
    I_d_4Hz = locs(ind);
    window_d_4Hz = ceil(width(ind));
    [~,locs,width,prominance] = findpeaks(spsth_d(inx_n(1),:));
    [~,ind] = max(prominance);
    I_d_n = locs(ind);
    window_d_n = ceil(width(ind));
    
end

RESP1_4Hz = sum(raszt_a{inx_4(1)}(:,max([1,I_a_n-window_a_n]):min([length(spsth_a),I_a_n+window_a_n])),2);
[psth_a,spsth_a,~,~,raszt_a] =  ultimate_psth(cellid1,'trial','StimulusOn',wind,'parts','All','sigma',sigma);
RESP1 = sum(raszt_a(newcues,max([1,I_a_n-window_a_n]):min([length(spsth_a),I_a_n+window_a_n])),2);
if numcells == 2
    RESP2_4Hz = sum(raszt_d{inx_4(1)}(:,max([1,I_d_n-window_d_n]):min([length(spsth_d),I_d_n+window_d_n])),2);
    [psth_d,spsth_d,~,~,raszt_d] =  ultimate_psth(cellid2,'trial','StimulusOn',wind,'parts','All','sigma',sigma);
    RESP2 = sum(raszt_d(newcues,max([1,I_d_n-window_d_n]):min([length(spsth_d),I_d_n+window_d_n])),2);
end

ach_b = RESP1/mean(RESP1_4Hz)';
ach_s = smoothdata(ach_b,smoothmet,meanwindow_s)';
ach = ach_s';
pe = performance';

if numcells == 2
da_b = RESP2/mean(RESP2_4Hz)';
da_s = smoothdata(da_b,smoothmet,meanwindow_s)';
da = da_s';
end



[R_ach, P_ach] = corr(ach, pe,'type','Spearman'); 
[R_ach_b, P_ach_b] = corr(ach_b, double(performance_bin)','type','Spearman'); 
[R_ach_d, P_ach_d] = corr(diff(ach),diff(pe),'type','Spearman');


mdl = fitglm(ach_b, performance_bin', ...
'Distribution', 'binomial');
[~, ~, ~, B] = perfcurve(performance_bin', predict(mdl, ach_b), 1);
B_p = mdl.Coefficients.pValue(2);


if numcells == 2
[R_da, P_da] = corr(da, pe,'type','Spearman');
[R_da_b, P_da_b] = corr(da_b, double(performance_bin)','type','Spearman');
[R_da_d, P_da_d] = corr(diff(da),diff(pe),'type','Spearman'); 

mdl_d = fitglm(da_b, performance_bin', ...
'Distribution', 'binomial');
[~, ~, ~, Bd] = perfcurve(performance_bin', predict(mdl_d, da_b), 1);
Bd_p = mdl_d.Coefficients.pValue(2);

X = [ach_b,da_b];
mdl_c = fitglm(X, performance_bin', ...
'Distribution', 'binomial');
[~, ~, ~, Bc] = perfcurve(performance_bin', predict(mdl_c, X), 1);
Bc_p = mdl_c.Coefficients.Estimate;

X = [ach_b,da_b,ach_b.*da_b];
mdl_i = fitglm(X, performance_bin', ...
'Distribution', 'binomial');  
B_pC(1) = LR_test(mdl_c,mdl);
B_pC(2) = LR_test(mdl_c, mdl_d);
B_pC(3) = LR_test(mdl_i, mdl_c);
%Bc_p = mdl_i.Coefficients.Estimate;
%Bc_p = mdl_i.Coefficients.pValue;
%Bc = 1 - mdl_i.Deviance / mdl_null.Deviance;

[R_par, P_par] = corr([ach,da,pe],'type','Spearman');
[R_par_c, P_par_c] = partialcorr([ach,da,pe],'type','Spearman');
R_par(2,1) = R_par_c(2,1);
R_par(3,1) = R_par_c(3,1);
R_par(3,2) = R_par_c(3,2);
P_par(2,1) = P_par_c(2,1);
P_par(3,1) = P_par_c(3,1);
P_par(3,2) = P_par_c(3,2);

end

if sum(newcues) >= 100
[linparams_ach,R2_ach,p_ach,pred_ach] = deltatlinmodfit(ach_s,pe',10);%-mean(ach_s),pe'-mean(ach_s),10);
if numcells == 2
[linparams_da,R2_da,p_da,pred_da] = deltatlinmodfit(da_s,pe',10);%-mean(da_s),pe'-mean(da_s),10);
[linparams_cross,R2_cross,p_cross,~] = deltatlinmodfit(ach_s-mean(ach_s),da_s-mean(da_s),10);%-mean(da_s),pe'-mean(da_s),10);
[linparams_cross_rev,R2_cross_rev,p_corss_rev,~] = deltatlinmodfit(da_s-mean(da_s),ach_s-mean(ach_s),10);%-mean(da_s),pe'-mean(da_s),10);
else
    
end
end

if isplot 
figure
subplot(1,2,1)
hold on
plot(pe,'r')
plot(ach_s,'c')
xlabel('Trials')
ylabel('responsivness')
borders = find(diff(VE.thrdSound(newcues))~=0);
hold on
if ~isempty(borders)
line([borders;borders],[zeros(length(borders),1),ones(length(borders),1)]')
end
scatter(1:length(performance_bin),performance_bin,'|')
setmyplot_balazs

if numcells == 2
plot(da_s,'m')
end
subplot(1,2,2)
hold on
plot(pred_ach,'c')
xlabel('Trials')
ylabel('responsivness')
if numcells == 2
plot(pred_da,'m')
end
plot(pe,'r')


pe_var = var(pe);
trialnum = length(pe);

if numcells == 2
title([cellid1,cellid2]);
else
title(cellid1); 
end

end
% ach = smoothdata(ach_b,smoothmet,meanwindow_s.*[1,0]);
% if numcells == 2
% da = smoothdata(da_b,smoothmet,meanwindow_s.*[1,0]);
% end
% pe = smoothdata(performance_bin,smoothmet,meanwindow_s.*[1,0])';
% [h.ap,GRANGER_pvalue.ap,GRANGER_F.ap] = gctest(ach,pe,'NumLags',numlags);
% 
% [h.pa,GRANGER_pvalue.pa,GRANGER_F.pa,cvalue] = gctest(pe,ach,'NumLags',numlags);
% if numcells == 2
% [h.dp,GRANGER_pvalue.dp,GRANGER_F.dp] = gctest(da,pe,'NumLags',numlags);
% 
% [h.pd,GRANGER_pvalue.pd,GRANGER_F.pd,cvalue] = gctest(pe,da,'NumLags',numlags);
% [h.da,GRANGER_pvalue.da,GRANGER_F.da,~] = gctest(da,ach,'NumLags',numlags);
% [h.ad,GRANGER_pvalue.ad,GRANGER_F.ad,cvalue] = gctest(ach,da,'NumLags',numlags);
% 
% 
% 
% [h,GRANGER_pvalue.apd,GRANGER_F.apd,cvalue] = gctest(ach,pe,da,'NumLags',numlags);
% [h,GRANGER_pvalue.dpa,GRANGER_F.dpa,cvalue] = gctest(da,pe,ach,'NumLags',numlags);
% [h,GRANGER_pvalue.adp,GRANGER_F.adp,cvalue] = gctest(ach,da,pe,'NumLags',numlags);
% [h,GRANGER_pvalue.dap,GRANGER_F.dap,cvalue] = gctest(da,ach,pe,'NumLags',numlags);
% [h,GRANGER_pvalue.pda,GRANGER_F.pda,cvalue] = gctest(pe,da,ach,'NumLags',numlags);
% [h,GRANGER_pvalue.pad,GRANGER_F.pad,cvalue] = gctest(pe,ach,da,'NumLags',numlags);
%end


function [estimated_parameters, R2, p, predy] = deltatlinmodfit(x, y, max_shift)

t = 1:length(y);

% Define model with fixed integer shift
model_with_fixed_shift = @(p1p3, shift, t) ...
    abs(p1p3(1)) * interp1(t, x, t + shift, 'linear', 'extrap') + p1p3(2);

% Initial guess for p1 (gain) and p3 (bias)
initial_guess = [1, 0];

% Range of integer shifts to try
shift_range = -max_shift:max_shift;

best_SSR = Inf;
best_params = [];

for shift = shift_range
    model = @(p1p3, t) model_with_fixed_shift(p1p3, shift, t);
    options = optimset('Display','off');
    
    % Fit p1 and p3 with fixed integer shift
    [p1p3_fit, resnorm] = lsqcurvefit(model, initial_guess, t, y, [], [], options);
    
    if resnorm < best_SSR

        best_SSR = resnorm;
        best_params = [p1p3_fit(1), shift, p1p3_fit(2)];
    end
end

estimated_parameters = best_params;
predy = abs(estimated_parameters(1)) * interp1(t, x, t + estimated_parameters(2), 'linear', 'extrap') + estimated_parameters(3);

% Compute R² and p-value
SSR = sum((y - predy).^2);
SST = sum((y - mean(y)).^2);
R2 = 1 - (SSR / SST);

n = length(y);
k = 3;
df = n - k;
F_stat = (R2 / k) / ((1 - R2) / df);
p  = 1 - fcdf(F_stat, k, df);

function p = LR_test(mdl1, mdl0)
    D1 = mdl1.Deviance;
    D0 = mdl0.Deviance;
    df = mdl0.DFE - mdl1.DFE;
    chi2 = D0 - D1;
    p = 1 - chi2cdf(chi2, df);
    

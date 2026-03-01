function [R, R_d, linparams, R2, R_par, linparams_cross, linparams_cross_rev] = changedyn2(cellids,isplot)
if nargin < 2
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
R = nan;
R_d = nan;
R_par = nan;
linparams = nan(1,3);
R2 = nan;
linparams_cross = nan(1,3);
linparams_cross_rev = nan(1,3);

smoother = 7;
meanwindow = [smoother,smoother];
smoothmet = "gaussian";

VE = loadcb(cellid1,'TrialEvents');   % load events
newcues = VE.Soundtype > 2;
if sum(newcues)==0
    'no new cue trials';
    return
end
resps = ~isnan(VE.Type);
performance_bin = resps(newcues);
performance = smoothdata(performance_bin,smoothmet,meanwindow);
pe_var = var(performance);

if pe_var < 0.005
    'no variance to explain in behaviour';
    return
end

ach_b = get_norm_spikeresp(cellid1,newcues);
ach = smoothdata(ach_b,smoothmet,meanwindow);
pe = performance';

if numcells == 2
da_b = get_norm_spikeresp(cellid2,newcues);
da = smoothdata(da_b,smoothmet,meanwindow);
end

R = corr(ach, pe,'type','Spearman'); 
R_d = corr(diff(ach),diff(pe),'type','Spearman');

if numcells == 2    
    R_par = corr([ach,da,pe],'type','Spearman');
    R_par_c = partialcorr([ach,da,pe],'type','Spearman');
    R_par(2,1) = R_par_c(2,1);
    R_par(3,1) = R_par_c(3,1);
    R_par(3,2) = R_par_c(3,2);
end

if sum(newcues) >= 100
    [linparams,R2,~,pred_ach] = deltatlinmodfit(ach,pe',10);
    if numcells == 2
        linparams_cross = deltatlinmodfit(ach-mean(ach),da'-mean(da),10);
        linparams_cross_rev = deltatlinmodfit(da-mean(da),ach'-mean(ach),10);
    end
end

if isplot 
figure
subplot(1,2,1)
hold on
plot(pe,'r')
plot(ach,'c')
if numcells == 2
   plot(da,'m')
end
xlabel('Trials')
ylabel('responsivness')
borders = find(diff(VE.thrdSound(newcues))~=0);
hold on
if ~isempty(borders)
line([borders;borders],[zeros(length(borders),1),ones(length(borders),1)]')
end
scatter(1:length(performance_bin),performance_bin,'|')
setmyplot_balazs

subplot(1,2,2)
hold on
plot(pred_ach,'c')
xlabel('Trials')
ylabel('responsivness')
plot(pe,'r')

if numcells == 2
title([cellid1,cellid2]);
else
title(cellid1); 
end

end

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

    

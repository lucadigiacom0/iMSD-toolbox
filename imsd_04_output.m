%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = optimoptions('lsqcurvefit','Display','off');

if exist('fitt2comp','var') && fitt2comp
%%% ==================== TWO-COMPONENT OUTPUT ====================
tvv = time_vec(:);
rmax = max(abs(a1));

%%% reliable mobile range: the fast MOBILE component saturates in the analysis
%%% window at long lag (sigma_D approaches the crop half-size). Keep only the
%%% initial part where sigma_D^2 is still well inside the window.
i_rel = find(s_G > (0.7*rmax)^2, 1) - 1;
if isempty(i_rel) || i_rel < 5, i_rel = numel(s_G); end
rel = (1:i_rel)';
xr = tvv(rel);  yr = s_G(rel);

%%% ----- MOBILE curve = FREE-DIFFUSION component of the decomposition -----
%%% Fit the LINEAR law  sigma_D^2(tau) = sigma0^2 + 4*D*tau  on the short
%%% reliable (linear) range -> D_mobile. 

nfit  = max(4, min(round(numel(tvv)*0.15), i_rel));
[plin,Slin] = polyfit(tvv(2:nfit), s_G(2:nfit), 1);
D_M   = plin(1)/4;
offs02= plin(2);
y_fit_lin = polyval(plin, tvv);
alpha = 1; L = NaN; tau_c = NaN; D_m = NaN; v = NaN; mob_type = 'diffusion';

%%% parameter errors (from the covariance of the linear fit) and R^2
Clin = (inv(Slin.R)*inv(Slin.R)') * Slin.normr^2 / max(Slin.df,1);
D_M_err    = sqrt(abs(Clin(1,1)))/4;
offs02_err = sqrt(abs(Clin(2,2)));
R2_mob = 1 - sum((yr - y_fit_lin(rel)).^2) / sum((yr - mean(yr)).^2);


% PSF_FWHM is provided in µm
% if ~exist('PSF_FWHM','var')
%     PSF_FWHM = 0.25;   % default FWHM (µm)
% end

PSF_var = (PSF_FWHM/2.355)^2;
if offs02 > PSF_var
    Size_nm = sqrt(offs02 - PSF_var)*1000;
else
    Size_nm = sqrt(max(offs02,0))*1000;
end

%%% ----- TRAPPED amplitude decay -> tau_T (Eq. A20) -----
tauD = sigma_T^2/(4*max(D_M,eps));
tlo  = max(3*tauD, tvv(max(2,round(numel(tvv)*0.1))));
mm   = tvv >= min(tlo, tvv(round(numel(tvv)/4)));
dec  = @(p,t) p(1)*exp(-t./p(2)) + p(3);
tauT = NaN; tauT_err = NaN; R2_dec = NaN; amp_fit = nan(numel(tvv),1);
try
    i1 = find(mm,1);
    [pT,~,rT,~,~,~,JT] = lsqcurvefit(dec, [A_trap(i1), tvv(end)/3, 0], tvv(mm), A_trap(mm), [0 1e-3 0], [Inf 50 Inf], opts);
    tauT = pT(2); amp_fit = dec(pT, tvv);
    try, ciT = nlparci(pT, rT, 'jacobian', JT); tauT_err = (ciT(2,2)-ciT(2,1))/2; catch, end
    R2_dec = 1 - sum((A_trap(mm)-dec(pT,tvv(mm))).^2) / sum((A_trap(mm)-mean(A_trap(mm))).^2);
catch
end

clc
fprintf('***********************************\niMSD output (two-component):\n\n');
fprintf('  tau_T (trapping time) = %.4f +/- %.4f s   (R^2 = %.3f)\n', tauT, tauT_err, R2_dec);
fprintf('  D_mobile              = %.4f +/- %.4f um^2/s   (R^2 = %.3f)\n', D_M, D_M_err, R2_mob);
fprintf('  sigma0^2 (intercept)  = %.4f +/- %.4f um^2\n', offs02, offs02_err);
fprintf('  Size (deconvolved)    = %.1f nm   (PSF_waist = %.0f nm)\n', Size_nm, PSF_waist*1000);
fprintf('  sigma_T (waist)       = %.4f um   (trapped iMSD = %.4f um^2, flat)\n', sigma_T, sigma_T^2);
fprintf('  tau_D = %.4f s  ->  tau_T/tau_D = %.1f\n', tauD, tauT/tauD);
if tauT < 3*tauD
    warning('tau_T is close to tau_D -> tau_T is UNRELIABLE (biased).');
end
fprintf('***********************************\n');

%%% ------------------- FIGURE "iMSD output" -------------------
figure('Name','iMSD output','NumberTitle','off');

%%% (1) trapped amplitude decay -> tau_T
subplot(2,2,1)
semilogy(tvv, A_trap, 'o', 'Color',[0.85 0.33 0.10], 'MarkerFaceColor',[0.85 0.33 0.10]); hold on
if ~isnan(tauT), plot(tvv, amp_fit, 'k-', 'LineWidth',1.5); end
hold off; grid on
xlabel('\tau (s)'); ylabel('trapped amplitude g_T(\tau)');
title(sprintf('Trapped amplitude decay  ->  \\tau_T = %.3f s', tauT));
set(gca,'fontsize',12)

%%% (2) the TWO iMSD curves (main panel)
ax = subplot(2,2,2);
p_tr = plot(tvv, s_G_trap, 'o-', 'Color',[0.85 0.33 0.10], 'MarkerFaceColor',[0.85 0.33 0.10], 'LineWidth',1.2); hold on
p_mo = plot(xr, yr, 's', 'Color',[0 0.45 0.74], 'MarkerFaceColor',[0 0.45 0.74]);
p_fi = plot(xr, y_fit_lin(rel), 'k--', 'LineWidth',1.5);
hold off; grid on
ylim([0 1.15*max([yr; s_G_trap])])
xlabel('\tau (s)'); ylabel('iMSD  \sigma^2 (\mum^2)'); title('iMSD plot: trapped (flat) + mobile');
legend([p_tr p_mo p_fi], {'trapped (flat)','mobile (data)','mobile linear fit'}, 'Location','northwest');
annT = sprintf(['MOBILE (free diffusion)\nD=%.3f\\pm%.3f \\mum^2/s (R^2=%.2f)\n\\sigma_0^2=%.3f \\mum^2\nSize=%.0f nm\n\\sigma_T=%.3f \\mum\n\\tau_T=%.3f\\pm%.3f s (R^2=%.2f)'], ...
    D_M, D_M_err, R2_mob, offs02, Size_nm, sigma_T, tauT, tauT_err, R2_dec);
text(0.98,0.05, annT, 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',9,'BackgroundColor','w','EdgeColor','k');
set(gca,'fontsize',12)

%%% (3) mobile amplitude
subplot(2,2,3)
plot(tvv, amp_G, '.', 'Color',[0 0.45 0.74]); grid on
xlabel('\tau (s)'); ylabel('mobile amplitude g_D(\tau)'); title('Mobile amplitude')
set(gca,'fontsize',12)

%%% (4) residuals of the mobile fit
subplot(2,2,4)
plot(xr, yr - y_fit_lin(rel), 'r.-', 'LineWidth',1); hold on
yline(0,'b-','LineWidth',1.2); hold off; grid on
xlabel('\tau (s)'); ylabel('mobile residual (\mum^2)'); title('Mobile fit residuals')
set(gca,'fontsize',12)

%%% ------------------- EXCEL -------------------
if export_value==1
    [out_dir,~,~] = fileparts(filename);
    export_filename = fullfile(out_dir, [filename '_export.xlsx']);
    if isfile(export_filename), delete(export_filename); end
    T_data = table(tvv, s_G_trap(:), s_G(:), s_G_err(:), y_fit_lin(:), ...
        A_trap(:), A_trap_err(:), amp_fit(:), amp_G(:), amp_G_err(:), R2_G(:), ...
        'VariableNames', {'time_s','iMSD_trapped_flat','iMSD_mobile','iMSD_mobile_err', ...
                          'iMSD_mobile_linfit','trapped_amplitude','trapped_amplitude_err', ...
                          'trapped_amp_fit','mobile_amplitude','mobile_amplitude_err','R2_perlag_fit'});
    writetable(T_data, export_filename, 'Sheet','iMSD_curves');
    Pn={'tau_T (s)';'tau_T_err (s)';'R2_trapped_decay';'sigma_T (um)'; ...
        'iMSD_trapped_flat (um^2)';'D_mobile (um^2/s)';'D_mobile_err (um^2/s)';'R2_mobile_fit'; ...
        'sigma0^2_mobile (um^2)';'sigma0^2_mobile_err (um^2)';'Size (nm)';'PSF_waist (nm)'; ...
        'tau_D (s)';'tau_T/tau_D';'mobile_fit_lags';'total_frames'};
    Pv=[tauT; tauT_err; R2_dec; sigma_T; sigma_T^2; D_M; D_M_err; R2_mob; ...
        offs02; offs02_err; Size_nm; PSF_waist*1000; tauD; tauT/tauD; nfit; N];
    writetable(table(Pn,Pv,'VariableNames',{'Parameter','Value'}), export_filename, 'Sheet','Parameters');
    fprintf('Exported results to: %s\n', export_filename);
end
fprintf('Done.\n');

else
%%% ==================== STANDARD single-curve output (original toolbox) =========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% power_law_iMSD
x = time_vec(:);
y = s_G(:);

modelFun = @(p, x) p(1) * x.^p(2) + p(3);   % p = [a0, alpha, c0]

p0 = [10^-2, 1, min(y)];
opts = optimoptions('lsqcurvefit','Display','off');
lb = [0, 0, 0];
ub = [Inf, 2, Inf];
[p_fit,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(modelFun, p0, x, y, lb, ub, opts);

a0 = p_fit(1);
alpha = p_fit(2);
c0 = p_fit(3);

% --- Compute R² ---
y_pred = modelFun(p_fit, x);
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res/SS_tot;

% --- Compute parameter standard errors ---
% ci = nlparci(p_fit, residual, 'jacobian', J);
% param_errors = (ci(:,2) - ci(:,1)) / (2*1.96);  % approximate standard errors (95% CI)

mse = resnorm / (length(residual) - length(p_fit));
R = chol(J'*J);
Cov = mse * inv(R) * inv(R)';
param_errors = full(sqrt(diag(Cov)));
param_errors = double(param_errors(:));

 clc
 fprintf('***********************************\niMSD output:\n\n');
 fprintf('Results for power-law fitting:\n');
 fprintf('alpha             = %.4f ± %.4f\n', alpha, param_errors(2));
 fprintf('K (mum^2/s)       = %.4f ± %.4f\n', a0, param_errors(1));
 fprintf('offset (mum^2)    = %.4f ± %.4f\n', c0, param_errors(3));
 fprintf('R²                = %.4f\n', R2);

x_fit=time_vec;
y_fit_alpha = modelFun(p_fit, x_fit);
residual_alpha=s_G-y_fit_alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alpha<1

modelFun = @(p, x) p(1) + 4*p(2)*x + (p(3)^2/3)*(1 - exp(-x/p(4)));
% p = [sigma0, D_M, L, tau_c]

% Initial guesses for parameters: [sigma0, D_M, L, tau_c]
p0 = [(min(y)), 10^-4, 0.05, 1.5];  % adjust based on your data
lb = [0, 0, 0, 0];
ub = [Inf, Inf, Inf, Inf];

% Fit using lsqcurvefit
options = optimoptions('lsqcurvefit','Display','off');
[p_fit,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(modelFun, p0, x(:), y(:), lb, ub, options);

% Calculate confidence intervals using the Jacobian
ci = nlparci(p_fit, residual, 'jacobian', J);
% Standard errors from 95% CI
% param_errors = (ci(:,2) - ci(:,1)) / (2*1.96);

mse = resnorm / (length(residual) - length(p_fit));
R = chol(J'*J);
Cov = mse * inv(R) * inv(R)';
param_errors = full(sqrt(diag(Cov)));
param_errors = double(param_errors(:));

% Compute R^2
y_fit_conf = modelFun(p_fit,x(:));
SSres = sum((y(:)-y_fit_conf).^2);
SStot = sum((y(:)-mean(y)).^2);
R2_conf = 1 - SSres/SStot;

offs02=p_fit(1);
offs02_err=param_errors(1);
D_M=p_fit(2);
D_M_error=param_errors(2);
L=p_fit(3);
tau_c=p_fit(4);
dL=param_errors(3);
dtau=param_errors(4);
D_m=p_fit(3).^2./(12*p_fit(4));
D_m_error = sqrt((2*L/(12*tau_c)*dL)^2 + (-L^2/(12*tau_c^2)*dtau)^2 );

%%% confinement strength Sconf = D_micro/D_macro
Sconf = (D_m)/max(D_M,eps);
% Sconf_error = sqrt( ...
%     (D_m_error./max(D_M,eps)).^2 + ...
%     (D_m.*D_M_error./max(D_M,eps).^2).^2 );

residual_conf=s_G-y_fit_conf;

% Display results
fprintf('\n')
fprintf('Results for confined-motion fitting\n');
fprintf('D_M (mum^2/s)     = %.4f ± %.4f\n', p_fit(2), param_errors(2));
fprintf('D_m (mum^2/s)     = %.4f ± %.4f\n', D_m, D_m_error);
fprintf('Sconf (Dm/DM)     = %.4f\n', Sconf);
fprintf('L (mum)           = %.4f ± %.4f\n', p_fit(3), param_errors(3));
fprintf('tau_c (s)         = %.4f ± %.4f\n', p_fit(4), param_errors(4));
fprintf('offset (mum^2)    = %.4f ± %.4f\n', p_fit(1), param_errors(1));
fprintf('R²                = %.4f\n', R2_conf);
fprintf('***********************************\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else   % alpha >= 1

modelFun = @(p,x) p(1) + 4*p(2)*x + (p(3)^2)*x.^2;
% p = [sigma0, D, v]

% Initial guesses
p0 = [min(y), 1e-2, 1];   % adjust if needed
lb = [0, 0, 0];
ub = [Inf, Inf, Inf];

% Fit using lsqcurvefit
options = optimoptions('lsqcurvefit','Display','off');
[p_fit,resnorm,residual,exitflag,output,lambda,J] = ...
    lsqcurvefit(modelFun, p0, x(:), y(:), lb, ub, options);

% Confidence intervals using the Jacobian
ci = nlparci(p_fit, residual, 'jacobian', J);

% Standard errors from 95% CI
param_errors = (ci(:,2) - ci(:,1)) / (2*1.96);

% Compute R^2
y_fit_conf = modelFun(p_fit, x(:));
SSres = sum((y(:) - y_fit_conf).^2);
SStot = sum((y(:) - mean(y)).^2);
R2_conf = 1 - SSres/SStot;

% Assign fitted parameters and errors
offs02      = p_fit(1);
offs02_err  = param_errors(1);

D_M         = p_fit(2);
D_M_error   = param_errors(2);

v           = p_fit(3);
v_error     = param_errors(3);

Sconf = NaN;   %%% not defined for directed motion

% Residuals
residual_conf = s_G - y_fit_conf;

fprintf('\n')
fprintf('Results for directed-motion fitting\n');
fprintf('D (mum^2/s)       = %.4f ± %.4f\n', D_M, D_M_error);
fprintf('v (mum/s)         = %.4f ± %.4f\n', v, v_error);
fprintf('offset (mum^2)    = %.4f ± %.4f\n', offs02, offs02_err);
fprintf('R²                = %.4f\n', R2_conf);
fprintf('***********************************\n');

end

% PSF_FWHM is provided in µm
% if ~exist('PSF_FWHM','var')
%     PSF_FWHM = 0.25;   % default FWHM (µm)
% end

PSF_var = (PSF_FWHM/2.355)^2;
if offs02 > PSF_var
    Size_nm = sqrt(offs02 - PSF_var)*1000;
else
    Size_nm = sqrt(max(offs02,0))*1000;
end

fprintf('Sconf                   = %.4f   (D_micro/D_macro; NaN if directed)\n', Sconf);
fprintf('Size (deconvolved)      = %.1f nm   (PSF_FWHM = %.0f nm)\n', Size_nm, PSF_FWHM*1000);
fprintf('***********************************\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

hFig = figure('Name','iMSD output','NumberTitle','off');

subplot (2,3,1)
plot1 = errorbar(time_vec,amp_G,amp_G_err,'.','LineWidth',1,'Color',[.5 .5 .5]);
xlabel('\tau (s)');
ylabel('Gaussian peak');
xlim([0 time_vec(end)+frame_time]);
ylim ([0.9*min(amp_G-amp_G_err) 1.1*max(amp_G+amp_G_err)])
set(gca,'fontsize',14)
grid on
title ('G_0 plot and corr. with 1/\sigma^2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = time_vec;
y = amp_G;
modelFun = @(p, x) p(1)*(1+x/p(2)).^(-1)+p(3);   % p(1)=a, p(2)=b
p0 = [amp_G(1), x(end), 0];
opts = optimoptions('lsqcurvefit','Display','off');
lb = [0, 0, 0];
ub = [Inf, Inf, Inf];
[p_fit,~,residual,~,~,~,J] = lsqcurvefit(modelFun, p0, x, y, lb, ub, opts);
a_stics = p_fit(1);
b_stics = p_fit(2);
c_stics = p_fit(3);
D_stics=mean(s_Gxy)^2/(4*b_stics);

y_pred_g0 = modelFun(p_fit, time_vec);
SS_res = sum((y - y_pred_g0).^2);
SS_tot = sum((y - mean(y)).^2);
R2_g0 = 1 - SS_res/SS_tot;

hold on
plot (time_vec,y_pred_g0,'r-')
text(0.05, 0.05, sprintf('D_{STICS}^{(Br.)} = %.2e \\mum^2/s',D_stics),...
    'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontSize', 12,'BackgroundColor', 'w', ...
    'EdgeColor', 'k');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax_inset=subplot (2,3,2);
scatter (1./s_G,amp_G,30,[.5 .5 .5],'filled','MarkerFaceAlpha',.5)
p = polyfit(1./s_G, amp_G, 1);
y_pred = polyval(p, 1./s_G);
c00=corr([1./s_G,amp_G]);
c_sG_amp_G=c00(1,2);
hold on
plot (1./s_G,y_pred,'r-')
hold off
text(0.05, 0.95, sprintf('Corr. = %.3f',c_sG_amp_G),...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10);
axis square
box on
grid off
set (gca,'xticklabels','')
set (gca,'yticklabels','')
xlabel ('1/\sigma^2')
ylabel ('G_0')
set (gca,'fontsize',12)

pos_inset = get(ax_inset, 'Position');  % current position [left bottom width height]
pos0_inset=pos_inset;
pos_inset(3) = pos0_inset(3) * 0.45;  % width
pos_inset(4) = pos0_inset(4) * 0.45;  % height
pos_inset(1) = pos_inset(1) - 0.17;  % move left
pos_inset(2) = (pos0_inset(2) -pos_inset(4)+ pos0_inset(4))-0.03;  % move down
set(ax_inset, 'Position', pos_inset)

subplot (2,3,4)
plot2x=errorbar(time_vec,x_G,x_G_err,'o','LineWidth',1,'Color',[.9 .6 0]);
hold on
plot2y=errorbar(time_vec,y_G,y_G_err,'s','LineWidth',1,'Color',[.6 .9 0]);
hold off

vx = polyfit(time_vec, x_G, 1);
yx_fit = polyval(vx, time_vec);
vy = polyfit(time_vec, y_G, 1);
yy_fit = polyval(vy, time_vec);

hold on
plot (time_vec,yx_fit,'r-','LineWidth',1)
plot (time_vec,yy_fit,'r-','LineWidth',1)
hold off
title ('Net flow plot')
legend ('\xi_0','\eta_0','Location','Best')

text(0.05, 0.05, sprintf('v_x=%.1f nm/s v_y=%.1f nm/s',vx(1)*10^3,vy(1)*10^3),...
    'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontSize', 12,'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

grid on
xlabel('\tau (s)');
ylabel('Gaussian location (\mum)');
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax=subplot (2,3,2);
plot3 = errorbar(time_vec,s_G,s_G_err,'o','LineWidth',1);
hold on
plot_fit_alpha=plot(x_fit,y_fit_alpha,'m--','LineWidth',1.5);


text(0.05, 0.95, sprintf(['Power-law \n\n' ...
    '\\alpha = %.3f\n\\kappa = %.3f \\mum^2/s  \n\\sigma_0^2 = %.3f \\mum^2 \nR^2=%.3f'], alpha,a0,c0,R2),...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10,'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

if alpha<1

plot_fit_conf=plot(x, y_fit_conf, 'r-','LineWidth',1.5);
hold off

annotationText = sprintf(['Confined-motion \n\n', ...
    'D_M = %.3f \\mum^2/s\n', ...
    'D_m = %.3f \\mum^2/s\n', ...
    'L = %.3f \\mum\n', ...
    '\\tau_c = %.3f s \n', ...
    '\\sigma_0^2 = %.3f \\mum^2\n', ...
    'R² = %.3f'], ...
...
    D_M, D_m,L, tau_c, offs02, R2_conf);

text(0.25, 0.95, annotationText, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

legend('Data','Power-law','Confinement','Location','Southeast')

else   % alpha >= 1

plot_fit_conf = plot(x, y_fit_conf, 'g-','LineWidth',1.5);
hold off

annotationText = sprintf(['Directed-motion \n\n', ...
    'D = %.3f \\mum^2/s\n', ...
    'v = %.3f \\mum/s\n', ...
    '\\sigma_0^2 = %.3f \\mum^2\n', ...
    'R² = %.3f'], ...
    D_M, v, offs02, R2_conf);

text(0.25, 0.95, annotationText, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

legend('Data','Power-law','Directed','Location','Southeast')

end

grid on
ylabel('Gaussian width \sigma^2 (\mum^2)');
title ('iMSD plot')
xlim([0 time_vec(end)+frame_time]);
xlabel('\tau (s)');
set(gca,'fontsize',14)
ylim ([0.9*min(s_G-s_G_err) 1.1*max(s_G+s_G_err)])

pos = get(ax, 'Position');  % current position [left bottom width height]
pos0=pos;
pos(3) = pos0(3) * 2.2;  % width
pos(4) = pos0(4) * 1.8;  % height
pos(1) = pos(1) + 0.025;  % move left
pos(2) = (pos0(2) -pos(4)+ pos0(4));  % move down
set(ax, 'Position', pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

ax2=subplot (4,2,8);
plot (time_vec,residual_alpha,'m.--','linewidth',1.5)
hold on
if alpha<1
plot (time_vec,residual_conf,'r.-','linewidth',1.5)
else   % alpha >= 1
plot(time_vec, residual_conf, 'g.-', 'LineWidth', 1.5)
end
yline (0,'b-','LineWidth',1.5)
hold off
set (gca,'fontsize',14)
grid on

pos_b = get(ax2, 'Position');  % current position [left bottom width height]
pos0_b=pos_b;
pos_b(3) = pos(3);  % width
pos_b(4) = pos0_b(4) * 0.75;  % height
 pos_b(1) = (pos0_b(1) -pos_b(3)+ pos0_b(3));  % move down
set(ax2, 'Position', pos_b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fprintf('Done');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXPORT RESULTS TO EXCEL

if export_value==1

[out_dir,~,~] = fileparts(filename);
export_filename = fullfile(out_dir, [filename '_export.xlsx']);
if isfile(export_filename), delete(export_filename); end

%% -------- Sheet 1: Data and Fits --------
T_data = table( ...
    time_vec(:), ...
    s_G(:), ...
    s_G_err(:), ...
    y_fit_alpha(:), ...
    residual_alpha(:), ...
    'VariableNames', { ...
        'time_s', ...
        'sigma2_data', ...
        'sigma2_err', ...
        'powerlaw_fit', ...
        'powerlaw_residual' ...
    });

if alpha < 1
    T_data.conf_or_dir_fit = y_fit_conf(:);
    T_data.conf_or_dir_residual = residual_conf(:);
    T_data.Properties.VariableNames(end-1:end) = ...
        {'confined_fit','confined_residual'};
else
    T_data.conf_or_dir_fit = y_fit_conf(:);
    T_data.conf_or_dir_residual = residual_conf(:);
    T_data.Properties.VariableNames(end-1:end) = ...
        {'directed_fit','directed_residual'};
end

writetable(T_data, export_filename, 'Sheet', 'Data_and_Fits');


%% -------- Sheet 2: Fit Results (Formatted) --------

ParamName  = {};
ParamValue = [];
ParamError = [];

% ----- Power-law results -----
ParamName  = [ParamName; {'Results for power-law fitting:'}];
ParamValue = [ParamValue; NaN];
ParamError = [ParamError; NaN];

ParamName  = [ParamName; {'alpha'}];
ParamValue = [ParamValue; alpha];
ParamError = [ParamError; param_errors(2)];

ParamName  = [ParamName; {'K (mum^2/s)'}];
ParamValue = [ParamValue; a0];
ParamError = [ParamError; param_errors(1)];

ParamName  = [ParamName; {'offset (mum^2)'}];
ParamValue = [ParamValue; c0];
ParamError = [ParamError; param_errors(3)];

ParamName  = [ParamName; {'R^2'}];
ParamValue = [ParamValue; R2];
ParamError = [ParamError; NaN];

% Empty line
ParamName  = [ParamName; {''}];
ParamValue = [ParamValue; NaN];
ParamError = [ParamError; NaN];

% ----- Confined or Directed results -----
if alpha < 1
    ParamName  = [ParamName; {'Results for confined-motion fitting'}];
    ParamValue = [ParamValue; NaN];
    ParamError = [ParamError; NaN];

    ParamName  = [ParamName; {'D_M (mum^2/s)'}];
    ParamValue = [ParamValue; D_M];
    ParamError = [ParamError; D_M_error];

    ParamName  = [ParamName; {'D_m (mum^2/s)'}];
    ParamValue = [ParamValue; D_m];
    ParamError = [ParamError; D_m_error];

    ParamName  = [ParamName; {'L (mum)'}];
    ParamValue = [ParamValue; L];
    ParamError = [ParamError; dL];

    ParamName  = [ParamName; {'tau_c (s)'}];
    ParamValue = [ParamValue; tau_c];
    ParamError = [ParamError; dtau];

    ParamName  = [ParamName; {'offset (mum^2)'}];
    ParamValue = [ParamValue; offs02];
    ParamError = [ParamError; offs02_err];

    ParamName  = [ParamName; {'R^2'}];
    ParamValue = [ParamValue; R2_conf];
    ParamError = [ParamError; NaN];

else
    ParamName  = [ParamName; {'Results for directed-motion fitting'}];
    ParamValue = [ParamValue; NaN];
    ParamError = [ParamError; NaN];

    ParamName  = [ParamName; {'D (mum^2/s)'}];
    ParamValue = [ParamValue; D_M];
    ParamError = [ParamError; D_M_error];

    ParamName  = [ParamName; {'v (mum/s)'}];
    ParamValue = [ParamValue; v];
    ParamError = [ParamError; v_error];

    ParamName  = [ParamName; {'offset (mum^2)'}];
    ParamValue = [ParamValue; offs02];
    ParamError = [ParamError; offs02_err];

    ParamName  = [ParamName; {'R^2'}];
    ParamValue = [ParamValue; R2_conf];
    ParamError = [ParamError; NaN];
end

% ----- confinement strength and deconvolved size (both fit types) -----
ParamName  = [ParamName; {'Sconf (D_micro/D_macro)'}];
ParamValue = [ParamValue; Sconf];
ParamError = [ParamError; NaN];

ParamName  = [ParamName; {'Size (nm)'}];
ParamValue = [ParamValue; Size_nm];
ParamError = [ParamError; NaN];

ParamName  = [ParamName; {'PSF_waist (nm)'}];
ParamValue = [ParamValue; PSF_waist*1000];
ParamError = [ParamError; NaN];

% Convert to table
T_results = table(ParamName, ParamValue, ParamError, ...
    'VariableNames', {'Parameter','Value','Error'});

% Write to Excel
writetable(T_results, export_filename, ...
    'Sheet', 'Fit_Parameters', 'WriteMode', 'overwrite');

fprintf('\nExported results to: %s\n', export_filename);

end
end

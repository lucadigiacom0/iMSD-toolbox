%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% power_law_iMSD
x = time_vec(:);
y = s_G(:);

modelFun = @(p, x) p(1) * x.^p(2) + p(3);   % p = [a0, alpha, c0]

p0 = [10^-2, 1, min(y)];   
opts = optimoptions('lsqcurvefit','Display','off');
lb = [0, 0, 0];
ub = [Inf, 2, Inf];
[p_fit,~,residual,~,~,~,J] = lsqcurvefit(modelFun, p0, x, y, lb, ub, opts);

a0 = p_fit(1);
alpha = p_fit(2);
c0 = p_fit(3);

% --- Compute R² ---
y_pred = modelFun(p_fit, x);
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res/SS_tot;

% --- Compute parameter standard errors ---
ci = nlparci(p_fit, residual, 'jacobian', J);
param_errors = (ci(:,2) - ci(:,1)) / (2*1.96);  % approximate standard errors (95% CI)

 clc
 fprintf('***********************************\niMSD output:\n\n');
 fprintf('Results for power-law fitting:\n');
 fprintf('alpha             = %.4f ± %.4f\n', alpha, param_errors(2));
 fprintf('K (mum^2/s)      = %.4f ± %.4f\n', a0/4, param_errors(1)/4);
 fprintf('offset (mum^2)    = %.4f ± %.4f\n', c0, param_errors(3));
 fprintf('R²                = %.4f\n', R2);

x_fit=time_vec;
y_fit_alpha = modelFun(p_fit, x_fit);
residual_alpha=s_G-y_fit_alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alpha<1

modelFun = @(p, x) p(1) + 4*p(2)*x + (p(3)^3/3)*(1 - exp(-x/p(4)));
% p = [sigma0, D_M, L, tau_c]

% Initial guesses for parameters: [sigma0, D_M, L, tau_c]
p0 = [(min(y)), 10^-2, 0.5, 1];  % adjust based on your data
lb = [0, 0, 0, 0];
ub = [Inf, Inf, Inf, Inf];

% Fit using lsqcurvefit
options = optimoptions('lsqcurvefit','Display','off');
[p_fit,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(modelFun, p0, x(:), y(:), lb, ub, options);

% Calculate confidence intervals using the Jacobian
ci = nlparci(p_fit, residual, 'jacobian', J);
% Standard errors from 95% CI
param_errors = (ci(:,2) - ci(:,1)) / (2*1.96);

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
D_m=p_fit(3).^3./(12*p_fit(4));
D_m_error = sqrt( (3*L^2/(12*tau_c)*dL)^2 + (-L^3/(12*tau_c^2)*dtau)^2 );

residual_conf=s_G-y_fit_conf;

% Display results
% fprintf('**********************************\niMSD output:\n\n');
fprintf('\n')
fprintf('Results for confined-motion fitting\n');
fprintf('D_M (mum^2/s)     = %.4f ± %.4f\n', p_fit(2), param_errors(2));
fprintf('D_m (mum^2/s)     = %.4f ± %.4f\n', D_m, D_m_error);
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
% axis square
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
% Increase width and height by 20%, for example
pos_inset(3) = pos0_inset(3) * 0.45;  % width
pos_inset(4) = pos0_inset(4) * 0.45;  % height
% Optionally shift it slightly if overlapping occurs
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
% axis square 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax=subplot (2,3,2);
plot3 = errorbar(time_vec,s_G,s_G_err,'o','LineWidth',1);
hold on
plot_fit_alpha=plot(x_fit,y_fit_alpha,'m--','LineWidth',1.5);


text(0.05, 0.95, sprintf(['Power-law \n\n' ...
    '\\alpha = %.3f\n\\kappa = %.3f \\mum^2/s  \n\\sigma_0^2 = %.3f \\mum^2 \nR^2=%.3f'], alpha,a0/4,c0,R2),...
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
% Increase width and height by 20%, for example
pos(3) = pos0(3) * 2.2;  % width
pos(4) = pos0(4) * 1.8;  % height
% Optionally shift it slightly if overlapping occurs
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
% Increase width and height by 20%, for example
pos_b(3) = pos(3);  % width
pos_b(4) = pos0_b(4) * 0.75;  % height
% Optionally shift it slightly if overlapping occurs
% pos_b(1) = pos_b(1) - 0.025;  % move left
 pos_b(1) = (pos0_b(1) -pos_b(3)+ pos0_b(3));  % move down
set(ax2, 'Position', pos_b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fprintf('Done');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXPORT RESULTS TO EXCEL

if export_value==1

export_filename = [filename '_export.xlsx'];

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
ParamValue = [ParamValue; a0/4];
ParamError = [ParamError; param_errors(1)/4];

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

% Convert to table
T_results = table(ParamName, ParamValue, ParamError, ...
    'VariableNames', {'Parameter','Value','Error'});

% Write to Excel
writetable(T_results, export_filename, ...
    'Sheet', 'Fit_Parameters', 'WriteMode', 'overwrite');

fprintf('\nExported results to: %s\n', export_filename);

end

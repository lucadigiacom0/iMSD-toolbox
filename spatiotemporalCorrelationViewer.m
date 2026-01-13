function spatiotemporalCorrelationViewer(time_vec,frame_time, px_size, ...
    a1, b1, nt, STCorr_a, Gaussb, amp_G, var_G, o_G, R2_G,...
    amp_G_err, var_G_err, o_G_err)

hFig = figure('Name','SpatioTemporal Correlation Function','NumberTitle','off');
set(hFig, 'Tag', char(java.util.UUID.randomUUID));

  vector1_for_limits=[min(min(min(STCorr_a))) 1.1*max(max(max(STCorr_a)))];
  vector2_for_limits=[min(min(min(Gaussb))) 1.1*max(max(max(Gaussb)))];
  vector_limits=mean([vector1_for_limits; vector2_for_limits]);

% --- Subplot Data ---
ax1 = subplot(2,4,1);
s1 = surf(a1,b1,STCorr_a(:,:,1));
axis square
set(gca,'fontsize',14)
box on
title('Data. Frame-lag 1')
zlabel ('G(\xi,\eta,\tau)')
xlabel ('\xi (\mum)')
ylabel ('\eta (\mum)')
zlim (vector_limits)
set(ax1, 'Projection','perspective');

% --- Subplot Fit ---
ax2 = subplot(2,4,2);
s2 = surf(a1,b1,Gaussb(:,:,1));
axis square
set(gca,'fontsize',14)
box on
title('Fit. Frame-lag 1')
zlabel ('G_{fit}(\xi,\eta,\tau)')
xlabel ('\xi (\mum)')
ylabel ('\eta (\mum)')
zlim (vector_limits)
set(ax2, 'Projection','perspective');

% --- Subplot Residual ---
ax3 = subplot(2,4,3);
pl3 = pcolor(a1-px_size/2,b1-px_size/2,-Gaussb(:,:,1)+STCorr_a(:,:,1));
set(pl3,'EdgeColor','none')
c = colorbar(ax3,'southoutside');
c.Label.String='G_{fit}(\xi,\eta,\tau)-G(\xi,\eta,\tau)';
c.Label.FontSize=14;
axis square
set(gca,'fontsize',14)
title('Residuals. Frame-lag 1')
colormap(jet)
clim ([min(min(min(Gaussb-STCorr_a))) max(max(max(Gaussb-STCorr_a)))])

 x_min = min(a1-px_size/2); 
    x_max = max(a1-px_size/2);
    y_min = min(b1-px_size/2);
    y_max = max(b1-px_size/2);
rectangle('Position',[x_min, y_min, x_max - x_min, y_max - y_min], ...
          'EdgeColor','k', 'LineWidth',.5, 'LineStyle','-');


% ------
ax4 = subplot(2,4,4);
hold on
idx_x = round(size(STCorr_a,2)/2);
idx_y = round(size(STCorr_a,1)/2);
hZx = plot(a1, STCorr_a(idx_y,:,1), 'o-','Color',[.9 .6 0], 'DisplayName','G(\xi,0,\tau)',...
    'MarkerSize',8,'LineWidth',1.5);
% hGx = plot(a1, Gaussb(idx_y,:,1), 'k-', 'DisplayName','Fit (x-sec)');
hZy = plot(b1, STCorr_a(:,idx_x,1), 's-','Color',[.6 .9 0], 'DisplayName','G(0,\eta,\tau)',...
    'MarkerSize',8,'LineWidth',1.5);
hGy = plot(b1, Gaussb(:,idx_x,1), 'k-', 'DisplayName','Fit','HandleVisibility','off',...
    'linewidth',2);
hGx = plot(a1, Gaussb(idx_y,:,1), 'k-', 'DisplayName','Fit','LineWidth',2);
ylim (vector_limits)
% xlabel('Coordinate (\mum)')
% ylabel('G(\xi,\eta,\tau)')
title('Sections. Frame-lag 1')
legend('Location','South')
grid on
box on
set(gca,'fontsize',14)

% ------
ax5 = subplot(2,4,5);
% p5 = plot(amp_G,'.-'); hold on; yline(mean(amp_G));
p5 = errorbar(time_vec,amp_G,amp_G_err,'.-','LineWidth',1.5); hold on;
xlabel('\tau (s)'); ylabel('Amplitude G_0'); xlim([0 time_vec(end)+frame_time]); set(gca,'fontsize',14)
ylim ([0.9*min(amp_G-amp_G_err) 1.1*max(amp_G+amp_G_err)])

% ------
ax6 = subplot(2,4,7);
p6 = errorbar(time_vec,var_G,var_G_err,'.-','LineWidth',1.5); hold on;
xlabel('\tau (s)'); ylabel('Variance \sigma^2 (\mum^2)'); xlim([0 time_vec(end)+frame_time]); set(gca,'fontsize',14)
ylim ([0.9*min(var_G-var_G_err) 1.1*max(var_G+var_G_err)])

% ------
ax7 = subplot(2,4,6);
p7 = errorbar(time_vec,o_G,o_G_err,'.-','LineWidth',1.5); hold on;
xlabel('\tau (s)'); ylabel('Offset G_\infty'); xlim([0 time_vec(end)+frame_time]); set(gca,'fontsize',14)
ylim (sort([0.9*min(o_G-o_G_err) 1.1*max(o_G+o_G_err)]))

% ------
ax8 = subplot(2,4,8);
p8 = plot(time_vec,R2_G,'.-','LineWidth',1.5); hold on; yline(mean(R2_G));
xlabel('\tau (s)'); ylabel('R^2'); xlim([0 time_vec(end)+frame_time]); set(gca,'fontsize',14)
ylim ([min(R2_G) max(R2_G)])

% ------
vline(1) = line(ax5,frame_time*[1 1],ylim(ax5),'Color','r','LineWidth',2,'HandleVisibility','off');
vline(2) = line(ax6,frame_time*[1 1],ylim(ax6),'Color','r','LineWidth',2,'HandleVisibility','off');
vline(3) = line(ax7,frame_time*[1 1],ylim(ax7),'Color','r','LineWidth',2,'HandleVisibility','off');
vline(4) = line(ax8,frame_time*[1 1],ylim(ax8),'Color','r','LineWidth',2,'HandleVisibility','off');

% ------
currentTau = 1;

% --- Slider ---
hSlider = uicontrol('Parent',hFig,'Style','slider','Units','normalized',...
    'Position',[0.13 0.48 0.57 0.025],...
    'Min',1,'Max',double(nt),'Value',currentTau,...
    'SliderStep',[1/(double(nt)-1), 10/(double(nt)-1)],...
    'Callback',@sliderCallback);

% --- Play/Pause button ---
hButton = uicontrol('Parent',hFig,'Style','togglebutton','Units','normalized',...
    'Position',[0.76 0.48 0.14 0.05],...
    'String','Play');

% --- Timer ---
t = timer('ExecutionMode','fixedRate','Period',0.1,'TimerFcn',@timerCallback);

% --- Slider callback ---
function sliderCallback(~,~)
    currentTau = round(get(hSlider,'Value'));           % aggiorna currentTau
    currentTau = max(1,min(nt,currentTau));            % clamp
    updateTau(currentTau);                             % aggiorna grafici
end

% --- Timer callback ---
function timerCallback(~,~)
    if get(hButton,'Value') == 1                        % solo se Play attivo
        currentTau = currentTau + 1;
        if currentTau > nt
            currentTau = 1;
        end
        set(hSlider,'Value',currentTau);               % aggiorna slider
        updateTau(currentTau);
    end
end

% ------
function updateTau(tau)
    % Aggiorna superfici
    set(s1,'ZData',STCorr_a(:,:,tau));
    set(s2,'ZData',Gaussb(:,:,tau));
    set(pl3,'CData',-Gaussb(:,:,tau)+STCorr_a(:,:,tau));

    % 
    set(hZx,'YData',STCorr_a(idx_y,:,tau));
    set(hGx,'YData',Gaussb(idx_y,:,tau));
    set(hZy,'YData',STCorr_a(:,idx_x,tau));
    set(hGy,'YData',Gaussb(:,idx_x,tau));

    % 
    title(ax1,sprintf('Data. Frame-lag %d',tau));
    title(ax2,sprintf('Fit. Frame-lag %d',tau));
    title(ax3,sprintf('Residuals. Frame-lag %d',tau));
    title(ax4,sprintf('Sections. Frame-lag %d',tau));

    % 
    for k = 1:4
        set(vline(k),'XData',[tau*frame_time tau*frame_time]);
    end
end

% --- Play/Pause callback ---
hButton.Callback = @(src,~) togglePlay(src);
function togglePlay(src)
    if get(src,'Value') == 1
        src.String = 'Pause';
        start(t);
    else
        src.String = 'Play';
        stop(t);
    end
end

% --- Cleanup timer ---
hFig.CloseRequestFcn = @closeFigure;
function closeFigure(~,~)
    stop(t); delete(t); delete(hFig);
end

end


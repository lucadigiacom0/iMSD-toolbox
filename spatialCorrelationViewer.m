function spatialCorrelationViewer(a1,b1,Z2,Gaussxy,px_size,a_Gxy,s_Gxy,o_Gxy,R2_Gxy)
    
    nFrames = size(Z2,3);
    
    vector1_for_limits=[min(min(min(Z2))) max(max(max(Z2)))];
    vector2_for_limits=[min(min(min(Z2))) max(max(max(Z2)))];
    vector_limits=mean([vector1_for_limits; vector2_for_limits]);

    hFig = figure('Name','Spatial Correlation Function','NumberTitle','off');
    set(hFig, 'Tag', char(java.util.UUID.randomUUID))

    % --- Subplot Data ---
    ax1 = subplot(2,4,1);
    s1 = surf(a1,b1,Z2(:,:,1));
    axis square
    set(gca,'fontsize',14)
    box on
    title('Data. Frame 1')
    zlabel ('g(\xi,\eta)')
    xlabel ('\xi (\mum)')
    ylabel ('\eta (\mum)')
    zlim (vector_limits)
    set(ax1, 'Projection','perspective');

    % --- Subplot Fit ---
    ax2 = subplot(2,4,2);
    s2 = surf(a1,b1,Gaussxy(:,:,1));
    axis square
    set(gca,'fontsize',14)
    box on
    title('Fit. Frame 1')
    zlabel ('g_{fit}(\xi,\eta)')
    xlabel ('\xi (\mum)')
    ylabel ('\eta (\mum)')
    zlim (vector_limits)
    set(ax2, 'Projection','perspective');

    % --- Subplot Residual ---
    ax3 = subplot(2,4,3);
    pl3 = pcolor(a1-px_size/2,b1-px_size/2,-Gaussxy(:,:,1)+Z2(:,:,1));
    set(pl3,'edgecolor','none')
    c = colorbar(ax3,'southoutside');
    c.Label.String='g_{fit}(\xi,\eta)-g(\xi,\eta)';
    c.Label.FontSize=14;
    box on
    axis square
    set(gca,'fontsize',14)
    title('Residuals. Frame 1')
    colormap(jet)
    set (gca,'xtick','')
    set (gca,'ytick','')
    x_min = min(a1-px_size/2); 
    x_max = max(a1-px_size/2);
    y_min = min(b1-px_size/2);
    y_max = max(b1-px_size/2);
rectangle('Position',[x_min, y_min, x_max - x_min, y_max - y_min], ...
          'EdgeColor','k', 'LineWidth',.5, 'LineStyle','-');

    %  xlabel ('\xi (px)')
    % ylabel ('\eta (px)')
    clim ([min(min(min(Gaussxy-Z2))) max(max(max(Gaussxy-Z2)))])
    % box on

    axSec = subplot(2,4,4);
    hold on
    idx_x = round(size(Z2,2)/2); 
    idx_y = round(size(Z2,1)/2); %
    hZx = plot(a1, Z2(idx_y,:,1), 'o-','Color',[.9 .6 0], 'DisplayName','g(\xi,0)','MarkerSize',8,'LineWidth',1.5);
    hZy = plot(b1, Z2(:,idx_x,1), 's-','Color',[.6 .9 0], 'DisplayName','g(0,\eta)','MarkerSize',8,'LineWidth',1.5);
    hGx = plot(a1, Gaussxy(idx_y,:,1), 'k-', 'DisplayName','g_{fit}','LineWidth',2);
    % hGy = plot(b1, Gaussxy(:,idx_x,1), 'm--', 'DisplayName','Fit (y-sec)');
    box on
    % ylabel('g')
    title('Sections (Frame 1)')
    legend('Location','South')
    grid on
    % axis square
    set(gca,'fontsize',14)
    


    % --- Subplot a_Gxy ---
    ax4 = subplot(2,4,5);
    p4 = plot(a_Gxy,'.-'); hold on
    yline(mean(a_Gxy));
    xlabel('Frame')
    ylabel('Amplitude')
    xlim([0 nFrames+1])
    set(gca,'fontsize',14)

    % --- Subplot s_Gxy^2 ---
    ax5 = subplot(2,4,7);
    p5 = plot(s_Gxy.^2,'.-'); hold on
    yline(mean(s_Gxy.^2));
    xlabel('Frame')
    ylabel('Variance (\mum^2)')
    xlim([0 nFrames+1])
    set(gca,'fontsize',14)

    % --- Subplot o_Gxy ---
    ax6 = subplot(2,4,6);
    p6 = plot(o_Gxy,'.-'); hold on
    yline(mean(o_Gxy));
    xlabel('Frame')
    ylabel('Offset')
    xlim([0 nFrames+1])
    set(gca,'fontsize',14)

    % --- Subplot R2_Gxy ---
    ax7 = subplot(2,4,8);
    p7 = plot(R2_Gxy,'.-'); hold on
    yline(mean(R2_Gxy));
    xlabel('# frame')
    ylabel('R^2')
    xlim([0 nFrames+1])
    set(gca,'fontsize',14)

    % --- Linee verticali rosse nei 4 subplot ---
    vline(1) = line(ax4,[1 1],ylim(ax4),'Color','r','LineWidth',2,'HandleVisibility','off');
    vline(2) = line(ax5,[1 1],ylim(ax5),'Color','r','LineWidth',2,'HandleVisibility','off');
    vline(3) = line(ax6,[1 1],ylim(ax6),'Color','r','LineWidth',2,'HandleVisibility','off');
    vline(4) = line(ax7,[1 1],ylim(ax7),'Color','r','LineWidth',2,'HandleVisibility','off');

   
% --- Slider ---
hSlider = uicontrol('Parent',hFig,'Style','slider','Units','normalized',...
    'Position',[0.13 0.48 0.57 0.025],...
    'Min',1,'Max',nFrames,'Value',1,...
    'SliderStep',[1/(nFrames-1) , 10/(nFrames-1)],...
    'Callback',@sliderCallback);

% --- Play/Pause button ---
hButton = uicontrol('Parent',hFig,'Style','togglebutton','Units','normalized',...
    'Position',[0.76 0.48 0.14 0.05],...
    'String','Play');

% --- Timer ---
t = timer('ExecutionMode','fixedRate','Period',0.1,'TimerFcn',@timerCallback);

% --- Callback slider ---
function sliderCallback(~,~)
    frame = round(get(hSlider,'Value'));
    frame = max(1, min(nFrames, frame)); % clamp per sicurezza
    updateFrame(frame);
end

% --- Callback timer ---
function timerCallback(~,~)
    frame = round(get(hSlider,'Value')) + 1;
    if frame > nFrames
        frame = 1;
    end
    set(hSlider,'Value',frame);
    updateFrame(frame);
end

% --- Funzione per aggiornare il frame ---
function updateFrame(frame)
    % Aggiorna superfici
    set(s1,'ZData',Z2(:,:,frame));
    set(s2,'ZData',Gaussxy(:,:,frame));
    set(pl3,'CData',-Gaussxy(:,:,frame)+Z2(:,:,frame));

    % Aggiorna linee verticali
    for k = 1:4
        set(vline(k),'XData',[frame frame],...
                     'YData',ylim(eval(['ax' num2str(3+k)])));
    end

    % Aggiorna titoli
    title(ax1,sprintf('Data. Frame %d',frame));
    title(ax2,sprintf('Fit. Frame %d',frame));
    title(ax3,sprintf('Residuals. Frame %d',frame));

    % Aggiorna sezioni centrali
    idx_x = round(size(Z2,2)/2);
    idx_y = round(size(Z2,1)/2);
    set(hZx,'YData',Z2(idx_y,:,frame));
    set(hGx,'YData',Gaussxy(idx_y,:,frame));
    set(hZy,'YData',Z2(:,idx_x,frame));
    % set(hGy,'YData',Gaussxy(:,idx_x,frame));

    % Aggiorna titolo subplot sezioni
    title(axSec,sprintf('Sections. Frame %d',frame));
end

% --- Play/Pause button callback ---
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
    stop(t);
    delete(t);
    delete(hFig);
end

end

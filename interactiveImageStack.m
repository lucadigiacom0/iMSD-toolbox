function interactiveImageStack(imgStack)
    % imgStack is assumed to be (rows x cols x frames)
    
    % figure (3)
    % Create figure
    hFig = figure('Name','Image Stack Viewer','NumberTitle','off');
    set(hFig, 'Tag', char(java.util.UUID.randomUUID))
    nFrames = size(imgStack,3);

    % Subplot for image
    ax1 = subplot(1,3,1);
    % hImg = imshow(imgStack(:,:,1),[], 'Parent',ax1);
    pl1=pcolor(imgStack(:,:,1));
    set (pl1,'edgecolor','none');
    colormap (jet)
    axis square 
    title(ax1,'Frame 1');
    set (gca,'fontsize',14)
    xlabel ('x (px)')
    ylabel ('y (px)')
    clim ([min(min(min(imgStack))) max(max(max(imgStack)))])

    % Subplot for histogram
    ax2 = subplot(1,3,2);
    hHist = histogram(imgStack(:,:,1),32);
    title(ax2,'Histogram. Frame 1');
    axis square 
    % set (gca,'xscale','log')
    set (gca,'fontsize',14)
    xlabel ('Intensity')
    ylabel ('Counts')

    ax3=subplot(1,3,3);
    pl3(1)=plot (1:size(imgStack,3),reshape(mean(mean(imgStack)),[],1),'.-','LineWidth',1);
    hold on
    pl3(2)=plot (1:size(imgStack,3),reshape(min(min(imgStack)),[],1),'.-','Color',[.6 .6 .6],'LineWidth',1);
    pl3(3)=plot (1:size(imgStack,3),reshape(max(max(imgStack)),[],1),'.-','Color',[.9 .6 0],'LineWidth',1);
    hold off
    legend(pl3,'Mean','Min','Max','Location','NorthEast')
    % set (gca,'yscale','log')
    set (gca,'fontsize',14)
    axis square
    ylim ([min(min(min(imgStack))) max(max(max(imgStack)))])
    ylabel ('Intensity')
    xlabel ('# frame')
    vline = line([1 1], ylim(ax3), 'Color','r','LineWidth',2,'HandleVisibility','off');
    title ('Intensity plot')

    % Slider to move through frames
    hSlider = uicontrol('Style','slider','Units','normalized',...
        'Position',[0.12 0.05 0.5 0.05],...
        'Min',1,'Max',nFrames,'Value',1,'SliderStep',[1/(nFrames-1) , 10/(nFrames-1)],...
        'Callback',@sliderCallback);

    % Play/Pause button
    hButton = uicontrol('Style','togglebutton','Units','normalized',...
        'Position',[0.7 0.05 0.2 0.05],...
        'String','Play');

    % Timer for animation
    t = timer('ExecutionMode','fixedRate','Period',0.1,'TimerFcn',@timerCallback);

    % Nested callbacks
    function sliderCallback(~,~)
        frame = round(get(hSlider,'Value'));
        updateFrame(frame);
    end

    function timerCallback(~,~)
        frame = round(get(hSlider,'Value')) + 1;
        if frame > nFrames
            frame = 1;
        end
        set(hSlider,'Value',frame);
        updateFrame(frame);
    end

    function updateFrame(frame)
        set(pl1,'CData',imgStack(:,:,frame));
        hHist.Data = imgStack(:,:,frame);
        title(ax1,sprintf('Frame %d',frame));
        title(ax2,sprintf('Histogram. Frame %d',frame));
        set(vline,'XData',[frame frame], 'YData',ylim(ax3));

    end

    % Button behavior
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

    % Clean up timer on figure close
        % Clean up timer on figure close
    hFig.CloseRequestFcn = @closeFigure;

    function closeFigure(~,~)
        stop(t);
        delete(t);
        delete(hFig);
    end

end
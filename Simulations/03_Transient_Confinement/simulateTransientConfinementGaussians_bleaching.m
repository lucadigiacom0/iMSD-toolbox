function [imgStack, traj, msd_all, msd_mean, msd_time] = ...
simulateTransientConfinementGaussians_bleaching( ...
    N, imgSize, nFrames, D,...
    L, P_cross,...
    meanSize, stdSize,...
    intensity, noiseLevel,...
    k_bleach, saveName)

% ------------------------------------------------------------
% simulateTransientConfinementGaussians_bleaching
%
% Brownian diffusion inside square compartments.
%
% Compartments are separated by infinitely thin barriers.
% Each barrier crossing is accepted with probability P_cross.
%
% P_cross = 1  --> free Brownian diffusion
% P_cross = 0  --> permanent confinement
%
% D       : px^2/frame
% L       : compartment size (px)
% P_cross : crossing probability
%
% ------------------------------------------------------------

rng('shuffle');

H = imgSize(1);
W = imgSize(2);

imgStack = zeros(H,W,nFrames);

[X,Y] = meshgrid(1:W,1:H);

%% ------------------------------------------------------------
% Initial particle positions

x = W*rand(N,1);
y = H*rand(N,1);

sigma = max(0.5,meanSize + stdSize*randn(N,1));

stepStd = sqrt(2*D);

%% ------------------------------------------------------------
% Store trajectories

traj.x = zeros(N,nFrames);
traj.y = zeros(N,nFrames);

%% ------------------------------------------------------------
% Simulation

for t = 1:nFrames

    frame = zeros(H,W);

    bleachFactor = exp(-k_bleach*(t-1));

    %% Draw particles

    for i = 1:N

        currentIntensity = intensity*bleachFactor;

        gauss = currentIntensity .* ...
            exp(-((X-x(i)).^2 + (Y-y(i)).^2)/(2*sigma(i)^2));

        frame = frame + gauss;

    end

    %% Add noise

    frame = frame + noiseLevel*randn(H,W);
    frame(frame<0)=0;

    imgStack(:,:,t)=frame;

    %% Save trajectories

    traj.x(:,t)=x;
    traj.y(:,t)=y;

       %% ------------------------------------------------------------
    %% ------------------------------------------------------------
% Brownian proposed step

dx = stepStd*randn(N,1);
dy = stepStd*randn(N,1);

%% ------------------------------------------------------------
% Apply compartment barriers

for i = 1:N

    x0 = x(i);
    y0 = y(i);

    dxi = dx(i);
    dyi = dy(i);

    %% ---------- Vertical barriers ----------

    cell_old = floor((x0-1)/L);
    cell_new = floor((x0+dxi-1)/L);

    if cell_old ~= cell_new

        % Crossing attempted

        if rand > P_cross

            % Barrier blocks the motion normal to itself
            dxi = 0;

        end

    end

    %% ---------- Horizontal barriers ----------

    cell_old = floor((y0-1)/L);
    cell_new = floor((y0+dyi-1)/L);

    if cell_old ~= cell_new

        if rand > P_cross

            dyi = 0;

        end

    end

    %% ---------- Update position

    x(i) = x0 + dxi;
    y(i) = y0 + dyi;

end

%% ------------------------------------------------------------
% Reflective borders of the image

x(x<1) = 2-x(x<1);
x(x>W) = 2*W-x(x>W);

y(y<1) = 2-y(y<1);
y(y>H) = 2*H-y(y>H);

end

%% ------------------------------------------------------------

imgStack = imgStack ./ max(imgStack(:));
%% ==========================================================
% Compute MSD

maxLag = nFrames-1;

msd_time = 1:maxLag;

msd_all = zeros(N,maxLag);

for lag = 1:maxLag

    dx = traj.x(:,1+lag:end) - traj.x(:,1:end-lag);
    dy = traj.y(:,1+lag:end) - traj.y(:,1:end-lag);

    msd_all(:,lag) = mean(dx.^2 + dy.^2,2);

end

msd_mean = mean(msd_all,1);

%% ==========================================================
% Normalize image stack

imgStack = imgStack ./ max(imgStack(:));

%% ==========================================================
% Save TIFF

if nargin>11 && ~isempty(saveName)

    fprintf('Saving image stack to %s...\n',saveName);

    for t=1:nFrames

        frame8 = im2uint8(imgStack(:,:,t));

        if t==1

            imwrite(frame8,...
                saveName,...
                'tif',...
                'Compression','none');

        else

            imwrite(frame8,...
                saveName,...
                'tif',...
                'WriteMode','append',...
                'Compression','none');

        end

    end

    fprintf('Saved %d frames.\n',nFrames);

end

%% ==========================================================
% Display

figure('Color','w');

subplot(2,3,1)

imagesc(imgStack(:,:,1))
axis image off
colormap hot
title('Frame 1')

subplot(2,3,2)

imagesc(imgStack(:,:,round(nFrames/2)))
axis image off
title('Middle frame')

subplot(2,3,3)

imagesc(imgStack(:,:,end))
axis image off
title('Last frame')

subplot(2,3,4)

plot(traj.x',traj.y','-')
axis equal
xlabel('x (px)')
ylabel('y (px)')
title('Particle trajectories')

%% Show compartment grid

hold on

for xx=L:L:W
    xline(xx,'k:')
end

for yy=L:L:H
    yline(yy,'k:')
end

subplot(2,3,[5 6])

plot(msd_time,msd_all',...
    'Color',[0.8 0.8 0.8])

hold on

plot(msd_time,msd_mean,...
    'k',...
    'LineWidth',2)

xlabel('Lag (\Delta t) [frames]')

ylabel('MSD [px^2]')

grid on

legend('Single particles',...
       'Mean MSD',...
       'Location','northwest')

title('Mean Square Displacement')

sgtitle(sprintf(['Transient confinement   ',...
    'L = %.1f px   ',...
    'P = %.3f   ',...
    'D = %.3f px^2/frame'],...
    L,P_cross,D))

%% ==========================================================
% Export variables

assignin('base','traj',traj);
assignin('base','msd_all',msd_all);
assignin('base','msd_mean',msd_mean);
assignin('base','msd_time',msd_time);

end
function [imgStack, traj, msd_all, msd_mean, msd_time] = ...
    simulateGuidedDiffusionGaussians_bleaching( ...
    N, imgSize, nFrames, ...
    Din, Dout, v, ...
    pathWidth, pathLength, Npaths, ...
    meanSize, stdSize, ...
    intensity, noiseLevel, ...
    k_bleach, saveName)

% ==============================================================
% SIMULATE GUIDED DIFFUSION
%
% Particles diffuse freely outside the guided regions (Dout)
% and diffuse with coefficient Din plus a constant drift v
% along randomly-oriented linear paths.
%
% ==============================================================

rng('shuffle');

%% -------------------------------------------------------------
% Image size

H = imgSize(1);
W = imgSize(2);

imgStack = zeros(H,W,nFrames);

[X,Y] = meshgrid(1:W,1:H);

%% -------------------------------------------------------------
% Particle initialization

x = W*rand(N,1);
y = H*rand(N,1);

sigma = max(0.5,meanSize + stdSize*randn(N,1));

traj.x = zeros(N,nFrames);
traj.y = zeros(N,nFrames);

%% -------------------------------------------------------------
% Generate random network of straight paths

paths = struct();

for k = 1:Npaths

    % random centre

    xc = rand*W;
    yc = rand*H;

    % random orientation

    theta = 2*pi*rand;

    ux = cos(theta);
    uy = sin(theta);

    % segment endpoints

    x1 = xc - pathLength/2*ux;
    y1 = yc - pathLength/2*uy;

    x2 = xc + pathLength/2*ux;
    y2 = yc + pathLength/2*uy;

    paths(k).x1 = x1;
    paths(k).y1 = y1;

    paths(k).x2 = x2;
    paths(k).y2 = y2;

    paths(k).ux = ux;
    paths(k).uy = uy;

end

%% -------------------------------------------------------------
% Brownian step amplitudes

stepIn  = sqrt(2*Din);

stepOut = sqrt(2*Dout);

%% -------------------------------------------------------------
% Start simulation

for t = 1:nFrames

    frame = zeros(H,W);

    bleachFactor = exp(-k_bleach*(t-1));

    traj.x(:,t)=x;
    traj.y(:,t)=y;

        %% ==========================================================
    % Move every particle

    for i = 1:N

        minDist = inf;

        guideDir = [0 0];

        %% ------------------------------------------------------
        % Search nearest path

        for k = 1:Npaths

            A = [paths(k).x1 paths(k).y1];
            B = [paths(k).x2 paths(k).y2];

            P = [x(i) y(i)];

            AB = B-A;

            AP = P-A;

            tproj = dot(AP,AB)/dot(AB,AB);

            tproj = max(0,min(1,tproj));

            Q = A + tproj*AB;

            dist = norm(P-Q);

            if dist < minDist

                minDist = dist;

                guideDir = [paths(k).ux paths(k).uy];

            end

        end

        %% ------------------------------------------------------
        % Choose diffusion mode

        if minDist <= pathWidth/2

            % Guided region

            dx = stepIn*randn + v*guideDir(1);

            dy = stepIn*randn + v*guideDir(2);

        else

            % Free diffusion

            dx = stepOut*randn;

            dy = stepOut*randn;

        end

        %% ------------------------------------------------------
        % Update position

        x(i)=x(i)+dx;

        y(i)=y(i)+dy;

        %% ------------------------------------------------------
        % Reflective image borders

        if x(i)<1
            x(i)=2-x(i);
        end

        if x(i)>W
            x(i)=2*W-x(i);
        end

        if y(i)<1
            y(i)=2-y(i);
        end

        if y(i)>H
            y(i)=2*H-y(i);
        end

        %% ------------------------------------------------------
        % Draw Gaussian spot

        gauss = bleachFactor*intensity*exp( ...
            -((X-x(i)).^2 + (Y-y(i)).^2) ...
            /(2*sigma(i)^2));

        frame = frame + gauss;

    end

    %% ----------------------------------------------------------
    % Add detector noise

    frame = frame + noiseLevel*randn(H,W);

    frame(frame<0)=0;

    imgStack(:,:,t)=frame;

    
    
end

%% ==============================================================
% Normalize image stack

imgStack = imgStack ./ max(imgStack(:));

%% ==============================================================
% Compute MSD

maxLag = nFrames - 1;

msd_time = 1:maxLag;

msd_all = zeros(N,maxLag);

for lag = 1:maxLag

    dx = traj.x(:,1+lag:end) - traj.x(:,1:end-lag);

    dy = traj.y(:,1+lag:end) - traj.y(:,1:end-lag);

    msd_all(:,lag) = mean(dx.^2 + dy.^2,2);

end

msd_mean = mean(msd_all,1);

%% ==============================================================
% Save TIFF stack

if nargin > 14 && ~isempty(saveName)

    fprintf('Saving image stack to %s ...\n',saveName);

    for t = 1:nFrames

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

    fprintf('Done.\n');

end

%% ==============================================================
% Display results

figure('Color','w');

subplot(2,3,1)

imagesc(imgStack(:,:,1));

axis image off

title('Frame 1')

colormap hot

subplot(2,3,2)

imagesc(imgStack(:,:,round(nFrames/2)));

axis image off

title('Middle frame')

subplot(2,3,3)

imagesc(imgStack(:,:,end));

axis image off

title('Last frame')

%% --------------------------------------------------------------
% Plot paths

subplot(2,3,4)

hold on

for k=1:Npaths

    plot([paths(k).x1 paths(k).x2],...
         [paths(k).y1 paths(k).y2],...
         'b-','LineWidth',2);

end

plot(traj.x',traj.y','r')

axis equal

xlim([1 W])

ylim([1 H])

title('Trajectories and guided paths')

xlabel('x (px)')

ylabel('y (px)')

%% --------------------------------------------------------------
% MSD

subplot(2,3,[5 6])

plot(msd_time,msd_all',...
    'Color',[0.8 0.8 0.8]);

hold on

plot(msd_time,...
    msd_mean,...
    'k',...
    'LineWidth',2)

xlabel('Lag (frames)')

ylabel('MSD (px^2)')

title('Mean Square Displacement')

grid on

%% ==============================================================
% Export variables

assignin('base','traj',traj);

assignin('base','msd_all',msd_all);

assignin('base','msd_mean',msd_mean);

assignin('base','msd_time',msd_time);

assignin('base','paths',paths);

end
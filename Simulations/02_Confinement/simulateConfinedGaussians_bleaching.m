function [imgStack, traj, msd_all, msd_mean, msd_time] = ...
    simulateConfinedGaussians_bleaching(N, imgSize, nFrames, L, tau_c, meanSize, stdSize,...
    intensity, noiseLevel, k_bleach, saveName)

% simulateConfinedGaussians_bleaching
%
% Ornstein-Uhlenbeck confined diffusion
%
% MSD(tau) = (L^2/3)*(1-exp(-tau/tau_c))
%
% INPUTS
% L        : confinement length (pixels)
% tau_c    : confinement time (frames)

%% --------------------------------------------------------
rng('shuffle');

H = imgSize(1);
W = imgSize(2);

imgStack = zeros(H,W,nFrames);

%% --------------------------------------------------------
% Confinement centres

xc = W*rand(N,1);
yc = H*rand(N,1);

%% --------------------------------------------------------
% OU parameters

alpha = exp(-1/tau_c);

% stationary variance = L^2/12
sigma_eq = L/sqrt(12);

% innovation variance
sigma_step = sigma_eq*sqrt(1-alpha^2);

%% --------------------------------------------------------
% Initial positions sampled from equilibrium

x = xc + sigma_eq*randn(N,1);
y = yc + sigma_eq*randn(N,1);

x = max(1,min(W,x));
y = max(1,min(H,y));

%% --------------------------------------------------------

sigma = max(0.5,meanSize + stdSize*randn(N,1));

[X,Y] = meshgrid(1:W,1:H);

traj.x = zeros(N,nFrames);
traj.y = zeros(N,nFrames);

%% =================== Simulation ==========================

for t=1:nFrames

    frame=zeros(H,W);

    bleachFactor = exp(-k_bleach*(t-1));

    for i=1:N

        currentIntensity = intensity*bleachFactor;

        gauss = currentIntensity .* ...
            exp(-((X-x(i)).^2 + (Y-y(i)).^2)/(2*sigma(i)^2));

        frame = frame + gauss;

    end

    frame = frame + noiseLevel*randn(H,W);

    frame(frame<0)=0;

    imgStack(:,:,t)=frame;

    traj.x(:,t)=x;
    traj.y(:,t)=y;

    %% -------- OU dynamics --------

    x = xc + alpha*(x-xc) + sigma_step*randn(N,1);
    y = yc + alpha*(y-yc) + sigma_step*randn(N,1);

    %% image boundaries

    x=max(1,min(W,x));
    y=max(1,min(H,y));

end

%% --------------------------------------------------------

imgStack = imgStack/max(imgStack(:));

%% ================= Compute MSD ===========================

maxLag = nFrames-1;

msd_time = 1:maxLag;

msd_all = zeros(N,maxLag);

for lag=1:maxLag

    dx = traj.x(:,1+lag:end)-traj.x(:,1:end-lag);

    dy = traj.y(:,1+lag:end)-traj.y(:,1:end-lag);

    msd_all(:,lag)=mean(dx.^2+dy.^2,2);

end

msd_mean = mean(msd_all,1);

%% ================= Save stack ============================

if nargin>10 && ~isempty(saveName)

    fprintf('Saving image stack to %s ...\n',saveName);

    for t=1:nFrames

        frame8 = im2uint8(imgStack(:,:,t));

        if t==1
            imwrite(frame8,saveName,'tif','Compression','none');
        else
            imwrite(frame8,saveName,'tif',...
                'WriteMode','append',...
                'Compression','none');
        end

    end

    fprintf('Saved %d frames to %s\n',nFrames,saveName);

end

%% ================= Display ===============================

figure('Name','Confined diffusion simulation','Color','w');

subplot(2,3,1)
imagesc(imgStack(:,:,1))
axis image off
title('Frame 1')
colormap hot

subplot(2,3,2)
imagesc(imgStack(:,:,round(nFrames/2)))
axis image off
title('Mid frame')

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

subplot(2,3,[5 6])

plot(msd_time,msd_all','Color',[0.75 0.75 0.75])

hold on

plot(msd_time,msd_mean,'k','LineWidth',2)

xlabel('Lag (\Delta t) [frames]')
ylabel('MSD [px^2]')

grid on

title('MSD')

legend('Single particles','Mean MSD','Location','northwest')

sgtitle(sprintf('Confined diffusion   L = %.1f px   \\tau_c = %.1f frames',L,tau_c))

%% --------------------------------------------------------

assignin('base','traj',traj);
assignin('base','msd_all',msd_all);
assignin('base','msd_mean',msd_mean);
assignin('base','msd_time',msd_time);

tau = msd_time;

msd_theory = (L^2/3)*(1-exp(-tau/tau_c));

% figure;
% plot(tau,msd_mean,'k','LineWidth',2)
% hold on
% plot(tau,msd_theory,'r--','LineWidth',2)
% legend('Simulation','Theory')
% grid on

end

%%%%%%%%%%%%%%%%%%%%


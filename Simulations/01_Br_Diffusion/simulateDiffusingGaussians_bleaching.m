function [imgStack, traj, msd_all, msd_mean, msd_time] = simulateDiffusingGaussians(...
    N, imgSize, nFrames, D, meanSize, stdSize, intensity, noiseLevel, k_bleach, saveName)

% simulateDiffusingGaussians WITH PHOTOBLEACHING
% k_bleach: photobleaching rate (per frame)

% ------------------- Setup -------------------
rng('shuffle');
H = imgSize(1);
W = imgSize(2);
imgStack = zeros(H, W, nFrames);

x = W * rand(N,1);
y = H * rand(N,1);
sigma = max(0.5, meanSize + stdSize*randn(N,1));
stepStd = sqrt(2*D);

[X, Y] = meshgrid(1:W, 1:H);

% Store trajectories
traj.x = zeros(N, nFrames);
traj.y = zeros(N, nFrames);

% ------------------- Simulation -------------------
for t = 1:nFrames
    frame = zeros(H,W);

    % ?Photobleaching factor (exponential decay)
    bleachFactor = exp(-k_bleach * (t-1));

    for i = 1:N
        % Apply bleaching to intensity
        currentIntensity = intensity * bleachFactor;

        gauss = currentIntensity * exp(-((X - x(i)).^2 + (Y - y(i)).^2) / (2*sigma(i)^2));
        frame = frame + gauss;
    end

    frame = frame + noiseLevel * randn(H,W);
    frame(frame < 0) = 0;
    imgStack(:,:,t) = frame;

    % Store positions
    traj.x(:,t) = x;
    traj.y(:,t) = y;

    % Diffusion step
    x = x + stepStd * randn(N,1);
    y = y + stepStd * randn(N,1);
    x = max(1, min(W, x));
    y = max(1, min(H, y));
end

% Normalize to [0,1]
imgStack = imgStack / max(imgStack(:));

% ------------------- Compute MSD -------------------
maxLag = nFrames - 1;
msd_time = 1:maxLag;
msd_all = zeros(N, maxLag);

for lag = 1:maxLag
    dx = traj.x(:,1+lag:end) - traj.x(:,1:end-lag);
    dy = traj.y(:,1+lag:end) - traj.y(:,1:end-lag);
    msd_all(:,lag) = mean(dx.^2 + dy.^2, 2);
end

msd_mean = mean(msd_all, 1);

% ------------------- Save image stack if requested -------------------
if nargin > 9 && ~isempty(saveName)
    fprintf('Saving image stack to %s ...\n', saveName);
    for t = 1:nFrames
        frame8 = im2uint8(imgStack(:,:,t));
        if t == 1
            imwrite(frame8, saveName, 'tif', 'Compression','none');
        else
            imwrite(frame8, saveName, 'tif', 'WriteMode','append', 'Compression','none');
        end
    end
    fprintf('Saved %d frames to %s\n', nFrames, saveName);
end

% ------------------- Display sample frames -------------------
figure('Name','Diffusing Gaussians Simulation','Color','w');
subplot(2,3,1);
imagesc(imgStack(:,:,1)); axis image off; title('Frame 1'); colormap hot;
subplot(2,3,2);
imagesc(imgStack(:,:,round(nFrames/2))); axis image off; title('Mid frame');
subplot(2,3,3);
imagesc(imgStack(:,:,end)); axis image off; title('Last frame');

subplot(2,3,4);
plot(traj.x', traj.y', '-', 'LineWidth', 1);
axis equal; title('Particle trajectories');
xlabel('x (px)'); ylabel('y (px)');

subplot(2,3,[5 6]);
plot(msd_time, msd_all', '-', 'Color',[0.7 0.7 0.7]); hold on;
plot(msd_time, msd_mean, 'k-', 'LineWidth', 2);
xlabel('Lag (\Delta t) [frames]'); ylabel('MSD [px^2]');
title('MSD per particle and ensemble average');
grid on;
legend('Single particles','Mean MSD','Location','northwest');

%  Show bleaching rate in title
sgtitle(sprintf('%d diffusing Gaussians (D=%.2f, k_{bleach}=%.3f)', N, D, k_bleach));

% ------------------- Assign to workspace -------------------
assignin('base','traj',traj);
assignin('base','msd_all',msd_all);
assignin('base','msd_mean',msd_mean);
assignin('base','msd_time',msd_time);

end
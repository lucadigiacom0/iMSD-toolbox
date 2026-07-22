clear; clc; close all;

%% ---------------- USER INPUT ----------------
dataDir      = 'calibration_tif';   % folder with TIFF stacks
pixelSize_nm = 50;                  % lateral pixel size (nm)
zStep_nm     = 100;                 % z-step (nm) [used for completeness]
satLevel     = 65535;               % 16-bit saturation
% minSNR       = 15;                % bead detection threshold
minSNR       = 50;  
minSigma_px  = 0.8;                 % lower bound (pixels)
maxSigma_px  = 4;                   % upper bound (pixels)
% --------------------------------------------

files = dir(fullfile(dataDir,'*.tif'));

sigmaX = [];
sigmaY = [];
FWHM  = [];

exampleSaved = false;

%% ================= MAIN LOOP =================
for f = 1:numel(files)

    fprintf('Processing %s\n', files(f).name);
    fname = fullfile(files(f).folder, files(f).name);

    info = imfinfo(fname);
    Nz = numel(info);

    % Load stack
    stack = zeros(info(1).Height, info(1).Width, Nz);
    for z = 1:Nz
        stack(:,:,z) = double(imread(fname, z));
    end

    %% Focal plane selection (maximum integrated intensity)
    zScore = squeeze(sum(sum(stack,1),2));
    [~, z0] = max(zScore);
    img = stack(:,:,z0);

    %% Background subtraction
    bg = medfilt2(img,[25 25],'symmetric');
    imgBS = img - bg;
    imgBS(imgBS < 0) = 0;

    %% Bead detection (Difference of Gaussians)
    dog = imgaussfilt(imgBS,1) - imgaussfilt(imgBS,3);
    thresh = mean(dog(:)) + minSNR*std(dog(:));
    peaks = imregionalmax(dog) & dog > thresh;

    [y,x] = find(peaks);

    %% Optional: show detected beads
    figure('Color','w')
    imagesc(img), colormap gray, axis image off
    hold on
    plot(x,y,'ro')
    title(sprintf('Detected beads – %s',files(f).name))

    %% Bead loop
    for i = 1:numel(x)

        xi = x(i);
        yi = y(i);

        rMax = 8; % fixed safe radius (no user tuning)

        if xi<=rMax || yi<=rMax || ...
           xi>=size(img,2)-rMax || yi>=size(img,1)-rMax
            continue
        end

        roi = imgBS(yi-rMax:yi+rMax, xi-rMax:xi+rMax);

        if max(roi(:)) >= satLevel
            continue
        end

        if std(roi(:)) < 1
            continue
        end

        %% 2D Gaussian fit
        [X,Y] = meshgrid(1:size(roi,2),1:size(roi,1));
        XY = [X(:),Y(:)];
        Z  = roi(:);

        A0 = max(Z);
        x0 = (size(roi,2)+1)/2;
        y0 = (size(roi,1)+1)/2;
        s0 = 2;
        B0 = median(Z);

        p0 = [A0,x0,y0,s0,s0,B0];

        lb = [0, x0-1, y0-1, minSigma_px, minSigma_px, 0];
        ub = [Inf, x0+1, y0+1, maxSigma_px, maxSigma_px, max(Z)];

        gauss2D = @(p,XY) ...
            p(1)*exp(-((XY(:,1)-p(2)).^2/(2*p(4)^2) + ...
                       (XY(:,2)-p(3)).^2/(2*p(5)^2))) + p(6);

        opts = optimoptions('lsqcurvefit','Display','off');

        try
            p = lsqcurvefit(gauss2D,p0,XY,Z,lb,ub,opts);
        catch
            continue
        end

        sx = p(4);
        sy = p(5);

        % Shape filter
        if sx/sy > 1.5 || sy/sx > 1.5
            continue
        end

        sigmaX(end+1) = sx;
        sigmaY(end+1) = sy;

        fwhm_nm = 2.355 * mean([sx sy]) * pixelSize_nm;
        FWHM(end+1) = fwhm_nm;

        %% Save one example bead fit
        if ~exampleSaved
            exampleSaved = true;

            Zfit = reshape(gauss2D(p,XY), size(roi));

            figure('Color','w','Position',[100 100 900 300])
            subplot(1,3,1)
            imagesc(roi), axis image off
            title('Bead ROI')

            subplot(1,3,2)
            imagesc(Zfit), axis image off
            title('Gaussian fit')

            subplot(1,3,3)
            imagesc(roi - Zfit), axis image off
            title('Residuals')

            sgtitle(sprintf('Example bead | FWHM = %.1f nm',fwhm_nm))
        end
    end
end

%% ================= RESULTS =================

% Convert to nm
sigmaX_nm = sigmaX * pixelSize_nm;
sigmaY_nm = sigmaY * pixelSize_nm;

% PSF variance (sigma0^2)
sigma0_sq = mean((sigmaX_nm.^2 + sigmaY_nm.^2)/2);

% Mean FWHM
meanFWHM = mean(FWHM);
stdFWHM  = std(FWHM);

fprintf('\nMean FWHM = %.1f nm\n',meanFWHM);
fprintf('Std  FWHM = %.1f nm\n',stdFWHM);
fprintf('Sigma0^2  = %.2e nm^2\n',sigma0_sq);

%% ================= FIGURES =================

% FWHM distribution
figure('Color','w')
histogram(FWHM,20,'Normalization','pdf')
xlabel('FWHM (nm)')
ylabel('Probability density')
title({ ...
    sprintf('PSF FWHM distribution'), ...
    sprintf('Mean FWHM = %.1f nm | \\sigma_0^2 = %.2e nm^2', ...
            meanFWHM, sigma0_sq) })
grid on

% FWHM per bead
figure('Color','w')
plot(FWHM,'ko','MarkerFaceColor','k')
xlabel('Bead index')
ylabel('FWHM (nm)')
title('FWHM per bead')
grid on

%% ================= EXPORT =================

T = table( ...
    sigmaX_nm(:), ...
    sigmaY_nm(:), ...
    FWHM(:), ...
    'VariableNames',{'SigmaX_nm','SigmaY_nm','FWHM_nm'});

writetable(T,'PSF_beads.csv');

summary = table( ...
    meanFWHM, stdFWHM, sigma0_sq, ...
    'VariableNames',{'MeanFWHM_nm','StdFWHM_nm','Sigma0_sq_nm2'});

writetable(summary,'PSF_summary.csv');

disp('PSF calibration completed successfully.');

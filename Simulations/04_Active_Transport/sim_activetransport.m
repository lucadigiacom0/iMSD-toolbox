clear all
clear all
clc

%% ============================================================
% Acquisition parameters

px_size    = 0.05;      % µm/pixel
frame_time = 0.075;     % s/frame

%% ============================================================
% Physical parameters

Din  = 0.0005;            % µm²/s
Dout = 0.05;            % µm²/s

v = 0.050;               % µm/s

%% ============================================================
% Geometry of guided paths

pathWidth  = 0.20;      % µm
pathLength = 4.00;      % µm

Npaths = 100;

%% ============================================================
% Imaging parameters

sigma0     = 0.09;      % µm
sigma0_err = 0.015;     % µm

N_spot   = 100;
H        = 256;
n_frames = 1000;

intens   = 160;
noise_l  = 8;

%% ============================================================
% Bleaching

tau_bl = 1200;
k_bl   = 1/tau_bl;

%% ============================================================

saveName = 'sim_guided_diffusion.tif';

%% ============================================================
% Unit conversion

Din_px  = Din /(px_size^2)*frame_time;
Dout_px = Dout/(px_size^2)*frame_time;

v_px = v/px_size*frame_time;

pathWidth_px  = pathWidth /px_size;
pathLength_px = pathLength/px_size;

sigma_px     = sigma0/px_size;
sigmaErr_px  = sigma0_err/px_size;

%% ============================================================

simulateGuidedDiffusionGaussians_bleaching( ...
    N_spot,...
    [H H],...
    n_frames,...
    Din_px,...
    Dout_px,...
    v_px,...
    pathWidth_px,...
    pathLength_px,...
    Npaths,...
    sigma_px,...
    sigmaErr_px,...
    intens,...
    noise_l,...
    k_bl,...
    saveName);
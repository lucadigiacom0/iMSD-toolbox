clear all
clear all
clc

%% ------------------------------------------------------------
% Acquisition parameters

px_size    = 0.05;      % µm/pixel
frame_time = 0.075;     % s/frame

%% ------------------------------------------------------------
% Physical parameters

D = 3e-3;             % µm^2/s
L = 0.300;               % µm (compartment side)
P_cross = 0.01;         % probability of crossing a barrier

sigma0 = 0.09;          % µm
sigma0_err = 0.015;     % µm

%% ------------------------------------------------------------
% Simulation parameters

N_spot   = 150;
H         = 256;
n_frames  = 1000;

intens   = 160;
noise_l  = 8;

%% ------------------------------------------------------------
% Photobleaching

tau_bl = 1000;
k_bl   = 1/tau_bl;

%% ------------------------------------------------------------

saveName = 'sim_transient_confinement.tif';

%% ------------------------------------------------------------
% Run simulation

simulateTransientConfinementGaussians_bleaching( ...
    N_spot,...
    [H H],...
    n_frames,...
    D/(px_size^2)*frame_time,...
    L/px_size,...
    P_cross,...
    sigma0/px_size,...
    sigma0_err/px_size,...
    intens,...
    noise_l,...
    k_bl,...
    saveName);
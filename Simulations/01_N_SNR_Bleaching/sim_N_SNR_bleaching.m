clear all
clear all
clc

px_size=0.05;
frame_time=0.075;

D=3.5*10^-3;
sigma0=0.09;
sigma0_err=0.015;

H=256;
N_spot=150;
n_frames=1000;

intens=160;
noise_l=8;
tau_bl=800;
k_bl=1/tau_bl;

saveName='sim_SNR_bleaching_01.tif';

simulateDiffusingGaussians_bleaching(N_spot, [H,H], n_frames, D/(px_size^2)*frame_time, sigma0/px_size,...
    sigma0_err/px_size, intens, noise_l, k_bl, saveName);
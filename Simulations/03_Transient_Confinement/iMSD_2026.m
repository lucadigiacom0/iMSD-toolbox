clear all
clc
close all

tic

filename='sim_transient_confinement.tif';
export_value = 0;

disp('     importing image-stack')
imsd_01_import
A = double(imgStack);

%%% input parameters (edit for your acquisition)
px_size    = 0.05;                 %%% pixel size (micron)
frame_time = 0.075;                %%% time between successive frames (s)
PSF_waist  = 0.25;                 %%% PSF waist from calibration
N = size(A,3);
P = size(A,1);
R = px_size*P;

SpatialLimit = 24;                 %%% linear extent of spatial domain (px, even)
TauLimit     = int64(N/8)+1;       %%% lag-time domain (default N/10)
fitt2comp = false;                 %%% true = two-component fit (two iMSD curves)

time_vec     = frame_time*(1:double(TauLimit))';

imsd_02_spatial
imsd_03_spatiotemporal
imsd_04_output

toc
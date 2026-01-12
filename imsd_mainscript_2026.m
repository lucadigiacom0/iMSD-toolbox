clear all
clear all
clc
close all

filename = 'test_imagestack.tif';
disp('     importing image-stack')

export_value=1;

imsd_01_import
A=double(imgStack);

%%% input parameters
px_size= 0.04824;                  %%% pixel size (micron)
frame_time=0.192;                  %%% time between two successive frames (s)

%%%%
N=size(A,3);                       %%% number of frames
P=size(A,1);                       %%% lateral size (px)
R=px_size*P;                       %%% lateral size (micron)

%%% spatial-temporal domain
SpatialLimit=32;                   %%% linear extent of spatial domain (px)
                                   %%% only even values
TauLimit=int64(N/5)+1;          %%% lag-time domain (default N/10)

time_vec=frame_time*(1:double(TauLimit))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imsd_02_spatial 
imsd_03_spatiotemporal
imsd_04_output
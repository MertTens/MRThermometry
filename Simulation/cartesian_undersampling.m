function [shifted_usmat] = cartesian_undersampling(dimx, dimy, num_PE, time_pts)
%
% This file takes in an image, and returns the undersampled Fourier data
% according to radial undersampling with a given number of spokes and

% at a given start angle
param.n   = num_PE;   % readouts per frame--can be changed retrospectively
param.FR  = time_pts;  % Number of frames--can be changed retrospectively
param.PE  = dimx;  % Size of PE grid
param.E   = 1;   % Number of encoding, E=1 for cine, E=2 for flow (phase-contrast MRI)
param.ir  = 1;   % ir = 1 or 2 for golden angle, ir > 2 for tiny golden angles; default value: 1
param.k   = 3;   % k>=1. k=1 uniform; k>1 variable density profile; larger k means flatter top (default: 3)
param.s   = 2;   % s>=0; % largers s means higher sampling density in the middle (default: 2, range: 0-10, precision: 0.1)
param.dsp = 0;   % Display patter after reach sample


[samp, PEInd] = cava_fun(param);
t = reshape(samp,1,dimx,time_pts);
t = repmat(t,dimy,1,1);
usmat = t;
shifted_usmat = usmat;
% shifted_usmat = circshift(usmat, dimx/2, 1);
% shifted_usmat = circshift(shifted_usmat, dimy/2, 2);
close all
end
%================================================================================
% GB24 - Generate initial pressure
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex22_vesselsGB;

close all;
clear all;

% Define Domain
Nx = 128;         % number of grid points in the x (row) direction
Ny = 128;         % number of grid points in the x (row) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

%==============================
% Build initial pressure
%==============================
% Veins phantom
inputIm = imread('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_2DRT/Phantoms/Veins_modified.jpg');
u0 = double(255-inputIm)/255;
u0 = u0(:, 31:Nx+30);
u0(1:15,:)        = 0;
u0(end-15:end, :) = 0;
u0(:, end-15:end) = 0;
u0(:, 1:15)       = 0;

% kWave filtering
dt = 1e-8;
tMax = 1e-5;
kgrid.t_array = 0:dt:tMax;
u0 = smooth(kgrid, u0, true);

% Filter 
binMat = zeros(Nx, Ny);
factorFilter = 1/8;
binMat(:) = 0;
binMat(floor(Nx*(1/2-factorFilter)):floor(Nx+(1/2+factorFilter)), floor(Ny*(1/2-factorFilter)):2*floor(Ny/3)) = 1;
fft_u0    = fft2(u0);
fft_u0    = fftshift(fft_u0);
fft_u0_LF = fft_u0.*binMat;
fft_u0_HF = fft_u0.*(1-binMat);

% Low frequency
fft_u0_LF = fftshift(fft_u0_LF);
fft_u0_HF = fftshift(fft_u0_HF);
u0_LF = real(ifft2(fft_u0_LF));
u0_HF = real(ifft2(fft_u0_HF));

%=========================================================================
% PLOT
%=========================================================================
figure;
imagesc(real(fft_u0));
title('Frequency Decomposition');

figure;
imagesc(u0);
title('Initial Pressure');

figure;
imagesc(u0_LF);
title('Initial Pressure - LF');

figure;
imagesc(u0_HF);
title('Initial Pressure - HF');

figure;
imagesc(u0_LF + u0_HF);
title('Initial Pressure - LF + HF');

% Save
save initial_pressure_128x128_broadBand u0;
save initial_pressure_128x128_LF u0_LF;
save initial_pressure_128x128_HF u0_HF;

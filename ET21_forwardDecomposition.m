%================================================================================
% ET21
% Formula (35) of "Fast Multiscale Gaussian Wavepacket Transforms" by Qian and Ying
% Decompose given initial pressure in Gaussian Beams
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex21_singleGB;

close all;
clear all;


%==================================================
%=======    PHANTOM          ======================
%==================================================
%load sensor_data_singleGB_GB;
load sensor_data_singleGB_kWave;
u0 = real(sensor_data);

%==================================================
%=======    DOMAIN           ======================
%==================================================
% Define Domain
Nx = size(u0, 1);
Ny = size(u0, 2);
N_boundF = 0;
dx = 1e-4;
dy = 3e-5; % dt*c

%===============================================================================================
%=======                     ===================================================================
%=======    BOXES            ===================================================================
%=======                     ===================================================================
%===============================================================================================
gb1 = GB(Nx, dx, Ny, dy);
x_axis = gb1.x_axis;
y_axis = gb1.y_axis;
[X, Y] = meshgrid(y_axis, x_axis);

% Number of coronae
L = 2;

% Box scales
gb1.parametrize_domain(L);

%===============================================================================================
%=======                          ==============================================================
%=======    FREQUENCY ANALYSIS    ==============================================================
%=======                          ==============================================================
%===============================================================================================
% Plot Initial Pressure
figure;
surf(X, Y, u0, 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
title('Initial Pressure to decompose - HF')
  
%========================================
% FFT (MANUAL)
%========================================
gb1.setU(u0);
gb1.fft_resize(N_boundF);
fU0_manual = gb1.fft_u0;

figure;
surf(gb1.fy_axis, gb1.fx_axis, real(fU0_manual), 'EdgeColor', 'none');
title('Real fft u0');
colorbar();
view(2);
figure;
surf(gb1.fy_axis, gb1.fx_axis, imag(fU0_manual), 'EdgeColor', 'none');
title('Imag fft u0');
colorbar();
view(2);


%========================================
% FILTERS
%========================================
gb1.filters();

%========================================
% COMPUTE F
%========================================
GBF = gb1.coefficients();
gb1.plot_coeffs()

% Plot
figure;
surf(y_axis, x_axis, real(GBF)/max(real(GBF(:))), 'EdgeColor', 'none');
%surf(x_axis, y_axis, real(GBF), 'EdgeColor', 'none');
view(2);
axis tight;
title('GB real - HF');
colorbar();
%saveas(gcf, 'GB21_GBapprox_real_HF', 'png');

figure;
surf(y_axis, x_axis, imag(GBF), 'EdgeColor', 'none');
view(2)
axis tight;
title('GB imag')
colorbar();

figure;
surf(y_axis, x_axis, real(u0)/max(real(u0(:)))-real(GBF)/max(real(GBF(:))), 'EdgeColor', 'none');
view(2);
axis tight;
title('Difference Approximation - HF')
colorbar();

%========================================
% SAVE PARAMETERS
%========================================
sensor_data_GBdecomposition = GBF;
save sensor_data_singleGB_GBparameters.mat sensor_data_GBdecomposition gb1;

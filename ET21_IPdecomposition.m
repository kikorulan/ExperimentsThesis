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
load singleGB;
u0 = singleGB;

%==================================================
%=======    DOMAIN           ======================
%==================================================
% Define Domain
Nx = size(u0, 1);
Ny = size(u0, 2);
N_boundF = 0;
dx = 1e-4;
dy = 1e-4;

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
surf(X, Y, real(u0), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
title('Initial Pressure to decompose - Single GB')
  
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

normGB = max(real(GBF(:)));
normIP = max(real(u0(:)));
% Plot
figure;
surf(y_axis, x_axis, real(GBF)/normGB, 'EdgeColor', 'none');
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
surf(y_axis, x_axis, real(u0)/normIP-real(GBF)/normGB, 'EdgeColor', 'none');
view(2);
axis tight;
title('Difference Approximation - HF')
colorbar();

singleGB_GBdecomposition = GBF;
save singleGB_GBdecomposition singleGB_GBdecomposition gb1;


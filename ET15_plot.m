% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex15_synth_recon3D_het;
cd /scratch0/NOT_BACKED_UP/frullan/ExperimentsThesis/Ex15_synth_recon3D_het;

clear all;
close all;

% Functions
[TV, D, DTV] = TVOperators(3, 'none');
norm_distance = @(x, y) sum((x(:) - y(:)).*(x(:) - y(:)));
obj_data = @(y0, y) 0.5*norm_distance(y0, y);
obj_reg  = @(lambda, u0) lambda*TV(u0);
obj_function = @(y0, y, lambda, u0) obj_data(y0, y) + obj_reg(lambda, u0);

%==================================================
% Dimensions
%==================================================
% Import dimensions 
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
% Objective value of function
OBJ_VAL = 0.93475;
rel_distance = @(curve, iter) (curve(iter) - OBJ_VAL)/(curve(1) - OBJ_VAL);

% LOAD DATA
load ./results/error_vectors/GD_error_lambda1em4;
%load ./results/error_vectors/SGD_error_lambda1em4_batch1800;
%load ./results/error_vectors/FISTA_error_lambda1em4;
%load ./results/error_vectors/PDHG_error_lambda1em4_sigma5em1;
%load ./results/error_vectors/SPDHG_error_lambda1em4_sigma1em1_batch100;

% Save results
saveResults = 0;

% Choose index
GD.plotIndex    = 2;
SGD.plotIndex   = 2;
FISTA.plotIndex = 2;
PDHG.plotIndex  = 2;
SPDHG.plotIndex = 3;

% Choose subiter
GD.subiter    = 18;
SGD.subiter   = 30;
FISTA.subiter = 16;
PDHG.subiter  = 19;
SPDHG.subiter = 8;
lim_iter = 0;

% Reconstruction plots
plot_reconGD    = 0;
plot_reconSGD   = 0;
plot_reconFISTA = 0;
plot_reconPDHG  = 0;
plot_reconSPDHG = 0;

% Auxiliary plots
plot_auxGD    = 1;
plot_auxSGD   = 0;
plot_auxFISTA = 0;
plot_auxPDHG  = 0;
plot_auxSPDHG = 0;

% Converge plots
plot_PSNR    = 0;
plot_primalE = 0;
plot_dualE   = 0;
plot_relDE   = 0;
plot_dataE   = 0;
plot_regE    = 0;
plot_relObj  = 0;

%========================================================================================================================
% PRIMAL AND DUAL DATA
%========================================================================================================================
param.dx = dx;

% Forward signal - full data
%%  time_signal = importdata(['./input_data/forwardSignal_reference_14400sensors.dat'], ' ', 0);
%%  y0 = time_signal(2:end, :);

% Forward signal - subsampled
time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
t_array = 1e6*time_signal(1, :);
y0 = time_signal(2:end, :);
figure;
imagesc(y0(1:60, :));
xticklabels = 0:2:8;
xticks = linspace(1, find(t_array > 8, 1), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
colorbar();
set(gca,'FontSize',15);
xlabel('t [\mus]')
if(saveResults)
    saveas(gcf, './figures/ET15_forwardSignal_noisy5', 'epsc');
end
  
% Load Initial pressure
pressure_adjoint = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
plot_projection_compact(pressure_adjoint, param);
if(saveResults)
    saveas(gcf, './figures/ET15_initialPressure', 'epsc');
end

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pressure_adjoint_RT.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
plot_projection_compact(pressure_adjoint, param);
if(saveResults)
    saveas(gcf, './figures/ET15_pixelPressure_adjoint', 'epsc');
end

%========================================================================================================================
% RECONSTRUCTION PLOTS
%========================================================================================================================

% Gradient Descent
if(plot_reconGD)
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{GD.plotIndex}, '_lambda', GD.lambda, '_iter', int2str(GD.subiter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, param);
if(saveResults)
    saveas(gcf, ['./figures/ET15_GD'], 'epsc');
end
end

% Stochastic Gradient Descent
if(plot_reconSGD)
pixelPressureMatrix = importdata(['./results/adjoint/SFB/pixelPressure_S-GD_tau', SGD.tau{SGD.plotIndex}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(SGD.subiter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, param);
if(saveResults)
    saveas(gcf, ['./figures/ET15_S-GD'], 'epsc');
end
end

% FISTA
if(plot_reconFISTA)
pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau{FISTA.plotIndex}, '_lambda', FISTA.lambda, '_iter', int2str(FISTA.subiter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, param);
if(saveResults)
    saveas(gcf, ['./figures/ET15_FISTA'], 'epsc');
end
end

% PDHG
if(plot_reconPDHG)
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{PDHG.plotIndex}, '_theta1_lambda', PDHG.lambda, '_iter', int2str(PDHG.subiter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, param);
if(saveResults)
    saveas(gcf, ['./figures/ET15_PDHG'], 'epsc');
end
end

% S-PDHG
if(plot_reconSPDHG)
pixelPressureMatrix = importdata(['./results/adjoint/SPDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{SPDHG.plotIndex}, '_theta1_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(SPDHG.subiter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, param);
if(saveResults)
    saveas(gcf, ['./figures/ET15_S-PDHG'], 'epsc');
end
end


%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================

%======================================================================
% Gradient descent
%======================================================================
if(plot_auxGD)
disp('GD');
L = length(GD.tau);
% Plot PSNR
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_psnr{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('GD PSNR')
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('GD primal')
% Plot dual
figure();
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
title('GD dual')
end

%======================================================================
% Stochastic Gradient descent
%======================================================================
if(plot_auxSGD)
disp('SGD');
L = length(SGD.tau);
% Plot PSNR
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}-1, SGD_error_psnr{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('SGD PSNR')
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}-1, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('SGD primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}-1, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('SGD dual')
end

%======================================================================
% FISTA
%======================================================================
if(plot_auxFISTA)
disp('FISTA');
L = length(FISTA.tau);
% Plot PSNR
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:FISTA.nIter{ii}-1, FISTA_error_psnr{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('FISTA PSNR')
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:FISTA.nIter{ii}-1, FISTA_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('FISTA primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:FISTA.nIter{ii}-1, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('FISTA dual')
end

%======================================================================
% PDHG
%======================================================================
if(plot_auxPDHG)
disp('PDHG');
L = length(PDHG.tau);
% Plot psnr
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:PDHG.nIter{ii}-1, PDHG_error_psnr{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('PDHG PSNR')
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:PDHG.nIter{ii}-1, PDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('PDHG primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:PDHG.nIter{ii}-1, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('PDHG dual')
end

%==============================
% SPDHG
%==============================
if(plot_auxSPDHG)
disp('SPDHG');
L = length(SPDHG.tau);
% Plot PSNR
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SPDHG.nIter{ii}-1, SPDHG_error_psnr{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('SPDHG PSNR')
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SPDHG.nIter{ii}-1, SPDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('SPDHG primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SPDHG.nIter{ii}-1, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('SPDHG dual')
end

%========================================================================================================================
% PLOT ALL
%========================================================================================================================

%======================================================================
% PLOT PSNR
%======================================================================
if(plot_PSNR)
figure();
if (lim_iter)
    semilogy(0:GD.subiter-1, GD_error_psnr{GD.plotIndex}(1:GD.subiter), 'Color', 'r', 'Linewidth', 1.5);
    hold on;       
    semilogy(0:SGD.subiter-1, SGD_error_psnr{SGD.plotIndex}(1:SGD.subiter), 'Color', 'g', 'Linewidth', 1.5);
    semilogy(0:FISTA.subiter-1, FISTA_error_psnr{FISTA.plotIndex}(1:FISTA.subiter), 'Color', 'b', 'Linewidth', 1.5);
    semilogy(0:PDHG.subiter-1, PDHG_error_psnr{PDHG.plotIndex}(1:PDHG.subiter), 'Color', 'm', 'Linewidth', 1.5);
    semilogy(0:SPDHG.subiter-1, SPDHG_error_psnr{SPDHG.plotIndex}(1:SPDHG.subiter), 'Color', 'c', 'Linewidth', 1.5);
else
    semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_psnr{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
    hold on;       
    semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_psnr{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
    semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_psnr{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
    semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_psnr{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
    semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_psnr{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
end
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 100 46 50]);
xlabel('iter/epoch');
ylabel('PSNR');
set(gca,'FontSize',15);
if(saveResults)
    saveas(gcf, './figures/ET15_PSNR', 'epsc');
end
end

%======================================================================
% PLOT PRIMAL
%======================================================================
if(plot_primalE)
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_pd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_pd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_pd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_pd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_pd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 40 1e4]);
xlabel('iter/epoch');
ylabel('primal error');
set(gca,'FontSize',15);
if(saveResults)
    saveas(gcf, './figures/ET15_pe', 'epsc');
end
end

%======================================================================
% PLOT DUAL AND RELATIVE
%======================================================================
if(plot_dualE)
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_dd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_dd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_dd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_dd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_dd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 9 16]);
xlabel('iter/epoch');
ylabel('dual error');
%set(gca,'FontSize',15);
if(saveResults)
    saveas(gcf, './figures/ET15_de', 'epsc');
end
end

if(plot_relDE)
% Relative error
figure();
semilogy(1:GD.nIter{GD.plotIndex}-1, (GD_error_dd{GD.plotIndex}(1:end-1)-GD_error_dd{GD.plotIndex}(2:end))./GD_error_dd{GD.plotIndex}(2:end), 'Color', 'r', 'Linewidth', 1.5)
hold on;       
semilogy(1:SGD.nIter{SGD.plotIndex}-1, (SGD_error_dd{SGD.plotIndex}(1:end-1)-SGD_error_dd{SGD.plotIndex}(2:end))./SGD_error_dd{SGD.plotIndex}(2:end), 'Color', 'g', 'Linewidth', 1.5)
semilogy(1:FISTA.nIter{FISTA.plotIndex}-1, (FISTA_error_dd{FISTA.plotIndex}(1:end-1)-FISTA_error_dd{FISTA.plotIndex}(2:end))./FISTA_error_dd{FISTA.plotIndex}(2:end), 'Color', 'b', 'Linewidth', 1.5)
semilogy(1:PDHG.nIter{PDHG.plotIndex}-1, (PDHG_error_dd{PDHG.plotIndex}(1:end-1)-PDHG_error_dd{PDHG.plotIndex}(2:end))./PDHG_error_dd{PDHG.plotIndex}(2:end), 'Color', 'm', 'Linewidth', 1.5)
semilogy(1:SPDHG.nIter{SPDHG.plotIndex}-1, (SPDHG_error_dd{SPDHG.plotIndex}(1:end-1)-SPDHG_error_dd{SPDHG.plotIndex}(2:end))./SPDHG_error_dd{SPDHG.plotIndex}(2:end), 'Color', 'c', 'Linewidth', 1.5)
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([1 30 5e-3 1]);
xlabel('iter/epoch');
ylabel('relative error');
%set(gca,'FontSize',15);
if(saveResult)
    saveas(gcf, './figures/ET15_relativeDualError', 'epsc');
end
end

%======================================================================
% PLOT REGULARIZATION AND DATA FIT
%======================================================================
% Data term
if(plot_dataE)
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_data{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_data{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_data{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_data{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_data{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('data term');
end

% Plot Regularization Term
if(plot_regE)
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_reg{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_reg{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_reg{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_reg{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_reg{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('regularization');
end

%======================================================================
% RELATIVE DISTANCE TO OBJECTIVE FUNCTION
%======================================================================
if (plot_relObj)
figure();
if (lim_iter)
    semilogy(0:GD.subiter-1, rel_distance(GD_error_dd{GD.plotIndex}, 1:GD.subiter), 'Color', 'r', 'Linewidth', 1.5);
    hold on;       
    semilogy(0:SGD.subiter-1, rel_distance(SGD_error_dd{SGD.plotIndex}, 1:SGD.subiter), 'Color', 'g', 'Linewidth', 1.5);
    semilogy(0:FISTA.subiter-1, rel_distance(FISTA_error_dd{SGD.plotIndex}, 1:FISTA.subiter), 'Color', 'b', 'Linewidth', 1.5);
    semilogy(0:PDHG.subiter-1, rel_distance(PDHG_error_dd{PDHG.plotIndex}, 1:PDHG.subiter), 'Color', 'm', 'Linewidth', 1.5);
    semilogy(0:SPDHG.subiter-1, rel_distance(SPDHG_error_dd{SPDHG.plotIndex}, 1:SPDHG.subiter), 'Color', 'c', 'Linewidth', 1.5);
else
    semilogy(0:GD.nIter{GD.plotIndex}-1, rel_distance(GD_error_dd{GD.plotIndex}, 1:GD.nIter{GD.plotIndex}), 'Color', 'r', 'Linewidth', 1.5);
    hold on;       
    semilogy(0:SGD.nIter{SGD.plotIndex}-1, rel_distance(SGD_error_dd{SGD.plotIndex}, 1:SGD.nIter{SGD.plotIndex}), 'Color', 'g', 'Linewidth', 1.5);
    semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, rel_distance(FISTA_error_dd{SGD.plotIndex}, 1:FISTA.nIter{FISTA.plotIndex}), 'Color', 'b', 'Linewidth', 1.5);
    semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, rel_distance(PDHG_error_dd{PDHG.plotIndex}, 1:PDHG.nIter{PDHG.plotIndex}), 'Color', 'm', 'Linewidth', 1.5);
    semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, rel_distance(SPDHG_error_dd{SPDHG.plotIndex}, 1:SPDHG.nIter{SPDHG.plotIndex}), 'Color', 'c', 'Linewidth', 1.5);
end
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 20 1e-4 1]);
xlabel('iter/epoch');
ylabel('relative distance to objective');
set(gca,'FontSize',15);
if(saveResults)
    saveas(gcf, './figures/ET15_relObjective', 'epsc');
end
end



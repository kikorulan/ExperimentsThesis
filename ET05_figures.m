cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex05_simulation3D_het;

clear all;
close all;

saveFigures = 0;
fontSize = 14;
%==================================================
% Dimensions
%==================================================
% Import dimensions 
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);
   
%==================================================
% Load data
%==================================================
% Forward data
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240_smooth.dat', ' ', 0);
sensor_data_kWave = importdata(['./input_data/forwardSignal_kWave.dat'], ' ', 0);
sensor_data_RT    = importdata(['./input_data/forwardSignal_RT.dat'], ' ', 0);

% Pressure
pressure_adjoint_kWave_matrix         = importdata(['./output_data/pressure_adjoint_kWave.dat'], ' ', 0);
pressure_adjoint_kWave_RT_data_matrix = importdata(['./output_data/pressure_adjoint_kWave_RT_data.dat'], ' ', 0);
pressure_adjoint_RT_matrix            = importdata(['./output_data/pressure_adjoint_RT.dat'], ' ', 0);
pressure_adjoint_RT_kWave_data_matrix = importdata(['./output_data/pressure_adjoint_RT_kWave_data.dat'], ' ', 0);

%==================================================
% Sound speed
%==================================================
param.dx = dx;
param.lim = true;
param.cMin = 1580*(1-0.06);
param.cMax = 1580*(1+0.06);
param.fontSize = fontSize;
% Sound speed
sound_speed = importdata(['./input_data/sound_speed.dat'], ' ', 0);
sound_speed = matrix2cube(sound_speed, Nz);
plot_slice_compact(sound_speed, param);
if saveFigures
saveas(gcf, 'Example05_soundSpeed', 'png');
saveas(gcf, 'Example05_soundSpeed.fig');
end

%========================================================================================================================
% INITIAL PRESSURE AND FORWARD SIGNAL
%========================================================================================================================
position     = [700 700 530 630];
positionY    = [700 700 550 630];
positionBar  = [700 700 620 630];
positionYBar = [700 700 620 630];

% kWave
figure;
surf(1e6*sensor_data_RT(1, :), 1:120, sensor_data_kWave(2:121, :), 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('t [\mus]');
ylabel('Sensor');
caxis([-0.073 0.073])
set(gca,'FontSize',16);
set(gcf, 'pos', positionY);
if saveFigures
saveas(gcf, 'Example05_kWave_f2D', 'png');
saveas(gcf, 'Example05_kWave_f2D.fig');
end

% RT
figure;
surf(1e6*sensor_data_RT(1, :), 1:120, sensor_data_RT(2:121, :), 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('t [\mus]');
caxis([-0.073 0.073])
colorbar();
set(gca,'FontSize',16);
set(gcf, 'pos', positionBar);
if saveFigures
saveas(gcf, 'Example05_RT_f2D', 'png');
saveas(gcf, 'Example05_RT_f2D.fig');
end


% RT
figure;
surf(1e6*sensor_data_RT(1, :), 1:120, sensor_data_RT(2:121, :) - sensor_data_kWave(2:121, :), 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('t [\mus]');
caxis([-0.0146 0.0146])
colorbar();
set(gca,'FontSize',16);
set(gcf, 'pos', positionBar);
if saveFigures
saveas(gcf, 'Example05_error_RT_f2D', 'png');
saveas(gcf, 'Example05_error_RT_f2D.fig');
end


ERT = sum(sensor_data_RT(:).^2)
EKW = sum(sensor_data_kWave(:).^2)
error = sensor_data_RT - sensor_data_kWave;
REE = sum(error(:).^2)/EKW
%========================================================================================================================
% ADJOINT PRESSURE
%========================================================================================================================
param.dx = dx;
param.lim = true;
param.cMin = 0;
param.cMax = 0.2;
param.fontSize = fontSize;
% Adjoint kWave
param.cb = false;
pressure_adjoint_kWave = matrix2cube(pressure_adjoint_kWave_matrix, Nz);
plot_projection_compact(pressure_adjoint_kWave, param);
if saveFigures
saveas(gcf, 'Example05_kWave_adjoint', 'png');
saveas(gcf, 'Example05_kWave_adjoint.fig');
end
% Adjoint RT
param.cb = true;
pressure_adjoint_RT = matrix2cube(pressure_adjoint_RT_matrix, Nz);
plot_projection_compact(pressure_adjoint_RT, param);
if saveFigures
saveas(gcf, 'Example05_RT_adjoint', 'png');
saveas(gcf, 'Example05_RT_adjoint.fig');
end
% Adjoint kWave - RT data
param.cb = false;
pressure_adjoint_kWave_RTdata = matrix2cube(pressure_adjoint_kWave_RT_data_matrix, Nz);
plot_projection_compact(pressure_adjoint_kWave_RTdata, param);
if saveFigures
saveas(gcf, 'Example05_kWave_adjoint_RT_data', 'png');
saveas(gcf, 'Example05_kWave_adjoint_RT_data.fig');
end
% Adjoint RT - kWave data
param.cb = true;
pressure_adjoint_RT_kWaveData = matrix2cube(pressure_adjoint_RT_kWave_data_matrix, Nz);
plot_projection_compact(pressure_adjoint_RT_kWaveData, param);
if saveFigures
saveas(gcf, 'Example05_RT_adjoint_kWave_data', 'png');
saveas(gcf, 'Example05_RT_adjoint_kWave_data.fig');
end

EKW = sum(pressure_adjoint_kWave(:).^2)
ERT = sum(pressure_adjoint_RT_kWaveData(:).^2)

error_adj = pressure_adjoint_kWave - pressure_adjoint_RT;
error_kW  = pressure_adjoint_kWave - pressure_adjoint_kWave_RTdata;
error_RT  = pressure_adjoint_kWave - pressure_adjoint_RT_kWaveData;

REE_adj = sum(error_adj(:).^2)/sum(pressure_adjoint_kWave(:).^2)
REE_kW  = sum(error_kW(:).^2)/sum(pressure_adjoint_kWave(:).^2)
REE_RT  = sum(error_RT(:).^2)/sum(pressure_adjoint_kWave(:).^2)

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex04_simulation3D_homo;

clear all;
close all;

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
pressure_adjoint_kWave_matrix        = importdata(['./output_data/pressure_adjoint_kWave.dat'], ' ', 0);
pressure_adjoint_kWave_RTdata_matrix = importdata(['./output_data/pressure_adjoint_kWave_RTdata.dat'], ' ', 0);
pressure_adjoint_RT_matrix           = importdata(['./output_data/pressure_adjoint_RT.dat'], ' ', 0);
%pressure_adjoint_RT_kWaveData_matrix = importdata(['./output_data/pressure_adjoint_RT_kWaveData.dat'], ' ', 0);

%========================================================================================================================
% INITIAL PRESSURE AND FORWARD SIGNAL
%========================================================================================================================
% Initial Pressure
u0 = matrix2cube(u0Matrix, Nz);
plot_projection_compact(u0, dx);

% kWave
figure;
imagesc(sensor_data_kWave(2:81, :));
colorbar();

% RT
figure;
imagesc(sensor_data_RT(2:81, :));
colorbar();

%========================================================================================================================
% ADJOINT PRESSURE
%========================================================================================================================
% Adjoint kWave
pressure_adjoint_kWave = matrix2cube(pressure_adjoint_kWave_matrix, Nz);
plot_projection_compact(pressure_adjoint_kWave, dx);

% Adjoint RT
pressure_adjoint_RT = matrix2cube(pressure_adjoint_RT_matrix, Nz);
plot_projection_compact(pressure_adjoint_RT, dx);

% Adjoint kWave - RT data
pressure_adjoint_kWave_RTdata = matrix2cube(pressure_adjoint_kWave_RTdata_matrix, Nz);
plot_projection_compact(pressure_adjoint_kWave_RTdata, dx);

% Adjoint RT - kWave data
%pressure_adjoint_RT_kWaveData = matrix2cube(pressure_adjoint_RT_kWaveData_matrix, Nz);
%plot_projection_compact(pressure_adjoint_RT_kWaveData, dx);


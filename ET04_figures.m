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
%%  u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240_smooth.dat', ' ', 0);
%%  sensor_data_kWave = importdata(['./input_data/forwardSignal_kWave.dat'], ' ', 0);
%%  sensor_data_RT    = importdata(['./input_data/forwardSignal_RT.dat'], ' ', 0);

% Pressure
%%  pressure_adjoint_kWave_matrix         = importdata(['./output_data/pressure_adjoint_kWave.dat'], ' ', 0);
%%  pressure_adjoint_kWave_RT_data_matrix = importdata(['./output_data/pressure_adjoint_kWave_RT_data.dat'], ' ', 0);
%%  pressure_adjoint_RT_matrix            = importdata(['./output_data/pressure_adjoint_RT.dat'], ' ', 0);
%%  pressure_adjoint_RT_kWave_data_matrix = importdata(['./output_data/pressure_adjoint_RT_kWave_data.dat'], ' ', 0);

%========================================================================================================================
% INITIAL PRESSURE AND FORWARD SIGNAL
%========================================================================================================================
%%  % Initial Pressure
%%  u0 = matrix2cube(u0Matrix, Nz);
%%  plot_projection_compact(u0, dx);
%%  
%%  % kWave
%%  figure;
%%  imagesc(sensor_data_kWave(2:121, :));
%%  colorbar();
%%  
%%  % RT
%%  figure;
%%  imagesc(sensor_data_RT(2:121, :));
%%  colorbar();
%%  
%%  % RT
%%  figure;
%%  imagesc(sensor_data_kWave(2:121, :)-sensor_data_RT(2:121, :));
%%  colorbar();

%========================================================================================================================
% ADJOINT PRESSURE
%========================================================================================================================
%%  % Adjoint kWave
%%  pressure_adjoint_kWave = matrix2cube(pressure_adjoint_kWave_matrix, Nz);
%%  plot_projection_compact(pressure_adjoint_kWave, dx);

%%  % Adjoint RT
%%  pressure_adjoint_RT = matrix2cube(pressure_adjoint_RT_matrix, Nz);
%%  plot_projection_compact(pressure_adjoint_RT, dx);
%%  
%%  % Adjoint kWave - RT data
%%  pressure_adjoint_kWave_RTdata = matrix2cube(pressure_adjoint_kWave_RT_data_matrix, Nz);
%%  plot_projection_compact(pressure_adjoint_kWave_RTdata, dx);

%%  % Adjoint RT - kWave data
%%  pressure_adjoint_RT_kWaveData = matrix2cube(pressure_adjoint_RT_kWave_data_matrix, Nz);
%%  plot_projection_compact(pressure_adjoint_RT_kWaveData, dx);


%========================================================================================================================
% Forward
%========================================================================================================================
sensor_KW = importdata(['input_data/forwardSignal_kWave_100.dat'], ' ', 0);
figure;
imagesc(sensor_KW);
colorbar()

sensor_RT = importdata(['input_data/forwardSignal_RT_100.dat'], ' ', 0);
figure;
imagesc(sensor_RT);
colorbar()

figure;
error = sensor_RT-sensor_KW;
imagesc(error);
colorbar()

ERT = sum(sensor_RT(:).^2)
EKW = sum(sensor_KW(:).^2)
EER = sum(error(:).^2)
ratio = sqrt(ERT/EKW)


% Adjoint kWave
pressure_adjoint_kWave_matrix = importdata(['./output_data/pressure_adjoint_kWave_100.dat'], ' ', 0);
pressure_adjoint_kWave = matrix2cube(pressure_adjoint_kWave_matrix, Nz);
plot_projection_compact(pressure_adjoint_kWave, dx);

% Adjoint RT
pressure_adjoint_RT_matrix = importdata(['./output_data/pressure_adjoint_RT_100.dat'], ' ', 0);
pressure_adjoint_RT = matrix2cube(pressure_adjoint_RT_matrix, Nz);
plot_projection_compact(pressure_adjoint_RT, dx);


figure;
error = pressure_adjoint_RT-pressure_adjoint_kWave;
plot_projection_compact(error, dx);


EKW = sum(pressure_adjoint_kWave(:).^2)
ERT = sum(pressure_adjoint_RT(:).^2)
EER = sum(error(:).^2)
ratio = sqrt(EKW/ERT)

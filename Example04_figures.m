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
   

%========================================================================================================================
% PRIMAL AND DUAL DATA
%========================================================================================================================
% Load Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240_smooth.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
plot_projection_compact(u0, dx);


% kWave
sensor_data = h5read('./output_data/Example04_forward_output.h5', '/p');
figure;
imagesc(sensor_data(1:80, :));

% RT
time_signal = importdata(['./output_data/ForwardSignal.dat'], ' ', 0);
y0 = time_signal(2:end, :);
figure;
imagesc(y0(1:80, :));

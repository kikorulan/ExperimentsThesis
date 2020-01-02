% Heterogeneous Propagation Medium Example
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex16_real_recon3D_DS1;
cd /scratch0/NOT_BACKED_UP/frullan/ExperimentsThesis/Ex16_real_recon3D_DS1;

clear all;
%close all;

load ./settings/timeseriesdata_c1574_offset10-500;
load ./settings/Full_scan1_temp@850nm_t0[0]_dx[106µm]_dy[106µm]_dt[17ns]_39s33m16h_10-05-18_avg1_2D_raw(QRS);

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 80;           % number of grid points in the x (row) direction
Ny = 240;           % number of grid points in the y (column) direction
Nz = 240;           % number of grid points in the y (column) direction
dx = 5.3e-5;        % grid point spacing in the x direction [m]
dy = 5.3e-5;        % grid point spacing in the y direction [m]
dz = 5.3e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1574.2;
medium.sound_speed = c0;
medium.density = 1;

% compute time
dt = 1.6667e-8;
Nt = 490;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

%==============================
% SUBSAMPLE
%==============================
% REDUCE DIMENSIONS
% array_y = 1:144 -> array_y = 11:130
% array_z = 1:133 -> array_z =  3:123
y_min = 11;
y_max = 130;
z_min = 3;
z_max = 122;

% Full data - SELECT
sensor_data = sensor_data(y_min:y_max, z_min:z_max, :);

% Reshape sensor data
sensor_data_reshape = reshape(sensor_data, [Ny*Nz/2/2, Nt]);
forward_signal = [kgrid.t_array; sensor_data_reshape];
dlmwrite('./input_data/forwardSignal_reference_14400sensors_490timesteps.dat', forward_signal, 'delimiter', ' ');
%%  % Subsample data
%%  subsample_factor = 2;
%%  sensor_data_sub = sensor_data(1:subsample_factor:end, 1:subsample_factor:end, :);
%%  sensor_data_sub = reshape(sensor_data_sub, [Ny*Nz/(2*subsample_factor)/(2*subsample_factor), Nt]);
%%  forward_signal = [kgrid.t_array; sensor_data_sub];
%%  dlmwrite('./input_data/forwardSignal_reference_3600sensors_490timesteps.dat', forward_signal, 'delimiter', ' ');
% Pixel Pressure
pixelPressure = zeros(Nx*Nz, Ny);
dlmwrite('./input_data/pixelPressure_0.dat', pixelPressure, 'delimiter', ' ');
% Sound speed
sound_speed = c0*ones(Nx*Nz, Ny);
dlmwrite('./input_data/sound_speed.dat', sound_speed, 'delimiter', ' ');

% Subsample data
rng(1);
nSensors = 120;
nSensors_total = nSensors*nSensors;
subsample_factor = 4;
nSensors_sub = floor(nSensors_total/subsample_factor);
index_sensors = zeros(nSensors_sub, 3);
index_sensors(:, 1) = sort(datasample(1:nSensors_total, nSensors_sub, 'Replace', false));
sensor_data_sub = zeros(nSensors_sub, Nt);
sensor_mask = zeros(nSensors, nSensors);
for ii = 1:nSensors_sub
    index_z = ceil(index_sensors(ii, 1)/nSensors);
    index_y = index_sensors(ii, 1)-nSensors*(index_z-1);
    index_sensors(ii, 2:3) = [index_y index_z];
    sensor_data_sub(ii, :) = sensor_data(index_y, index_z, :);
    sensor_mask(index_y, index_z) = 1;
end
sensor_data = sensor_data_sub;
% Save sensor data
save input_data/sensor_data_3600.mat sensor_data;
% Save forward signal
forward_signal = [kgrid.t_array; sensor_data];
dlmwrite('input_data/forwardSignal_random_3600sensors_490timesteps.dat', forward_signal, 'delimiter', ' ');
% Save indices
dlmwrite('input_data/subsampleIndices_3600sensors.dat', index_sensors, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';

%==============================
% ADJOINT - GENERATED PARAMETERS
%==============================
% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Consider all sensors
source_adj.p_mask = zeros(Nx, Ny, Nz);
source_adj.p_mask(1, 1:2:end, 1:2:end) = sensor_mask;
source_adj.p = fliplr(sensor_data);
source_adj.p_mode = 'additive';
%source_adj.p_mode = 'dirichlet';
% Sensor
sensor_adj.mask = ones(Nx, Ny, Nz);
sensor_adj.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adj, sensor_adj, input_args{:}, 'SaveToDisk', 'input_data/ET16_adjoint_input_3600sensors.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i input_data/ET16_adjoint_input_3600sensors.h5 -o output_data/ET16_adjoint_output_3600sensors.h5 --p_final');

%=========================================================================
% PLOT
%=========================================================================
% ADJOINT - GENERATED PARAMETERS
PML_size = 10;
p0_adj_PML = h5read('./output_data/ET16_adjoint_output_3600sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
plot_projection(p0_adj, dx);

%=========================================================================
% PLOT SENSOR MASK
%=========================================================================
figure;
imagesc(sensor_mask);

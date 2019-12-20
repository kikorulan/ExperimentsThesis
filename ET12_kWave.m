% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex12_adjoint3D;

clear all;
close all;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
Nx = 256;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
Nz = 128;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
dz = 1e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
c0 = 1500;
c = c0*ones(Nx, Ny, Nz);
% Load sound speed
medium.sound_speed = c;
medium.density = 1;
c_matrix  = cube2matrix(c);
dlmwrite('sound_speed_1500.dat', c_matrix, 'delimiter', ' ');

% compute time
dt = 2e-8;
tMax = 2e-5;
kgrid.t_array = 0:dt:tMax;

% Load Initial Pressure
[Y, X, Z] = meshgrid(kgrid.y_vec, kgrid.x_vec, kgrid.z_vec);
S1 = (X-3*floor(Nx/10)*dx).^2 + Y.^2 + Z.^2;
S2 = (X-1*floor(Nx/10)*dx).^2 + Y.^2 + Z.^2;
S3 = (X+1*floor(Nx/10)*dx).^2 + Y.^2 + Z.^2;
sigma = floor(Nx/30)*dx;
u0_low  = exp(-S1/sigma/sigma);
u0_mid  = exp(-S2/sigma/sigma);
u0_high  = exp(-S3/sigma/sigma);

% Convert to cube
u0_low_matrix  = cube2matrix(u0_low);
u0_mid_matrix  = cube2matrix(u0_mid);
u0_high_matrix = cube2matrix(u0_high);
dlmwrite('u0_low.dat', u0_low_matrix, 'delimiter', ' ');
dlmwrite('u0_mid.dat', u0_mid_matrix, 'delimiter', ' ');
dlmwrite('u0_high.dat', u0_high_matrix, 'delimiter', ' ');
plot_projection(u0_low + u0_mid + u0_high, dx);

% Sources
source_low.p0 = u0_low;
source_mid.p0 = u0_mid;
source_high.p0 = u0_high;

%=========================================================================
% SIMULATION
%=========================================================================
% Set Library
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';

%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, floor(Ny/2), floor(Nz/2)) = 1;

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
PML_size = 40;
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Save to disk
kspaceFirstOrder3D(kgrid, medium, source_low, sensor, input_args{:}, 'SaveToDisk', 'Example12_forward_input_low.h5');
kspaceFirstOrder3D(kgrid, medium, source_mid, sensor, input_args{:}, 'SaveToDisk', 'Example12_forward_input_mid.h5');
kspaceFirstOrder3D(kgrid, medium, source_high, sensor, input_args{:}, 'SaveToDisk', 'Example12_forward_input_high.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_forward_input_low.h5 -o Example12_forward_output_low.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_forward_input_mid.h5 -o Example12_forward_output_mid.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_forward_input_high.h5 -o Example12_forward_output_high.h5');

%==============================
% Adjoint
%==============================
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};

% Read results
sensor_data_low = h5read('Example12_forward_output_low.h5', '/p');
sensor_data_mid = h5read('Example12_forward_output_mid.h5', '/p');
sensor_data_high = h5read('Example12_forward_output_high.h5', '/p');

% Consider all sensors
source_adjoint_low.p_mask = sensor.mask;
source_adjoint_mid.p_mask = sensor.mask;
source_adjoint_high.p_mask = sensor.mask;
source_adjoint_low.p = fliplr(sensor_data_low);
source_adjoint_mid.p = fliplr(sensor_data_mid);
source_adjoint_high.p = fliplr(sensor_data_high);

% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint_low, sensor_adjoint, input_args{:}, 'SaveToDisk', 'Example12_adjoint_input_low.h5');
kspaceFirstOrder3D(kgrid, medium, source_adjoint_mid, sensor_adjoint, input_args{:}, 'SaveToDisk', 'Example12_adjoint_input_mid.h5');
kspaceFirstOrder3D(kgrid, medium, source_adjoint_high, sensor_adjoint, input_args{:}, 'SaveToDisk', 'Example12_adjoint_input_high.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_adjoint_input_low.h5 -o Example12_adjoint_output_low.h5 --p_final');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_adjoint_input_mid.h5 -o Example12_adjoint_output_mid.h5 --p_final');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example12_adjoint_input_high.h5 -o Example12_adjoint_output_high.h5 --p_final');

%=========================================================================
% EXTRACT TO RAY TRACING FORMAT
%=========================================================================
% Read results - Forward
sensor_data_low = h5read('Example12_forward_output_low.h5', '/p');
sensor_data_mid = h5read('Example12_forward_output_mid.h5', '/p');
sensor_data_high = h5read('Example12_forward_output_high.h5', '/p');

save sensor_data_1600 sensor_data_low sensor_data_mid sensor_data_high kgrid;

% Read results - Adjoint
adjoint_pressure_low_PML = h5read('Example12_adjoint_output_low.h5', '/p_final');
adjoint_pressure_mid_PML = h5read('Example12_adjoint_output_mid.h5', '/p_final');
adjoint_pressure_high_PML = h5read('Example12_adjoint_output_high.h5', '/p_final');
adjoint_pressure_low  = adjoint_pressure_low_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
adjoint_pressure_mid  = adjoint_pressure_mid_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
adjoint_pressure_high = adjoint_pressure_high_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);

save adjoint_pressure_1600 adjoint_pressure_low adjoint_pressure_mid adjoint_pressure_high kgrid;


%==================================================
% INTERPOLATE
%==================================================
load sensor_data_1600;
tMax = 2e-5;
tForward = 0:5e-9:tMax;

% DT1
spline_low = spline(kgrid.t_array, sensor_data_low, tForward);
spline_mid  = spline(kgrid.t_array, sensor_data_mid, tForward);
spline_high = spline(kgrid.t_array, sensor_data_high, tForward);

% Forward
sensor_data_low_time  = [tForward; spline_low];
sensor_data_mid_time  = [tForward; spline_mid];
sensor_data_high_time = [tForward; spline_high];
dlmwrite('forwardSignal_kWave_1600_low_5e-9.dat', sensor_data_low_time, 'delimiter', ' ');
dlmwrite('forwardSignal_kWave_1600_mid_5e-9.dat', sensor_data_mid_time, 'delimiter', ' ');
dlmwrite('forwardSignal_kWave_1600_high_5e-9.dat', sensor_data_high_time, 'delimiter', ' ');

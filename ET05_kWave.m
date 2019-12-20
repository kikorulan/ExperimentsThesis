% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex05_simulation3D_het;

clear all;
close all;

runForward = 1;
runAdjoint = 1;
runAdjoint_RTData = 0;
saveData = 0;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
Nx = 80;           % number of grid points in the x (row) direction
Ny = 240;           % number of grid points in the y (column) direction
Nz = 240;           % number of grid points in the y (column) direction
dx = 5.3e-5;        % grid point spacing in the x direction [m]
dy = 5.3e-5;        % grid point spacing in the y direction [m]
dz = 5.3e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%===============================
% BUILD SOUND SPEED
%===============================
rng(1);
gridR = gridRT(Nx, 1, Ny, 1, Nz, 1);
% Set sound speed
c0 = 1580.0;
c = c0*ones(Nx, Ny, Nz);
gridR.setCMatrix(c);
inc = 0.06;
gridR.randomiseC(.15, inc);

% Plot sound speed
plot_cube(gridR.c);
% Save data
c0_matrix = cube2matrix(gridR.c);
dlmwrite('./input_data/sound_speed.dat', c0_matrix, 'delimiter', ' ');


%==============================
% define the properties of the propagation medium
%==============================
cMatrix = importdata('./input_data/sound_speed.dat', ' ', 0);
c = matrix2cube(cMatrix, Nz);
% Load sound speed
medium.sound_speed = c;
medium.density = 1;

% compute time
dt = 1.5e-8;
Nt = 539;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

% Load Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, u0, true);
source.p0 = max(0, source.p0);
initial_pressure_veins_smooth = source.p0;
u0_smooth = cube2matrix(initial_pressure_veins_smooth);
dlmwrite('./input_data/initial_pressure_veins_80x240x240_smooth.dat', u0_smooth, 'delimiter', ' ');
save input_data/initial_pressure_veins_smooth initial_pressure_veins_smooth;
plot_projection(initial_pressure_veins_smooth, dx);

%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
if(runForward)
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, 1:2:end, 1:2:end)   = 1;
%sensor.mask(end, :, :) = 1;
%sensor.mask(:, 1, :)   = 1;
%sensor.mask(:, end, :) = 1;
%sensor.mask(:, :, 1)   = 1;
%sensor.mask(:, :, end) = 1;

% Number of sensors
numberSensors = sum(sensor.mask(:))

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
save input_data/sensor_data_veins.mat kgrid medium source sensor input_args;

% Save to disk
filename = 'input_data/Example05_forward_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i ./input_data/Example05_forward_input.h5 -o ./output_data/Example05_forward_output.h5');
end
%==============================
% Adjoint
%==============================
if(runAdjoint)
% Read results
sensor_data = h5read('output_data/Example05_forward_output.h5', '/p');
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);

% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example05_adjoint_input.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i ./input_data/Example05_adjoint_input.h5 -o ./output_data/Example05_adjoint_output.h5 --p_final');
end

%==============================
% Adjoint - RT DATA
%==============================
if(runAdjoint_RTData)
% Read results
sensor_data_RT = importdata(['./input_data/forwardSignal_RT.dat'], ' ', 0);
sensor_data_RT = sensor_data_RT(2:end, :);
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);

% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data_RT);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example05_adjoint_RTdata_input.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i input_data/Example05_adjoint_RTdata_input.h5 -o output_data/Example05_adjoint_RTdata_output.h5 --p_final');
end

%=========================================================================
% EXTRACT TO RAY TRACING FORMAT
%=========================================================================
if(saveData)
% kWave
sensor_data = h5read('output_data/Example05_forward_output.h5', '/p');
rt_data = [kgrid.t_array; sensor_data];
dlmwrite('input_data/forwardSignal_kWave.dat', rt_data, 'delimiter', ' ');


% Adjoint 
pressure_adjoint_PML = h5read('./output_data/Example05_adjoint_output.h5', '/p_final');
PML_size = 10;
pressure_adjoint = pressure_adjoint_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
pressure_adjoint_cube = cube2matrix(pressure_adjoint);
dlmwrite('output_data/pressure_adjoint_kWave.dat', pressure_adjoint_cube, 'delimiter', ' ');

% Adjoint - RT data
pressure_adjoint_PML = h5read('./output_data/Example05_adjoint_RTdata_output.h5', '/p_final');
PML_size = 10;
pressure_adjoint = pressure_adjoint_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
pressure_adjoint_cube = cube2matrix(pressure_adjoint);
dlmwrite('output_data/pressure_adjoint_kWave_RT_data.dat', pressure_adjoint_cube, 'delimiter', ' ');
end

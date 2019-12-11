% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex08_adjoint2D;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 512;           % number of grid points in the x (row) direction
Ny = 1024;           % number of grid points in the y (column) direction
dx = 2.5e-5;        % grid point spacing in the x direction [m]
dy = 2.5e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
medium.sound_speed = c0*ones(Nx, Ny);
medium.density = 1;
    
% compute time
dt = 2e-9;
tMax = 2e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Axis
[Y, X] = meshgrid(kgrid.y_vec, kgrid.x_vec);

% Build domain
v1 = -0.3;
Z1 = X.^2 + (Y+floor(Ny/10)*dy).^2;
Z2 = X.^2 + (Y-floor(Ny/10)*dy).^2;
Z3 = X.^2 + (Y-floor(3*Ny/10)*dy).^2;
sigma = floor(Nx/32)*dx;
u0_low  = exp(-Z1/sigma/sigma);
u0_mid  = exp(-Z2/sigma/sigma);
u0_high = exp(-Z3/sigma/sigma);

% smooth the initial pressure distribution and restore the magnitude
source_low.p0 = u0_low;
source_mid.p0 = u0_mid;
source_high.p0 = u0_high;

%=========================================================================
% SIMULATION - FORWARD
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(floor(Nx/2), 1) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data_low = kspaceFirstOrder2D(kgrid, medium, source_low, sensor, input_args{:});
sensor_data_mid = kspaceFirstOrder2D(kgrid, medium, source_mid, sensor, input_args{:});
sensor_data_high = kspaceFirstOrder2D(kgrid, medium, source_high, sensor, input_args{:});

save sensor_data.mat kgrid sensor source_low source_mid source_high medium c0 dt sensor_data_low sensor_data_mid sensor_data_high input_args;

%=========================================================================
% SIMULATION - ADJOINT
%=========================================================================
% Input arguments
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% Sensors
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
% Source
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(floor(Nx/2), 1) = 1;

% Simulations
source_adjoint.p = fliplr(sensor_data_low);
adjoint_pressure_low = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});
source_adjoint.p = fliplr(sensor_data_mid);
adjoint_pressure_mid = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});
source_adjoint.p = fliplr(sensor_data_high);
adjoint_pressure_high = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

save adjoint_data.mat kgrid adjoint_pressure_low adjoint_pressure_mid adjoint_pressure_high input_args;

%=========================================================================
% VISUALISATION
%=========================================================================
figure;
plot(kgrid.t_array, sensor_data_low, 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid, 'Color', 'g');
plot(kgrid.t_array, sensor_data_high, 'Color', 'b');
legend('IP low', 'IP mid', 'IP high');
xlabel('t (s)');
ylabel('amplitude');


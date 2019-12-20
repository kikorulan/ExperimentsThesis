% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex11_forward3D_het;

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
medium.sound_speed = gridR.c;
medium.density = 1;
dlmwrite('sound_speed.dat', c0_matrix, 'delimiter', ' ');

% compute time
dt = 1e-8;
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
input_args = {'PMLSize', 40, 'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Save to disk
kspaceFirstOrder3D(kgrid, medium, source_low, sensor, input_args{:}, 'SaveToDisk', 'Example11_forward_input_low.h5');
kspaceFirstOrder3D(kgrid, medium, source_mid, sensor, input_args{:}, 'SaveToDisk', 'Example11_forward_input_mid.h5');
kspaceFirstOrder3D(kgrid, medium, source_high, sensor, input_args{:}, 'SaveToDisk', 'Example11_forward_input_high.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example11_forward_input_low.h5 -o Example11_forward_output_low.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example11_forward_input_mid.h5 -o Example11_forward_output_mid.h5');
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i Example11_forward_input_high.h5 -o Example11_forward_output_high.h5');

%=========================================================================
% EXTRACT TO RAY TRACING FORMAT
%=========================================================================
% Read results
sensor_data_low = h5read('Example11_forward_output_low.h5', '/p');
sensor_data_mid = h5read('Example11_forward_output_mid.h5', '/p');
sensor_data_high = h5read('Example11_forward_output_high.h5', '/p');

save sensor_data sensor_data_low sensor_data_mid sensor_data_high kgrid;

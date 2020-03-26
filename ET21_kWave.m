%================================================================================
% ET21
% k-Wave simulation for single GB
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex21_singleGB;
close all;
clear all;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
N  = 128;           % number of grid points in the x (row) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
kgrid = makeGrid(N, dx, N, dx);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
medium.sound_speed = c0*ones(N, N);
medium.density = 1;
    
% compute time
dt = 2e-8;
tMax = 1e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

%==============================
% Build initial pressure
%==============================
% Gaussian beam Phantom
load singleGB;
sourceKW.p0 = real(singleGB);
% smooth the initial pressure distribution and restore the magnitude
%sourceKW.p0 = smooth(kgrid, sourceKW.p0, true);

%=========================================================================
% SIMULATION - LINEAR ARRAY SENSOR
%=========================================================================
%%%% Forward
% Define the sensors
sensor.mask = zeros(N, N);
sensor.mask(:, 1) = 1;
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, sourceKW, sensor, input_args{:});
save sensor_data_singleGB_kWave.mat kgrid sensor sourceKW medium sensor_data input_args;
%%%% Adjoint
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(:, 1) = 1;
source_adjoint.p = fliplr(sensor_data);
adjoint_singleGB_kWave = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});
save adjoint_singleGB_kWave adjoint_singleGB_kWave;

%=========================================================================
% VISUALISATION
%=========================================================================
X = 0:dx:(N-1)*dx;
Y = X;

%==============================
% Pressure Time Series
%==============================
load sensor_data_singleGB_kWave.mat;
figure;
surf(kgrid.t_array, X, real(sensor_data), 'EdgeColor', 'none');
view(2);
colorbar();
axis tight;

%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex08_adjoint2D;

close all;
clear all;

load sensor_data.mat;

run colourMap;
saveData = 0;
%========================================
% Grid definition
%========================================
Nx = 512;           % number of grid points in the x (row) direction
Ny = 1024;           % number of grid points in the y (column) direction
dx = 2.5e-5;        % grid point spacing in the x direction [m]
dy = 2.5e-5;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Make matrix symmetric
factorResize = 1;
c = imresizen(medium.sound_speed, factorResize);
gridR.setCMatrix(c);

%========================================
% Impulse Response
%========================================
% Set time
dt = 8e-9;
%dt = min(gridR.dx, gridR.dy)/c0/2;
tMax = 2e-5;
gridR.setTime(dt, tMax);
% Compute impulse response
gridR.impulse_additive('load');
save gridRT_impulse.mat gridR c0 dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;

gridR.setTime(dt, tMax);
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 8000;% 800
nSources = 1;%256

% Parametrisation
tStep = dt

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-2)/2*gridR.dx, 0);

% Sources
source(1) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);
source(2) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);
source(3) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);
%source(4) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);

% Set initial pressure
gridR.setUMatrix(imresizen(source_low.p0, factorResize));
gridR.computeHamil(source(1), 'p');
gridR.setUMatrix(imresizen(source_mid.p0, factorResize));
gridR.computeHamil(source(2), 'p');
gridR.setUMatrix(imresizen(source_high.p0, factorResize));
gridR.computeHamil(source(3), 'p');
%gridR.setUMatrix(imresizen(source_low.p0 + source_mid.p0 + source_high.p0, factorResize));
%gridR.computeHamil(source(4), 'p');

%========================================
% Compute reverse signal
%========================================
gridR.inverse_filter(100);
sensor_low = spline(kgrid.t_array, sensor_data_low(1, :), 0:dt:tMax);
sensor_mid = spline(kgrid.t_array, sensor_data_mid(1, :), 0:dt:tMax);
sensor_high = spline(kgrid.t_array, sensor_data_high(1, :), 0:dt:tMax);

source(1).setForwardSignal(sensor_low);
source(2).setForwardSignal(sensor_mid);
source(3).setForwardSignal(sensor_high);
%source(4).setForwardSignal(sensor_low + sensor_mid + sensor_high);
for n = 1:3
    gridR.inverse_beam(source(n));
end
%Rgrid.computeAdjointParallel(source);

adjoint_pressure_low_RT  = source(1).pixelAReverse;
adjoint_pressure_mid_RT  = source(2).pixelAReverse;
adjoint_pressure_high_RT = source(3).pixelAReverse;
%adjoint_pressure_all_RT  = source(4).pixelAReverse;

if(saveData)
%save adjoint_pressure_RT_2e-9 adjoint_pressure_low_RT adjoint_pressure_mid_RT adjoint_pressure_high_RT adjoint_pressure_all_RT gridR factorResize;
save adjoint_pressure_RT_2e-9 adjoint_pressure_low_RT adjoint_pressure_mid_RT adjoint_pressure_high_RT gridR factorResize;
end
%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);


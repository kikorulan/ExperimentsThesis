%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex03_simulation2D_het;

close all;
%clear all;

load sensor_data_kWave.mat;
load initial_pressure;
%========================================
% Rgrid definition
%========================================
Nx = 128;           % number of Rgrid points in the x (row) direction
Ny = 256;           % number of Rgrid points in the y (column) direction
dx = 2e-4;        % Rgrid point spacing in the x direction [m]
dy = 2e-4;        % Rgrid point spacing in the y direction [m]
Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
Rgrid.setCMatrix(medium.sound_speed);
Rgrid.setUMatrix(initial_pressure);

%========================================
% Impulse Response
%========================================
% Set time
dt = 5e-8;
tMax = 4e-5;
Rgrid.setTime(dt, tMax);
% Compute impulse response
Rgrid.impulse_additive('IV');

save gridRT_impulse.mat Rgrid;

%========================================
% Ray Shooting Parameters
%========================================
load gridRT_impulse;
cMax = max(Rgrid.c(:));
dt = 5e-8;
%dt = min(Rgrid.dx, Rgrid.dy)/cMax/2;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 2000;% 800
nSources = 764;%256

% Parametrisation
tMax = sqrt((Rgrid.Nx*Rgrid.dx)^2+(Rgrid.Ny*Rgrid.dy)^2)/cMax;
tStep = dt;

%========================================
% Forward Problem
%========================================
% Sources locations
clear x;
for n = 1:Rgrid.Nx
    x{n} = cat(3, (n-1)*Rgrid.dx, 0);
end
for n = 1:Rgrid.Ny-2
    x{2*n-1+Rgrid.Nx}   = cat(3,                   0, n*Rgrid.dy);
    x{2*n  +Rgrid.Nx}   = cat(3, (Rgrid.Nx-1)*Rgrid.dx, n*Rgrid.dy);
end
for n = 1:Rgrid.Nx
    x{n+Rgrid.Nx+2*(Rgrid.Ny-2)} = cat(3, (n-1)*Rgrid.dx, (Rgrid.Ny-1)*Rgrid.dy);
end

%%  for n = 1:nSources
%%      disp(['Source ', int2str(n)]);
%%      source = Rgrid.newSource(x{n}, 0, 2*pi-0.01, nRays, tStep, tMax);
%%      Rgrid.computeHamil(source, 'p');
%%      signalRT(n, :) = source.aForward;
%%  end

source = Rgrid.computeForwardParallel(x, 0, 2*pi-0.01, nRays, tStep, tMax, 'p', true);
signalRT = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    signalRT(n, :) = source(n).aForward;
end

save signalRT signalRT;

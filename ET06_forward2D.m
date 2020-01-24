%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex06_forward2D;

close all;
%clear all;

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
nRays = 4000;% 800
nSources = 1;%256

% Parametrisation
tStep = dt;

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

% Sources
%source(1) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);
source(1) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(2) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(3) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);

% Set initial pressure
gridR.setUMatrix(imresizen(source_low.p0, factorResize));
gridR.computeHamil(source(1), 'p');
gridR.setUMatrix(imresizen(source_mid.p0, factorResize));
gridR.computeHamil(source(2), 'p');
gridR.setUMatrix(imresizen(source_high.p0, factorResize));
gridR.computeHamil(source(3), 'p');

%%  clear source;
%%  source = gridR.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

%==============================
% Save results
%==============================
normRT = max(real(source(1).aForward));
sensor_RT_low = source(1).aForward;
sensor_RT_mid = source(2).aForward;
sensor_RT_high = source(3).aForward;

if(saveData)
save sensor_data_RT_8e-9.mat factorResize gridR normRT sensor_RT_low sensor_RT_mid sensor_RT_high;
end

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);


%==============================
% Plot ray trajectories
%==============================
position = [700 500 320 600];
source(1).plot_rays(gridR, 100, 'mm');
pbaspect([1 2 1]);
set(gcf, 'pos', position);
title('')
set(gca,'FontSize',13);
saveas(gcf, 'Example06_rayTrajectories', 'epsc');
saveas(gcf, 'Example06_rayTrajectories.fig');





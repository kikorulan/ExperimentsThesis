%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex01_caustics;

close all;
%clear all;

% Choose which caustics to run
caustic1 = 1;
caustic2 = 1;

%========================================
% Grid definition
%========================================
Nx = 256;           % number of grid points in the x (row) direction
Ny = 768;           % number of grid points in the y (column) direction
dx = 5e-5;        % grid point spacing in the x direction [m]
dy = 5e-5;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

%========================================
% Impulse Response
%========================================
% c = 1500;
% gridR.setCMatrix(c*ones(Nx, Ny));
% % Set time
% dt = 4e-8;
% %dt = min(gridR.dx, gridR.dy)/c0/2;
% tMax = 2e-5;
% gridR.setTime(dt, tMax);
% % Compute impulse response
% gridR.impulse_additive('IV');
% save gridRT_impulse.mat gridR c0 dt tMax;


% Load data
load gridRT_impulse;

%========================================
% Global parameters
%========================================
dt = 1e-8;
tMax = 6e-5;
gridR.setTime(dt, tMax);

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

position = [700 500 320 900];
positionYBar = [700 700 375 900];
positionBar = [700 700 375 900];
set(0,'DefaultFigurePaperPositionMode','auto');
fontSize = 16;
%==========================================================================================
%========                      ============================================================
%========    CAUSTIC TYPE 1    ============================================================
%========                      ============================================================
%==========================================================================================
if(caustic1)

%==========================================
% SOUND SPEED
%==========================================
c0 = 1500;
v1 = -0.3;
kernelSize = 20;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
vx = kernelSize + floor(Nx/2);
vy = floor(Ny/4);
focus = floor(16*Ny/60);
thr = 10;
c = addParabola(c, vx, vy, focus, c0*v1, thr);
%c(:, Ny/10:end) = c(:, Ny/10:end)*(1 + v1);
% Kernel convolution
K = ones(kernelSize);
%cConv = c;
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
gridR.setCMatrix(medium.sound_speed);

%==========================================
% RAY SHOOTING
%==========================================
% Measure computational time
tic;
start_time = clock;

% Sources
nRays = 200;
source(1) = gridR.newSource(x{1}, pi/2-pi/4, pi/2+pi/4, nRays, dt, tMax);
% Set initial pressure
gridR.computeHamil(source(1), 'p');
  
%==========================================
% Plot results
%==========================================

% Sound speed
gridR.plot_soundSpeed();
set(gcf, 'pos', positionYBar);
colorbar();
box on;
set(gca, 'FontSize', fontSize);
saveas(gcf, 'Example01_2D_caustic1_C', 'png'); 
saveas(gcf, 'Example01_2D_caustic1_C.fig'); 

% Ray trajectories
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
nColours = nRays;
colourList = summer(nColours);
[nRays nSteps dim] = size(source(1).x);
for j = 1:1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(1e3*source(1).x(j, :, 1), 1e3*source(1).x(j, :, 2), 'Color', colourList(colourNum, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
% Red mark
%xR = 6.3;
%yR = 15.6;
%plot(xR, yR, 'LineWidth', 3, 'Color', 'r', 'Marker', 'o', 'MarkerSize', 12);
x1 = 6.3:0.1:13;
y1 = 38 - 0.16*(x1-18).^2;
x2 = 0:0.1:6.3;
y2 = 38 - 0.16*(x2+(18-6.3-6.3)).^2;
plot(x1, y1, 'LineWidth', 2, 'Color', 'r');
plot(x2, y2, 'LineWidth', 2, 'Color', 'r');
% Save
set(gca, 'FontSize', fontSize);
saveas(gcf, 'Example01_2D_caustic1_rays', 'epsc');
saveas(gcf, 'Example01_2D_caustic1_rays.fig');


% Measure time
end_time = clock;
disp(['  total computation time for Caustic 1: ' num2str(etime(end_time, start_time))]);
end

%==========================================================================================
%========                      ============================================================
%========    CAUSTIC TYPE 2    ============================================================
%========                      ============================================================
%==========================================================================================
if(caustic2)
%==========================================
% SOUND SPEED
%==========================================
c0 = 1500;
v1 = 0.3;
kernelSize = 20;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
slope = 10;
x0 = floor(Nx/2);
y0 = floor(Ny/2);
c = addSlope(c, slope, x0, y0, c0*v1);
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
gridR.setCMatrix(medium.sound_speed);

%==========================================
% RAY SHOOTING
%==========================================
% Measure computational time
tic;
start_time = clock;
nRays = 200;
% Sources
%source(1) = gridR.newSource(x{1}, pi/2, pi/2+pi/10, nRays, dt, tMax);
source(1) = gridR.newSource(x{1}, 0, pi, nRays, dt, tMax);
% Set initial pressure
gridR.computeHamil(source(1), 'p');
  
%==========================================
% Plot results
%==========================================

%==============================
% Sound speed
%==============================
gridR.plot_soundSpeed();
set(gcf, 'pos', positionYBar);
colorbar();
set(gca, 'FontSize', fontSize);
saveas(gcf, 'Example01_2D_caustic2_C', 'png'); 
saveas(gcf, 'Example01_2D_caustic2_C.fig'); 

%==============================
% Ray trajectories
%==============================
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
nColours = nRays;
colourList = summer(nColours);
[nRays nSteps dim] = size(source(1).x);
for j = 1:1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(1e3*source(1).x(j, :, 1), 1e3*source(1).x(j, :, 2), 'Color', colourList(colourNum, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
% Red mark
xR = 4.6:1e-1:9;
yR = 9*(xR-5.95) + 20;
plot(xR, yR, 'LineWidth', 2, 'Color', 'r');

% Save
set(gca, 'FontSize', fontSize);
saveas(gcf, 'Example01_2D_caustic2_rays', 'epsc'); 
saveas(gcf, 'Example01_2D_caustic2_rays.fig'); 


% Measure time
end_time = clock;
disp(['  total computation time for Caustic 2: ' num2str(etime(end_time, start_time))]);
end

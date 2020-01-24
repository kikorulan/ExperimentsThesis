%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex02_simulation2D_homo;

close all;
%clear all;

load sensor_data_kWave.mat;

%========================================
% Rgrid definition
%========================================
Nx = 128;           % number of Rgrid points in the x (row) direction
Ny = 256;           % number of Rgrid points in the y (column) direction
dx = 1e-4;        % Rgrid point spacing in the x direction [m]
dy = 1e-4;        % Rgrid point spacing in the y direction [m]
Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
Rgrid.setCMatrix(medium.sound_speed);

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
% Sensor Selection
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

% Sources locations
n1 = round(Rgrid.Nx + 2*(round(Rgrid.Ny/3) - 2)   + 2);     % 1st sensor: Nx + 2*(Ny/3-2)   + 2
n2 = round(Rgrid.Nx + 2*(round(2*Rgrid.Ny/3) - 2) + 1);     % 2nd sensor: Nx + 2*(2*Ny/3-2) + 1
n3 = round(Rgrid.Nx + 2*(Rgrid.Ny - 2) + round(Rgrid.Nx/2)); % 3rd sensor: Nx + 2*(Ny-2)     + Nx/2
sensor1 = x{n1};
sensor2 = x{n2};
sensor3 = x{n3};
sourceSel(1) = Rgrid.newSource(sensor1, pi/2, 3*pi/2, nRays, tStep, tMax);
sourceSel(2) = Rgrid.newSource(sensor2, -pi/2, pi/2, nRays, tStep, tMax);
sourceSel(3) = Rgrid.newSource(sensor3, pi, 2*pi, nRays, tStep, tMax);
Rgrid.computeHamil(sourceSel(1), 'p');
Rgrid.computeHamil(sourceSel(2), 'p');
Rgrid.computeHamil(sourceSel(3), 'p');

%==================================================================================
% Plot results
%==================================================================================

position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];

fontSize = 13;

%==============================
% SENSORS
%==============================
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
axisGrid = [0 1e3*X(end) 0 1e3*Y(end)];
nRaysPlot = 140;
colorList = cool(nRaysPlot);
% Sensor 1
figure;
hold on;
axis(axisGrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(1).x(index, :, 1), 1e3*sourceSel(1).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor1(1), 1e3*sensor1(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','r');
box on;
pbaspect([1 2 1]);
set(gcf, 'pos', position);
xlabel('x [mm]');
ylabel('y [mm]');
set(gca,'FontSize',fontSize);
saveas(gcf, 'Example02_Sensor1_rays', 'epsc');
saveas(gcf, 'Example02_Sensor1_rays.fig'); 

% Sensor 2
figure;
hold on;
axis(axisGrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(2).x(index, :, 1), 1e3*sourceSel(2).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor2(1), 1e3*sensor2(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','g');
box on;
pbaspect([1 2 1]);
set(gcf, 'pos', position);
xlabel('x [mm]');
set(gca,'FontSize',fontSize)
saveas(gcf, 'Example02_Sensor2_rays', 'epsc');
saveas(gcf, 'Example02_Sensor2_rays.fig'); 


% Sensor 3
figure;
hold on;
axis(axisGrid);
for n = 1:nRaysPlot
    index = n*floor(nRays/nRaysPlot);
    plot(1e3*sourceSel(3).x(index, :, 1), 1e3*sourceSel(3).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
end
plot(1e3*sensor3(1), 1e3*sensor3(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','b');
box on;
pbaspect([1 2 1]);
set(gcf, 'pos', position);
xlabel('x [mm]');
set(gca,'FontSize',fontSize)
saveas(gcf, 'Example02_Sensor3_rays', 'epsc');
saveas(gcf, 'Example02_Sensor3_rays.fig'); 




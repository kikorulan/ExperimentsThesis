%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex18_caustics;

%clear all;
close all;

load sensor_data.mat;

runSim = 0;
saveFigures = 1;
%========================================
% Grid definition
%========================================
if(runSim)
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Make matrix symmetric
c = medium.sound_speed;
gridR.setCMatrix(c);


%========================================
% Impulse Response
%========================================
% Set time
dt = 1e-8;
tMax = 2e-5;
gridR.setTime(dt, tMax);
% Compute impulse response
gridR.impulse_additive('load');
save gridRT_impulse.mat gridR c0 dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 2000;% 800
nSources = 1;%256

% Parametrisation
tStep = dt;

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

% Sources
source(1) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(2) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(3) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
% Set initial pressure
gridR.setUMatrix(source_low.p0);
gridR.computeHamil(source(1), 'p');
gridR.setUMatrix(source_mid.p0);
gridR.computeHamil(source(2), 'p');
gridR.setUMatrix(source_high.p0);
gridR.computeHamil(source(3), 'p');
end
 
%==================================================================================
% Plot results
%==================================================================================

position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar  = [700 700 375 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Initial Pressure
%==============================
pressure = source_low.p0 + source_mid.p0 + source_high.p0;
figure;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
hold on;
surf(1e3*gridR.xAxis, 1e3*gridR.yAxis, pressure', 'EdgeColor', 'none');
view(2);
xlabel('x [mm]');
ylabel('');
set(gca,'FontSize',13);
%set(gca, 'ytick', []);
%set(gca, 'yticklabel', []);
colorbar();
set(gcf, 'pos', positionYBar);
if(saveFigures)
saveas(gcf, 'Example18_U', 'png'); 
saveas(gcf, 'Example18_U.fig'); 
end

%==============================
% Sound speed
%==============================
gridR.plot_soundSpeed();
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
%title('Initial Pressure');
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'ytick', []);
%set(gca, 'yticklabel', []);
colorbar();
box on;
set(gca,'FontSize',13);
set(gcf, 'pos', positionYBar);
pbaspect([1 2 1]);
if(saveFigures)
saveas(gcf, 'Example18_SoundSpeed', 'png');
saveas(gcf, 'Example18_SoundSpeed.fig');
end

%==============================
% Isotime curve
%==============================
%%  nCurves = 20;
%%  xIso = gridR.findIsoTime(source(1), 20, 1e-4, 1.2e-4);
%%  % Figure
%%  figure;
%%  flipGrid = 1./gridR.c';
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  hold on;
%%  surf(gridR.xAxis, gridR.yAxis, flipGrid, 'EdgeColor', 'none');
%%  hold on;
%%  for n = 1:nCurves
%%      plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(1501, [1 size(xIso, 2)]), '-m');
%%  end
%%  for j = 1:20:nRays
%%      plot3(source(2).x(j, :, 1), source(2).x(j, :, 2), repmat(1501, [1 source(2).nPoints]), '-g');
%%  end
%%  view(2);
%%  legend('Isocurves for times tMin = 1e-4s, tMax = 1.2e-4s');

%==============================
% Sound speed + rays
%==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = nRays;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nRays
%%      colourNum = floor(nColours*(j-1)/nRays) + 1;
%%      plot3(source(1).x(j, :, 1), source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  %surf(gridR.xAxis, gridR.yAxis, 1./gridR.c', 'EdgeColor', 'none');
%%  imagesc(gridR.xAxis, gridR.yAxis, 1./gridR.c');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_rays_soundSpeed.fig');

%==============================
% Ray trajectories
%==============================
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
nColours = nRays;
colourList = cool(nColours);
[nRays nSteps dim] = size(source(1).x);
zVec = 1600*ones(1, nSteps);
for j = 1:20:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot3(1e3*source(1).x(j, :, 1), 1e3*source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
end
xlabel('x [mm]');
%ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
set(gca,'FontSize',13);
if(saveFigures)
saveas(gcf, 'Example18_rays', 'png');
saveas(gcf, 'Example18_rays.fig');
end

%==============================
% Amplitude
%==============================
nSR = 500;
factor_sub = 40;
factor_time = 1e6;
colours = winter(nSR);

% Q
figure;
subplot(2, 1, 1);
for n = 1:factor_sub:nSR
    plot(factor_time*gridR.tForward, real(source(1).q(n+floor((nRays-nSR*2)/2), :)), 'Color', colours(n, :));
    hold on;
end;
grid on;
box on;
xlabel('t [\mus]');
ylabel('q');
pbaspect([2.5 1 1]);
set(gca, 'FontSize', 13);
% Amplitude
subplot(2, 1, 2);
for n = 1:factor_sub:nSR
    semilogy(factor_time*gridR.tForward, real(source(1).amplitude(n+floor((nRays-nSR*2)/2), :)), 'Color', colours(n, :));
    hold on;
end;
grid on;
box on;
xlabel('t [\mus]');
ylabel('Amplitude');
pbaspect([2.5 1 1])
set(gca, 'FontSize', 13);
set(gcf, 'pos', [700 700 600 600]);
if(saveFigures)
saveas(gcf, 'Example18_QA', 'epsc');
saveas(gcf, 'Example18_QA.fig');
end

% Trajectories
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
colourList = winter(nSR);
[nRays nSteps dim] = size(source(1).x);
for n = 1:factor_sub:nSR
    plot(1e3*source(1).x(n+floor((nRays-nSR*2)/2), :, 1), 1e3*source(1).x(n+floor((nRays-nSR*2)/2), :, 2), 'Color', colourList(n, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
set(gca,'FontSize',13);
pbaspect([1 2 1])
if(saveFigures)
saveas(gcf, 'Example18_subrays', 'epsc');
saveas(gcf, 'Example18_subrays.fig');
end


%==============================
% Time Signals
%==============================
% Norms
factor = 1e6;
normKWave = max(sensor_data_low(1, :));
normRT    = normKWave;
%normRT    = max(real(source(1).aForward));
%normRT_2 = max(real(sourceRT(1).aForward));
figure;
hold on;
axis([0 20 -1.5 2]);
grid on;
box on;
plot(factor*gridR.tForward, source(1).aForward/normRT, 'Color', 'r', 'LineWidth', 2);
plot(factor*kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(factor*gridR.tForward, source(2).aForward/normRT, 'Color', 'g', 'LineWidth', 2);
plot(factor*kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(factor*gridR.tForward, source(3).aForward/normRT, 'Color', 'b', 'LineWidth', 2);
plot(factor*kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('HG       - bottom', 'kWave - bottom', 'HG       - middle', 'kWave - middle', 'HG       - top', 'kWave - top');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gcf, 'pos', [700 700 600 500]);
set(gca, 'FontSize', 13);
if(saveFigures)
saveas(gcf, 'Example18_signals', 'epsc');
saveas(gcf, 'Example18_signals.fig');
end

figure;
hold on;
axis([0 20 -1 1]);
grid on;
box on;
plot(factor*gridR.tForward, source(1).aForward/normRT - sensor_data_low(1, :)/normKWave, 'Color', 'r', 'LineWidth', 2);
plot(factor*gridR.tForward, source(2).aForward/normRT - sensor_data_mid(1, :)/normKWave, 'Color', 'g', 'LineWidth', 2);
plot(factor*gridR.tForward, source(3).aForward/normRT - sensor_data_high(1, :)/normKWave, 'Color', 'b', 'LineWidth', 2);
legend('Bottom', 'Middle', 'Top');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gcf, 'pos', [700 700 600 500]);
set(gca, 'FontSize', 13);
if(saveFigures)
saveas(gcf, 'Example18_error', 'epsc');
saveas(gcf, 'Example18_error.fig');
end

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

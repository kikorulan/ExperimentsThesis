%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex06_basic2D;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Grid definition
%========================================
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
%%  % Set time
%%  dt = 2e-8;
%%  %dt = min(gridR.dx, gridR.dy)/c0/2;
%%  tMax = 2e-5;
%%  gridR.setTime(dt, tMax);
%%  % Compute impulse response
%%  gridR.impulse_additive('IV');
%%  save gridRT_impulse.mat gridR c0 dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800
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

%%  clear source;
%%  source = gridR.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

% Save results
%save gridRT.mat gridR nRays nSources x -v7.3;
  
%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar = [700 700 375 600];
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
%title('Initial Pressure');
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'ytick', []);
%set(gca, 'yticklabel', []);
colorbar();
set(gcf, 'pos', positionYBar);
%saveas(gcf, 'Example46_U', 'png'); 
%saveas(gcf, 'Example46_U.fig'); 

%==============================
% Amplitude
%==============================
%%  % Imagesc
%%  figure;
%%  surf(real(log(real(source(1).amplitude))), 'EdgeColor', 'none');
%%  view(2);
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Amplitudes (log scale)');
%%  colorbar();
%%  hold on;
%%  %plot(1:1000, 350*ones(1, 1000), 'Color', 'r');
%%  %plot(1:1000, 450*ones(1, 1000), 'Color', 'r');
%%  %%  saveas(gcf, 'Example46_surf_amplitudeDin_log', 'png');
%%  %%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, real(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic (Gaussian Beam)');
%%  xlabel('t (s)');
%%  saveas(gcf, 'Example46_amplitudeGB_log', 'epsc');
%%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, imag(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic - imag');
%%  xlabel('t (s)');
%%  %saveas(gcf, 'Example46_amplitudeODE_log_caustic.fig');


%==============================
% Time Signals
%==============================
%source = sourceRT;
%source = sourceGB;
% Norms
normRT = max(real(source(1).aForward));
normDRT = max(real(source(4).aForward));
%normRT_2 = max(real(sourceRT(1).aForward));
normKWave = max(sensor_data_low(1, :));
figure;
% Subplot 1
subplot(2, 1, 1);
hold on;
axis([0 2e-5 -1.5 2]);
grid on;
box on;
plot(gridR.tForward, source(1).aForward/normRT, 'Color', 'r', 'LineWidth', 2);
plot(gridR.tForward, source(4).aForward/normDRT, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
plot(kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(gridR.tForward, source(2).aForward/normRT, 'Color', 'g', 'LineWidth', 2);
plot(gridR.tForward, source(5).aForward/normDRT, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');
plot(kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(gridR.tForward, source(3).aForward/normRT, 'Color', 'b', 'LineWidth', 2);
plot(gridR.tForward, source(6).aForward/normDRT, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
plot(kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - bottom', 'DRT - bottom', 'kWave - bottom', 'RT - middle', 'DRT - middle', 'kWave - middle', 'RT - top', 'DRT - top', 'kWave - top');
xlabel('t (s)');
ylabel('Amplitude');
title('Forward data - RT vs DRT vs k-Wave');
% Subplot 2
subplot(2, 1, 2);
hold on;
axis([0 2e-5 -1 1]);
grid on;
box on;
plot(gridR.tForward, source(1).aForward/normRT - sensor_data_low(1, :)/normKWave, 'Color', 'r', 'LineWidth', 2);
plot(gridR.tForward, source(4).aForward/normDRT - sensor_data_low(1, :)/normKWave, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
plot(gridR.tForward, source(2).aForward/normRT - sensor_data_mid(1, :)/normKWave, 'Color', 'g', 'LineWidth', 2);
plot(gridR.tForward, source(5).aForward/normDRT - sensor_data_mid(1, :)/normKWave, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');
plot(gridR.tForward, source(3).aForward/normRT - sensor_data_high(1, :)/normKWave, 'Color', 'b', 'LineWidth', 2);
plot(gridR.tForward, source(6).aForward/normDRT - sensor_data_high(1, :)/normKWave, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
legend('Error RT - bottom', 'Error DRT - bottom', 'Error RT - middle', 'Error DRT - middle', 'Error RT - top', 'Error DRT - top');
xlabel('t (s)');
ylabel('Amplitude');
title('Error - RT vs DRT');
set(gcf, 'pos', [700 700 1200 800]);
%saveas(gcf, 'Example46_signalsRTvsDRT_error', 'png');
%saveas(gcf, 'Example46_signalsRTvsDRT_error.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);





cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex07_noCaustic2D;

close all;
%clear all;

load sensor_data.mat;

%==================================================================================
% Plot results
%==================================================================================

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar = [700 700 375 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% kWave and RT
%==============================
load sensor_data_RT_8e-9;
% Norms
normKWave = normRT; %max(sensor_data_low(1, :));
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Subplot 1
figure;
hold on;
axis([0 2e-5 -.6 1.1]);
grid on;
box on;
plot(gridR.tForward, sensor_RT_low/normRT, 'Color', 'r', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(gridR.tForward, sensor_RT_mid/normRT, 'Color', 'g', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(gridR.tForward, sensor_RT_high/normRT, 'Color', 'b', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - bottom ball', 'kWave - bottom ball','RT - middle ball', 'kWave - middle ball', 'RT - top ball', 'kWave - top ball');
xlabel('t (s)');
ylabel('Amplitude');
set(gcf, 'pos', [700 700 1200 400]);
%saveas(gcf, 'Example07_RTvsKWAVE', 'epsc');
%saveas(gcf, 'Example07_RTvsKWAVE.fig');

%==============================
% ERROR PLOT
% kWave and RT
%==============================
load sensor_data_RT_8e-9;
t_array_dt1 = gridR.tForward;
normKWave = normRT; %max(sensor_data_low(1, :));
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Error vectors
dif_low_dt1 = sensor_RT_low/normRT - spline_low;
dif_mid_dt1 = sensor_RT_mid/normRT - spline_mid;
dif_high_dt1 = sensor_RT_high/normRT - spline_high;

error_low_dt1 = sum(dif_low_dt1.^2)/sum(spline_low.^2);
error_mid_dt1 = sum(dif_mid_dt1.^2)/sum(spline_mid.^2);
error_high_dt1 = sum(dif_high_dt1.^2)/sum(spline_high.^2);

figure;
hold on;
%axis([0 2e-5 -.031 .031]);
grid on;
box on;
plot(gridR.tForward, dif_low_dt1, 'Color', 'r', 'LineWidth', 2);
plot(gridR.tForward, dif_mid_dt1, 'Color', 'g', 'LineWidth', 2);
plot(gridR.tForward, dif_high_dt1, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t (s)');
ylabel('Error');
set(gcf, 'pos', [700 700 1200 400]);
%saveas(gcf, 'Example07_errorRTvsKWAVE_8e-9', 'epsc');
%saveas(gcf, 'Example07_errorRTvsKWAVE_8e-9.fig');

%%  %==============================
%%  % ERROR PLOT
%%  % kWave and RT
%%  %==============================
%%  load sensor_data_RT_4e-9;
%%  t_array_dt2 = gridR.tForward;
%%  normKWave = normRT; %max(sensor_data_low(1, :));
%%  spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
%%  spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
%%  spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);
%%  
%%  % Error vectors
%%  dif_low_dt2 = sensor_RT_low/normRT - spline_low;
%%  dif_mid_dt2 = sensor_RT_mid/normRT - spline_mid;
%%  dif_high_dt2 = sensor_RT_high/normRT - spline_high;
%%  
%%  error_low_dt2 = sum(dif_low_dt2.^2)/sum(spline_low.^2);
%%  error_mid_dt2 = sum(dif_mid_dt2.^2)/sum(spline_mid.^2);
%%  error_high_dt2 = sum(dif_high_dt2.^2)/sum(spline_high.^2);
%%  
%%  figure;
%%  hold on;
%%  axis([0 2e-5 -.031 .031]);
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, dif_low_dt2, 'Color', 'r', 'LineWidth', 2);
%%  plot(gridR.tForward, dif_mid_dt2, 'Color', 'g', 'LineWidth', 2);
%%  plot(gridR.tForward, dif_high_dt2, 'Color', 'b', 'LineWidth', 2);
%%  legend('Bottom ball', 'Middle ball', 'Top ball');
%%  xlabel('t (s)');
%%  ylabel('Error');
%%  set(gcf, 'pos', [700 700 1200 400]);
%%  saveas(gcf, 'Example06_errorRTvsKWAVE_4e-9', 'epsc');
%%  saveas(gcf, 'Example06_errorRTvsKWAVE_4e-9.fig');
%%  
%%  %==============================
%%  % ERROR PLOT
%%  % kWave and RT
%%  %==============================
%%  load sensor_data_RT_2e-9;
%%  t_array_dt3 = gridR.tForward;
%%  normKWave = normRT; %max(sensor_data_low(1, :));
%%  spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
%%  spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
%%  spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);
%%  
%%  % Error vectors
%%  dif_low_dt3 = sensor_RT_low/normRT - spline_low;
%%  dif_mid_dt3 = sensor_RT_mid/normRT - spline_mid;
%%  dif_high_dt3 = sensor_RT_high/normRT - spline_high;
%%  
%%  error_low_dt3 = sum(dif_low_dt3.^2)/sum(spline_low.^2);
%%  error_mid_dt3 = sum(dif_mid_dt3.^2)/sum(spline_mid.^2);
%%  error_high_dt3 = sum(dif_high_dt3.^2)/sum(spline_high.^2);
%%  
%%  figure;
%%  hold on;
%%  axis([0 2e-5 -.031 .031]);
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, dif_low_dt3, 'Color', 'r', 'LineWidth', 2);
%%  plot(gridR.tForward, dif_mid_dt3, 'Color', 'g', 'LineWidth', 2);
%%  plot(gridR.tForward, dif_high_dt3, 'Color', 'b', 'LineWidth', 2);
%%  legend('Bottom ball', 'Middle ball', 'Top ball');
%%  xlabel('t (s)');
%%  ylabel('Error');
%%  set(gcf, 'pos', [700 700 1200 400]);
%%  saveas(gcf, 'Example06_errorRTvsKWAVE_2e-9', 'epsc');
%%  saveas(gcf, 'Example06_errorRTvsKWAVE_2e-9.fig');
%%  
%%  %==============================
%%  % CONVERGENCE PLOT
%%  %==============================
%%  conv_t = [8e-9 4e-9 2e-9];
%%  conv_low  = [error_low_dt1 error_low_dt2 error_low_dt3];
%%  conv_mid  = [error_mid_dt1 error_mid_dt2 error_mid_dt3];
%%  conv_high = [error_high_dt1 error_high_dt2 error_high_dt3];
%%  
%%  figure;
%%  semilogy(conv_t, conv_low, 'Color', 'r', 'LineWidth', 2);
%%  hold on;
%%  grid on;
%%  box on;
%%  axis([0 1e-8 2e-5 2e-3]);
%%  ylabel('Error convergence');
%%  semilogy(conv_t, conv_mid, 'Color', 'g', 'LineWidth', 2);
%%  semilogy(conv_t, conv_high, 'Color', 'b', 'LineWidth', 2);
%%  legend('Bottom ball', 'Middle ball', 'Top ball', 'Location', 'southeast');
%%  saveas(gcf, 'Example06_ErrorConvergence', 'epsc');
%%  saveas(gcf, 'Example06_ErrorConvergence.fig');
%%  
%%  %==============================
%%  % Initial Pressure
%%  %==============================
%%  pressure = imresize(source_low.p0 + source_mid.p0 + source_high.p0, factorResize);
%%  figure;
%%  axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
%%  hold on;
%%  surf(1e3*gridR.xAxis, 1e3*gridR.yAxis, pressure', 'EdgeColor', 'none');
%%  view(2);
%%  %title('Initial Pressure');
%%  xlabel('x [mm]');
%%  ylabel('y [mm]');
%%  %set(gca, 'ytick', []);
%%  %set(gca, 'yticklabel', []);
%%  colorbar();
%%  set(gcf, 'pos', positionYBar);
%%  saveas(gcf, 'Example06_InitialPressure', 'png');
%%  saveas(gcf, 'Example06_InitialPressure.fig');



cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex06_forward2D;

close all;
%clear all;

load sensor_data.mat;

saveFigures = 0;
loadData = 1;
%==================================================================================
% Plot results
%==================================================================================

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 400 600];
positionBar = [700 700 375 600];
positionHorizontal = [700 700 1000 400];
set(0,'DefaultFigurePaperPositionMode','auto');
fontSize = 16;
scaleFactor = 1e6;

%==============================
% kWave and RT
%==============================
if(loadData)
load sensor_data_RT_8e-9;
end
% Norms
normKWave = max(sensor_data_low(1, :));
normRT = normKWave;
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Subplot 1
figure;
hold on;
axis([0 20 -.6 1.1]);
grid on;
box on;
plot(scaleFactor*gridR.tForward, sensor_RT_low/normRT, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(scaleFactor*gridR.tForward, sensor_RT_mid/normRT, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(scaleFactor*gridR.tForward, sensor_RT_high/normRT, 'Color', 'b', 'LineWidth', 2);
plot(scaleFactor*kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('HG - bottom', 'kWave - bottom','HG - middle', 'kWave - middle', 'HG - top', 'kWave - top');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gcf, 'pos', positionHorizontal);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
set(gca,'FontSize',fontSize);
if(saveFigures)
saveas(gcf, 'Example06_RTvsKWAVE', 'epsc');
saveas(gcf, 'Example06_RTvsKWAVE.fig');
end
%==============================
% ERROR PLOT
% kWave and RT
%==============================
if(loadData)
load sensor_data_RT_8e-9;
end
t_array_dt1 = gridR.tForward;
normKWave = normRT; %max(sensor_data_low(1, :));
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Error vectors
dif_low_dt1 = sensor_RT_low/normRT - spline_low;
dif_mid_dt1 = sensor_RT_mid/normRT - spline_mid;
dif_high_dt1 = sensor_RT_high/normRT - spline_high;

% Rel error
error_low_dt1 = sum(dif_low_dt1.^2)/sum(spline_low.^2);
error_mid_dt1 = sum(dif_mid_dt1.^2)/sum(spline_mid.^2);
error_high_dt1 = sum(dif_high_dt1.^2)/sum(spline_high.^2);

% Max error
max_error_low_dt1 = max(abs(dif_low_dt1))*normKWave;
max_error_mid_dt1 = max(abs(dif_mid_dt1))*normKWave;
max_error_high_dt1 = max(abs(dif_high_dt1))*normKWave;

figure;
hold on;
axis([0 20 -.05 .05]);
grid on;
box on;
plot(scaleFactor*gridR.tForward, dif_low_dt1, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_mid_dt1, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_high_dt1, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
if(saveFigures)
saveas(gcf, 'Example06_errorRTvsKWAVE_8e-9', 'epsc');
saveas(gcf, 'Example06_errorRTvsKWAVE_8e-9.fig');
end
%==============================
% ERROR PLOT
% kWave and RT
%==============================
if(loadData)
load sensor_data_RT_4e-9;
end
t_array_dt2 = gridR.tForward;
normKWave = max(sensor_data_low(1, :));
normRT = normKWave; %
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Error vectors
dif_low_dt2 = sensor_RT_low/normRT - spline_low;
dif_mid_dt2 = sensor_RT_mid/normRT - spline_mid;
dif_high_dt2 = sensor_RT_high/normRT - spline_high;

% Rel error
error_low_dt2 = sum(dif_low_dt2.^2)/sum(spline_low.^2);
error_mid_dt2 = sum(dif_mid_dt2.^2)/sum(spline_mid.^2);
error_high_dt2 = sum(dif_high_dt2.^2)/sum(spline_high.^2);

% Max error
max_error_low_dt2 = max(abs(dif_low_dt2))*normKWave;
max_error_mid_dt2 = max(abs(dif_mid_dt2))*normKWave;
max_error_high_dt2 = max(abs(dif_high_dt2))*normKWave;

figure;
hold on;
axis([0 20 -.03 .03]);
grid on;
box on;
plot(scaleFactor*gridR.tForward, dif_low_dt2, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_mid_dt2, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_high_dt2, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
if(saveFigures)
saveas(gcf, 'Example06_errorRTvsKWAVE_4e-9', 'epsc');
saveas(gcf, 'Example06_errorRTvsKWAVE_4e-9.fig');
end
%==============================
% ERROR PLOT
% kWave and RT
%==============================
if(loadData)
load sensor_data_RT_2e-9;
end
t_array_dt3 = gridR.tForward;
normKWave = max(sensor_data_low(1, :));
normRT = normKWave;
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, gridR.tForward);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, gridR.tForward);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, gridR.tForward);

% Error vectors
dif_low_dt3 = sensor_RT_low/normRT - spline_low;
dif_mid_dt3 = sensor_RT_mid/normRT - spline_mid;
dif_high_dt3 = sensor_RT_high/normRT - spline_high;

% Rel error
error_low_dt3 = sum(dif_low_dt3.^2)/sum(spline_low.^2);
error_mid_dt3 = sum(dif_mid_dt3.^2)/sum(spline_mid.^2);
error_high_dt3 = sum(dif_high_dt3.^2)/sum(spline_high.^2);

% Max error
max_error_low_dt3 = max(abs(dif_low_dt3))*normKWave;
max_error_mid_dt3 = max(abs(dif_mid_dt3))*normKWave;
max_error_high_dt3 = max(abs(dif_high_dt3))*normKWave;

figure;
hold on;
axis([0 20 -.02 .02]);
grid on;
box on;
plot(scaleFactor*gridR.tForward, dif_low_dt3, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_mid_dt3, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*gridR.tForward, dif_high_dt3, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
if(saveFigures)
saveas(gcf, 'Example06_errorRTvsKWAVE_2e-9', 'epsc');
saveas(gcf, 'Example06_errorRTvsKWAVE_2e-9.fig');
end
%==============================
% CONVERGENCE PLOT
%==============================
scaleFactorDelta = 1e9;
conv_t = [8e-9 4e-9 2e-9];
conv_low  = [error_low_dt1 error_low_dt2 error_low_dt3];
conv_mid  = [error_mid_dt1 error_mid_dt2 error_mid_dt3];
conv_high = [error_high_dt1 error_high_dt2 error_high_dt3];

figure;
semilogy(scaleFactorDelta*conv_t, conv_low, 'Color', 'r', 'LineWidth', 2);
hold on;
grid on;
box on;
axis([0 10 2e-5 2e-3]);
set(gca, 'ytick', [2e-5 1e-4 1e-3])
set(gca, 'yticklabel', {'2\cdot 10^{-5}', '10^{-4}', '10^{-3}'})
ylabel('REE');
xlabel('\Delta t [ns]');
semilogy(scaleFactorDelta*conv_t, conv_mid, 'Color', 'g', 'LineWidth', 2);
semilogy(scaleFactorDelta*conv_t, conv_high, 'Color', 'b', 'LineWidth', 2);
set(gca,'FontSize',fontSize);
legend('Bottom ball', 'Middle ball', 'Top ball', 'Location', 'southeast');
if(saveFigures)
saveas(gcf, 'Example06_ErrorConvergence', 'epsc');
saveas(gcf, 'Example06_ErrorConvergence.fig');
end

%==============================
% CONVERGENCE PLOT - MAX ERROR
%==============================
scaleFactorDelta = 1e9;
conv_t = [8e-9 4e-9 2e-9];
conv_low  = [max_error_low_dt1 max_error_low_dt2 max_error_low_dt3];
conv_mid  = [max_error_mid_dt1 max_error_mid_dt2 max_error_mid_dt3];
conv_high = [max_error_high_dt1 max_error_high_dt2 max_error_high_dt3];

figure;
semilogy(scaleFactorDelta*conv_t, conv_low, 'Color', 'r', 'LineWidth', 2);
hold on;
box on;
axis([0 10 3e-4 5e-3]);
grid on;
ylabel('L infinity error');
xlabel('\Delta t [ns]');
set(gca, 'ytick', [3e-4 1e-3 5e-3])
set(gca, 'yticklabel', {'3\cdot 10^{-4}', '10^{-3}', '5\cdot 10^{-3}'})
h = gca;
h.YRuler.MinorTick =  [3e-4:1e-4:9e-4 2e-3:1e-3:5e-3];
%set(gca, 'yticklabels', [3e-4 1e-3 5e-3])
semilogy(scaleFactorDelta*conv_t, conv_mid, 'Color', 'g', 'LineWidth', 2);
semilogy(scaleFactorDelta*conv_t, conv_high, 'Color', 'b', 'LineWidth', 2);
set(gca,'FontSize',fontSize);
legend('Bottom ball', 'Middle ball', 'Top ball', 'Location', 'southeast');
if(saveFigures)
saveas(gcf, 'Example06_ErrorConvergence_max', 'epsc');
saveas(gcf, 'Example06_ErrorConvergence_max.fig');
end
%==============================
% Initial Pressure
%==============================
pressure = imresize(source_low.p0 + source_mid.p0 + source_high.p0, factorResize);
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
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionYBar);
pbaspect([1 2 1])
if(saveFigures)
saveas(gcf, 'Example06_InitialPressure', 'png');
saveas(gcf, 'Example06_InitialPressure.fig');
end

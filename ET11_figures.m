cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex11_forward3D_het;

clear all;
%close all;

saveFigures = 1;
%==================================================
% Dimensions
%==================================================
% Import dimensions 
dim = importdata('dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);
   
scaleFactor = 1e6;
fontSize = 16;
positionHorizontal = [700 700 1000 400];

%==================================================
% Load data
%==================================================
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

% Forward data
sound_speed_matrix = importdata('sound_speed.dat', ' ', 0);
% Initial Pressure
sound_speed = matrix2cube(sound_speed_matrix, Nz);
figure;
imagesc(1e3*y_axis, 1e3*x_axis, sound_speed(:, :, floor(Nz/2)));
view(2);
box on;
axis tight;
xlabel('y [mm]');
ylabel('x [mm]');
colorbar();
set(gca,'FontSize',fontSize);
positionYBar = [700 700 375 600];
set(gcf, 'pos', positionYBar);
pbaspect([1 2 1])
if(saveFigures)
saveas(gcf, 'Example11_SoundSpeed', 'png');
saveas(gcf, 'Example11_SoundSpeed.fig');
end



load sensor_data;
sensor_RT_low_dt1 = importdata(['forwardSignal_RT_low_4e-8.dat'], ' ', 0);
tForward_dt1 = sensor_RT_low_dt1(1, :);
sensor_RT_low_dt1  = sum(sensor_RT_low_dt1(2:end, :), 1)/(size(sensor_RT_low_dt1, 1)-1);
sensor_RT_mid_dt1  = importdata(['forwardSignal_RT_mid_4e-8.dat'], ' ', 0);
sensor_RT_mid_dt1  = sum(sensor_RT_mid_dt1(2:end, :), 1)/(size(sensor_RT_mid_dt1, 1)-1);
sensor_RT_high_dt1 = importdata(['forwardSignal_RT_high_4e-8.dat'], ' ', 0);
sensor_RT_high_dt1 = sum(sensor_RT_high_dt1(2:end, :), 1)/(size(sensor_RT_high_dt1, 1)-1);

sensor_RT_low_dt2 = importdata(['forwardSignal_RT_low_2e-8.dat'], ' ', 0);
tForward_dt2 = sensor_RT_low_dt2(1, :);
sensor_RT_low_dt2  = sum(sensor_RT_low_dt2(2:end, :), 1)/(size(sensor_RT_low_dt2, 1)-1);
sensor_RT_mid_dt2  = importdata(['forwardSignal_RT_mid_2e-8.dat'], ' ', 0);
sensor_RT_mid_dt2  = sum(sensor_RT_mid_dt2(2:end, :), 1)/(size(sensor_RT_mid_dt2, 1)-1);
sensor_RT_high_dt2 = importdata(['forwardSignal_RT_high_2e-8.dat'], ' ', 0);
sensor_RT_high_dt2 = sum(sensor_RT_high_dt2(2:end, :), 1)/(size(sensor_RT_high_dt2, 1)-1);

sensor_RT_low_dt3 = importdata(['forwardSignal_RT_low_1e-8.dat'], ' ', 0);
tForward_dt3 = sensor_RT_low_dt3(1, :);
sensor_RT_low_dt3  = sum(sensor_RT_low_dt3(2:end, :), 1)/(size(sensor_RT_low_dt3, 1)-1);
sensor_RT_mid_dt3  = importdata(['forwardSignal_RT_mid_1e-8.dat'], ' ', 0);
sensor_RT_mid_dt3  = sum(sensor_RT_mid_dt3(2:end, :), 1)/(size(sensor_RT_mid_dt3, 1)-1);
sensor_RT_high_dt3 = importdata(['forwardSignal_RT_high_1e-8.dat'], ' ', 0);
sensor_RT_high_dt3 = sum(sensor_RT_high_dt3(2:end, :), 1)/(size(sensor_RT_high_dt3, 1)-1);

%========================================================================================================================
% INITIAL PRESSURE AND FORWARD SIGNAL
%========================================================================================================================

normKWave = max(sensor_data_high(1, :));
normRT = normKWave;% max(sensor_RT_high_dt1);
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt1);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt1);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt1);
% Subplot 1
figure;
hold on;
axis([0 20 -1.2 1.2]);
grid on;
box on;
plot(scaleFactor*tForward_dt1, sensor_RT_low_dt1/normRT, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt1, sensor_RT_mid_dt1/normRT, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt1, sensor_RT_high_dt1/normRT, 'Color', 'b', 'LineWidth', 2);
plot(scaleFactor*kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(scaleFactor*kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(scaleFactor*kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - bottom', 'RT - middle', 'RT - top', 'kWave - middle', 'kWave - bottom', 'kWave - top');
xlabel('t [\mus]');
ylabel('Amplitude');
set(gcf, 'pos', positionHorizontal);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
set(gca,'FontSize',fontSize);
if(saveFigures)
saveas(gcf, 'Example11_RTvsKWAVE_4e-8', 'epsc');
saveas(gcf, 'Example11_RTvsKWAVE_4e-8.fig');
end

normKWave = max(sensor_data_high(1, :));
normRT = normKWave; %max(sensor_RT_high_dt2);
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt2);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt2);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt2);
% Subplot 1
figure;
hold on;
%axis([0 20 -.6 1.1]);
grid on;
box on;
plot(scaleFactor*tForward_dt2, sensor_RT_low_dt2/normRT, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt2, sensor_RT_mid_dt2/normRT, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt2, sensor_RT_high_dt2/normRT, 'Color', 'b', 'LineWidth', 2);
%plot(scaleFactor*kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
%plot(scaleFactor*kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
%plot(scaleFactor*kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - bottom', 'RT - middle', 'RT - top', 'kWave - middle', 'kWave - bottom', 'kWave - top');
%xlabel('t [\mus]');


normKWave = max(sensor_data_high(1, :));
normRT = normKWave; %max(sensor_RT_high_dt3);
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt3);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt3);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt3);
% Subplot 1
figure;
hold on;
%axis([0 20 -.6 1.1]);
grid on;
box on;
plot(scaleFactor*tForward_dt3, sensor_RT_low_dt3/normRT, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt3, sensor_RT_mid_dt3/normRT, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt3, sensor_RT_high_dt3/normRT, 'Color', 'b', 'LineWidth', 2);
%plot(scaleFactor*kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
%plot(scaleFactor*kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
%plot(scaleFactor*kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - bottom', 'RT - middle', 'RT - top', 'kWave - middle', 'kWave - bottom', 'kWave - top');

%==============================
% ERROR PLOT
% kWave and RT
%==============================
normRT = normKWave; %max(sensor_RT_high_dt1);
% Spline kWave
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt1);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt1);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt1);

% Error vectors
dif_low_dt1 = sensor_RT_low_dt1/normRT - spline_low;
dif_mid_dt1 = sensor_RT_mid_dt1/normRT - spline_mid;
dif_high_dt1 = sensor_RT_high_dt1/normRT - spline_high;

error_low_dt1 = sum(dif_low_dt1.^2)/sum(spline_low.^2);
error_mid_dt1 = sum(dif_mid_dt1.^2)/sum(spline_mid.^2);
error_high_dt1 = sum(dif_high_dt1.^2)/sum(spline_high.^2);

figure;
hold on;
axis([0 20 -.15 .15]);
grid on;
box on;
plot(scaleFactor*tForward_dt1, dif_low_dt1, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt1, dif_mid_dt1, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt1, dif_high_dt1, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
set(gcf, 'pos', positionHorizontal);
if(saveFigures)
saveas(gcf, 'Example11_errorRTvsKWAVE_4e-8', 'epsc');
saveas(gcf, 'Example11_errorRTvsKWAVE_4e-8.fig');
end

%==============================
% ERROR PLOT
% kWave and RT
%==============================
normRT = normKWave; %max(sensor_RT_high_dt2);
% Spline kWave
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt2);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt2);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt2);

% Error vectors
dif_low_dt2 = sensor_RT_low_dt2/normRT - spline_low;
dif_mid_dt2 = sensor_RT_mid_dt2/normRT - spline_mid;
dif_high_dt2 = sensor_RT_high_dt2/normRT - spline_high;

error_low_dt2 = sum(dif_low_dt2.^2)/sum(spline_low.^2);
error_mid_dt2 = sum(dif_mid_dt2.^2)/sum(spline_mid.^2);
error_high_dt2 = sum(dif_high_dt2.^2)/sum(spline_high.^2);

figure;
hold on;
axis([0 20 -.15 .15]);
grid on;
box on;
plot(scaleFactor*tForward_dt2, dif_low_dt2, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt2, dif_mid_dt2, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt2, dif_high_dt2, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
set(gcf, 'pos', positionHorizontal);
if(saveFigures)
saveas(gcf, 'Example11_errorRTvsKWAVE_2e-8', 'epsc');
saveas(gcf, 'Example11_errorRTvsKWAVE_2e-8.fig');
end

%==============================
% ERROR PLOT
% kWave and RT
%==============================
normRT = normKWave; %max(sensor_RT_high_dt3);
% Spline kWave
spline_low = spline(kgrid.t_array, sensor_data_low(1, :)/normKWave, tForward_dt3);
spline_mid = spline(kgrid.t_array, sensor_data_mid(1, :)/normKWave, tForward_dt3);
spline_high = spline(kgrid.t_array, sensor_data_high(1, :)/normKWave, tForward_dt3);

% Error vectors
dif_low_dt3 = sensor_RT_low_dt3/normRT - spline_low;
dif_mid_dt3 = sensor_RT_mid_dt3/normRT - spline_mid;
dif_high_dt3 = sensor_RT_high_dt3/normRT - spline_high;

error_low_dt3 = sum(dif_low_dt3.^2)/sum(spline_low.^2);
error_mid_dt3 = sum(dif_mid_dt3.^2)/sum(spline_mid.^2);
error_high_dt3 = sum(dif_high_dt3.^2)/sum(spline_high.^2);

figure;
hold on;
axis([0 20 -.15 .15]);
grid on;
box on;
plot(scaleFactor*tForward_dt3, dif_low_dt3, 'Color', 'r', 'LineWidth', 2);
plot(scaleFactor*tForward_dt3, dif_mid_dt3, 'Color', 'g', 'LineWidth', 2);
plot(scaleFactor*tForward_dt3, dif_high_dt3, 'Color', 'b', 'LineWidth', 2);
legend('Bottom ball', 'Middle ball', 'Top ball');
xlabel('t [\mus]');
ylabel('Error');
set(gcf, 'pos', positionHorizontal);
set(gca,'FontSize',fontSize);
set(gca, 'xtick', [0 2 4 6 8 10 12 14 16 18 20]);
set(gcf, 'pos', positionHorizontal);
if(saveFigures)
saveas(gcf, 'Example11_errorRTvsKWAVE_1e-8', 'epsc');
saveas(gcf, 'Example11_errorRTvsKWAVE_1e-8.fig');
end


%==============================
% CONVERGENCE PLOT
%==============================
scaleFactorDelta = 1e9;
conv_t = [1e-8 2e-8 4e-8];
conv_low  = [error_low_dt3 error_low_dt2 error_low_dt1];
conv_mid  = [error_mid_dt3 error_mid_dt2 error_mid_dt1];
conv_high = [error_high_dt3 error_high_dt2 error_high_dt1];

figure;
semilogy(scaleFactorDelta*conv_t, conv_low, 'Color', 'r', 'LineWidth', 2);
hold on;
grid on;
box on;
axis([0 50 1e-4 2e-1]);
ylabel('Error convergence');
xlabel('t [ns]');
semilogy(scaleFactorDelta*conv_t, conv_mid, 'Color', 'g', 'LineWidth', 2);
semilogy(scaleFactorDelta*conv_t, conv_high, 'Color', 'b', 'LineWidth', 2);
set(gca,'FontSize',fontSize);
legend('Bottom ball', 'Middle ball', 'Top ball', 'Location', 'southeast');
if(saveFigures)
saveas(gcf, 'Example11_ErrorConvergence', 'epsc');
saveas(gcf, 'Example11_ErrorConvergence.fig');
end



%=========================================================================
% TRAJECTORIES
%=========================================================================
% Import data
filenameData = 'output_data/Trajectory0.dat';
trajectories = importdata(filenameData, ' ', 0);

% Read number of rays and steps
[nSteps nRays] = size(trajectories);
%nRays = floor(nRays/2);
xCoord = trajectories(:, 1:3:nRays);
yCoord = trajectories(:, 2:3:nRays);
zCoord = trajectories(:, 3:3:nRays);

nRays = floor(nRays/3);
factorScale = 1e3;
% Plot the figure
figure;
ax = gca;
ax.GridAlpha = 0.2;
grid on;
axis(factorScale*[0 Ny*dy 0 Nz*dz 0 Nx*dx]);
hold on;
colours = winter(nRays);
for n = 1:1:nRays
    plot3(factorScale*yCoord(:, n), factorScale*zCoord(:, n), factorScale*xCoord(:, n), 'Color', colours(n, :));
end
xlabel('y [mm]');
ylabel('z [mm]');
zlabel('x [mm]');
pbaspect([1 1 2]);
box on;
positionYBar = [700 700 375 600];
set(gca, 'Zdir', 'reverse');
set(gca, 'FontSize', fontSize);
set(gcf, 'pos', positionYBar);
view(0, 0);
if(saveFigures)
saveas(gcf, 'Example11_rays', 'epsc');
saveas(gcf, 'Example11_rays.fig');
end

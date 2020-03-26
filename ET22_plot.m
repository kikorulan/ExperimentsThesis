%================================================================================
% ET21 - PLOT
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex22_vesselsGB;

close all;
clear all;

load initial_pressure_128x128_broadBand;
load initial_pressure_128x128_HF;
load vessels_GBdecomposition;

load sensor_data_kWave;
load sensor_data_GB;
load sensor_data_GBdecomposition.mat
load adjoint_kWave;
load adjoint_GB;
%===================================================
% PARAMETERS
%===================================================
saveResults  = 0;
plot_ip      = 1; 
plot_forward = 1;
plot_adjoint = 1;

% Position 
position    = [700 700 450 420];
positionHor = [700 700 600 420];

% Dimensions
N = 128;
dx = 1e-4;
x_axis = 0:dx:(N-1)*dx;
y_axis = x_axis;
[Y, X] = meshgrid(x_axis, y_axis);
factorS = 1e3;
factorT = 1e6;
%=====================================================================================================
% INITIAL PRESSURE
%=====================================================================================================
if(plot_ip)
figure;
surf(factorS*x_axis, factorS*y_axis, u0, 'EdgeColor', 'none');
view(2);
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if (saveResults)
saveas(gcf, 'ET22_initialPressure_broadband', 'png')
end

figure;
surf(factorS*x_axis, factorS*y_axis, u0_HF, 'EdgeColor', 'none');
view(2);
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if (saveResults)
saveas(gcf, 'ET22_initialPressure_HF', 'png')
end
end

%=====================================================================================================
% MULTIPLE GB APPROX - FORWARD
%=====================================================================================================
if (plot_forward)
GBF = vessels_GBdecomposition;
normGB = max(real(GBF(:)));
normIP = max(real(u0_HF(:)));

%==================================================
%=======    DECOMPOSITION INITIAL PRESSURE
%==================================================
% GB approximation - Real part
figure;
surf(factorS*y_axis, factorS*x_axis, real(GBF)/normGB, 'EdgeColor', 'none');
view(2)
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET22_vessels_GBapprox_real', 'png')
end

% GB approximation - Imaginary part
figure;
surf(factorS*y_axis, factorS*x_axis, imag(GBF)/normGB, 'EdgeColor', 'none');
view(2)
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET22_vessels_GBapprox_imag', 'png')
end

% GB approximation - Error
figure;
surf(factorS*y_axis, factorS*x_axis, real(u0_HF)/normIP-real(GBF)/normGB, 'EdgeColor', 'none');
view(2)
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET22_vessels_GBapprox_error', 'png')
end

%==================================================
%=======    FORWARD DATA
%==================================================
pos = [0 0 1000 800];
normKW = max(sensor_data(:));
normGB = max(real(sensor_data_GB(:)));


figure;
surf(factorT*t_array_GB, factorS*x_axis, real(sensor_data), 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1.45 1 1])
xlabel('t [\mus]')
ylabel('y [mm]')
box on;
set(gca, 'FontSize', 15)
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_kWave', 'png')
end

% Full set of sensors
figure;
surf(factorT*t_array_GB, factorS*x_axis, real(sensor_data_GB)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1.45 1 1])
xlabel('t [\mus]')
ylabel('y [mm]')
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_GB', 'png')
end


% Difference
figure;
surf(factorT*t_array_GB, factorS*x_axis, real(sensor_data)-real(sensor_data_GB)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1.45 1 1])
xlabel('t [\mus]')
ylabel('y [mm]')
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_errorGB', 'png')
end

%==================================================
%=======    PLOT BY SENSOR
%==================================================
% SENSOR 20
figure;
plot(factorT*t_array_GB, real(sensor_data_GB(20, :))/normGB*normKW, 'LineWidth', 2);
hold on;
plot(factorT*t_array_GB, real(sensor_data(20, :)), 'LineWidth', 2)
box on;
grid on;
axis([0 10 -0.08 0.08]);
pbaspect([1.45 1 1]);
legend('GB', 'k-Wave');
xlabel('t [\mus]')
ylabel('Amplitude')
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_sensor20', 'png')
end

% SENSOR 40
figure;
plot(factorT*t_array_GB, real(sensor_data_GB(40, :))/normGB*normKW, 'LineWidth', 2);
hold on;
plot(factorT*t_array_GB, real(sensor_data(40, :)), 'LineWidth', 2)
box on;
grid on;
axis([0 10 -0.08 0.08])
pbaspect([1.45 1 1]);
legend('GB', 'k-Wave');
xlabel('t [\mus]')
ylabel('Amplitude')
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_sensor40', 'png')
end

% SENSOR 60
figure;
plot(factorT*t_array_GB, real(sensor_data_GB(60, :))/normGB*normKW, 'LineWidth', 2);
hold on;
plot(factorT*t_array_GB, real(sensor_data(60, :)), 'LineWidth', 2)
box on;
grid on;
axis([0 10 -0.08 0.08])
legend('GB', 'k-Wave');
xlabel('t [\mus]')
ylabel('Amplitude')
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_sensor60', 'png')
end

%==================================================
%=======    ERROR BY SENSOR
%==================================================
error_sensor = zeros(size(x_axis));
energy_sensor = zeros(size(x_axis));
for ii = 1:N
    error_sensor(ii) = sum(real(sensor_data(ii, :)/normKW - sensor_data_GB(ii, :)/normGB).^2);
    energy_sensor(ii) = sum(real(sensor_data(ii, :)/normKW).^2);
end
% Plot error and energy
figure;
semilogy(factorS*x_axis, error_sensor, 'LineWidth', 2);
hold on;
%semilogy(x_axis, error_sensor./max_energy_sensor, 'LineWidth', 2);
semilogy(factorS*x_axis, error_sensor./max(energy_sensor, 1e-20), 'LineWidth', 2);
semilogy(factorS*x_axis, energy_sensor, 'LineWidth', 2);
legend('Error', 'REE', 'Energy');
grid on;
axis([0 factorS*(N-1)*dx 1e-2 100])
xlabel('y [mm]')
ylabel('Error')
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET22_sensor_data_errorSensor', 'png')
end


end


%=====================================================================================================
% MULTIPLE GB APPROX - ADJOINT
%=====================================================================================================
if (plot_adjoint)
%==================================================
%=======    DECOMPOSITION FORWARD DATA
%==================================================
GBF = sensor_data_GBdecomposition;
normKW = max(real(sensor_data(:)));
normGB = max(real(GBF(:)));
% Plot
figure;
surf(factorT*t_array_GB, factorS*x_axis, real(GBF)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1.45 1 1]);
xlabel('t [\mus]');
ylabel('y [mm]');
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if (saveResults)
saveas(gcf, 'ET22_sensor_data_GBdecomposition', 'png');
end

figure;
surf(factorT*t_array_GB, factorS*x_axis, real(sensor_data)-real(GBF)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1.45 1 1]);
xlabel('t [\mus]')
ylabel('y [mm]')
box on
set(gca, 'FontSize', 15);
set(gcf, 'pos', positionHor);
if (saveResults)
saveas(gcf, 'ET22_sensor_data_GBdecomposition_error', 'png');
end

%==================================================
%=======    ADJOINT
%==================================================
adjoint_kWave = adjoint_kWave.p_final;
normGB = max(real(adjoint_GB(:)));
normKW = max(real(adjoint_kWave(:)));

figure;
surf(factorS*y_axis, factorS*x_axis, real(adjoint_GB)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
box on;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
pbaspect([1 1 1]);
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
caxis([-0.16 0.16]);
if (saveResults)
saveas(gcf, 'ET22_adjoint_GB', 'png');
end

figure;
surf(factorS*y_axis, factorS*x_axis, real(adjoint_kWave), 'EdgeColor', 'none');
view(2);
axis tight;
box on;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
pbaspect([1 1 1]);
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
caxis([-0.16 0.16]);
if (saveResults)
saveas(gcf, 'ET22_adjoint_kWave', 'png');
end

figure;
surf(factorS*y_axis, factorS*x_axis, real(adjoint_kWave)-real(adjoint_GB)/normGB*normKW, 'EdgeColor', 'none');
view(2);
axis tight;
box on;
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
pbaspect([1 1 1]);
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
caxis([-0.16 0.16]);
if (saveResults)
saveas(gcf, 'ET22_adjoint_error', 'png');
end
end

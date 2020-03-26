%================================================================================
% ET21 - PLOT
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex21_singleGB;

close all;
clear all;

load singleGB;
load singleGB_approx;
load singleGB_GBdecomposition;
load sensor_data_singleGB_kWave;
load sensor_data_singleGB_approx;
load sensor_data_singleGB_GB;
load sensor_data_singleGB_GBparameters.mat
load adjoint_singleGB_kWave;
load adjoint_singleGB_GB;
%===================================================
% PARAMETERS
%===================================================
saveResults   = 1;
plot_IP       = 0;
plot_singleGB = 0;
plot_forward  = 0;
plot_adjoint  = 1;
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
%===================================================
% INITIAL PRESSURE
%===================================================
if(plot_IP)
figure;
surf(factorS*x_axis, factorS*y_axis, real(singleGB), 'EdgeColor', 'none');
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
saveas(gcf, 'ET21_singleGB', 'png')
end
end
%=====================================================================================================
% SINGLE GB APPROX
%=====================================================================================================
if(plot_singleGB)
figure;
surf(factorS*x_axis, factorS*y_axis, real(2*singleGB_approx), 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
box on;
set(gca, 'FontSize', 15)
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET21_singleGB_approx', 'png')
end


figure;
surf(factorS*x_axis, factorS*y_axis, real(singleGB-2*singleGB_approx), 'EdgeColor', 'none');
view(2);
axis tight;
colorbar();
pbaspect([1 1 1])
xlabel('x [mm]')
ylabel('y [mm]')
box on;
set(gca, 'FontSize', 15)
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET21_singleGB_approx_error', 'png')
end

%==================================================
%=======    FORWARD DATA                           
%==================================================
GBF = sensor_data_singleGB_approx;

figure;
surf(factorT*t_array, factorS*x_axis, real(sensor_data), 'EdgeColor', 'none');
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
saveas(gcf, 'ET21_singleGB_sensor_data_kWave', 'png')
end

figure;
surf(factorT*t_array, factorS*x_axis, real(GBF), 'EdgeColor', 'none');
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
saveas(gcf, 'ET21_singleGB_sensor_data_approx', 'png')
end

figure;
surf(factorT*t_array, factorS*x_axis, real(sensor_data-GBF), 'EdgeColor', 'none');
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
saveas(gcf, 'ET21_singleGB_sensor_data_errorApprox', 'png')
end

%==================================================
%=======    ERROR BY SENSOR            
%==================================================
error_sensor = zeros(size(x_axis));
energy_sensor = zeros(size(x_axis));
max_energy_sensor = zeros(size(x_axis));
for ii = 1:N
    error_sensor(ii) = sum(real(sensor_data(ii, :) - GBF(ii, :)).^2);
    energy_sensor(ii) = sum(real(sensor_data(ii, :)).^2);
    max_energy_sensor(:) = max(energy_sensor(:));
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
axis([0 factorS*(N-1)*dx 1e-3 100])
xlabel('y [mm]')
ylabel('Error')
set(gca, 'FontSize', 15);
set(gca, 'xTick', 0:2:12);
if(saveResults)
saveas(gcf, 'ET21_singleGB_sensor_data_approx_errorSensor', 'png')
end

end

%=====================================================================================================
% MULTIPLE GB APPROX - FORWARD
%=====================================================================================================
if (plot_forward)
GBF = singleGB_GBdecomposition;
normGB = max(real(GBF(:)));
normIP = max(real(singleGB(:)));

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
saveas(gcf, 'ET21_singleGB_GBapprox_real', 'png')
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
saveas(gcf, 'ET21_singleGB_GBapprox_imag', 'png')
end

% GB approximation - Error
figure;
surf(factorS*y_axis, factorS*x_axis, real(singleGB)/normIP-real(GBF)/normGB, 'EdgeColor', 'none');
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
saveas(gcf, 'ET21_singleGB_GBapprox_error', 'png')
end

% GB approximation - Coefficients
[h1, h2] = gb1.plot_coeffs();
set(0, 'CurrentFigure', h1);
view(2)
colorbar();
pbaspect([1 1 1])
xlabel('index x')
ylabel('index y')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET21_singleGB_GBapprox_coeffs1', 'png')
end

set(0, 'CurrentFigure', h2);
view(2)
colorbar();
caxis([-2 2]);
pbaspect([1 1 1])
xlabel('index x')
ylabel('index y')
axis tight;
box on;
set(gca, 'FontSize', 15);
set(gcf, 'pos', position);
if(saveResults)
saveas(gcf, 'ET21_singleGB_GBapprox_coeffs2', 'png')
end


%==================================================
%=======    FORWARD DATA
%==================================================
pos = [0 0 1000 800];
normKW = max(sensor_data(:));
normGB = max(real(sensor_data_GB(:)));
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
set(gca, 'FontSize', 15)
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET21_singleGB_sensor_data_GB', 'png')
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
set(gca, 'FontSize', 15)
set(gcf, 'pos', positionHor);
if(saveResults)
saveas(gcf, 'ET21_singleGB_sensor_data_errorGB', 'png')
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
axis([0 factorS*(N-1)*dx 1e-3 100])
xlabel('y [mm]')
ylabel('Error')
set(gca, 'FontSize', 15);
set(gca, 'xTick', 0:2:12);
if(saveResults)
saveas(gcf, 'ET21_singleGB_sensor_data_GB_errorSensor', 'png')
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
saveas(gcf, 'ET21_singleGB_sensor_data_GBdecomposition', 'png');
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
saveas(gcf, 'ET21_singleGB_sensor_data_GBdecomposition_error', 'png');
end

%==================================================
%=======    ADJOINT
%==================================================
adjoint_kWave = adjoint_singleGB_kWave.p_final;
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
if (saveResults)
saveas(gcf, 'ET21_singleGB_adjoint_GB', 'png');
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
if (saveResults)
saveas(gcf, 'ET21_singleGB_adjoint_kWave', 'png');
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
if (saveResults)
saveas(gcf, 'ET21_singleGB_adjoint_error', 'png');
end



end

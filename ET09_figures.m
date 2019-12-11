
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex09_adjoint2D_het;

close all;
%clear all;

load sensor_data.mat;
load adjoint_data.mat;

saveFigures = 1;

%==================================================================================
% Plot results
%==================================================================================

position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar = [700 700 375 600];
positionY = [700 700 340 600];
positionHorizontal = [700 700 1000 400];
set(0,'DefaultFigurePaperPositionMode','auto');
fontSize = 16;
scaleFactor = 1e3;

figure;
%================================================================================
% ERROR PLOT
% kWave and RT
%================================================================================

load adjoint_pressure_RT_8e-9;
p_kWave_low = imresizen(adjoint_pressure_low.p_final', factorResize);
p_kWave_mid = imresizen(adjoint_pressure_mid.p_final', factorResize);
p_kWave_high = imresizen(adjoint_pressure_high.p_final', factorResize);

%==============================
% DT = 8e-9
%==============================
load adjoint_pressure_RT_8e-9;
p_RT_low = adjoint_pressure_low_RT';
p_RT_mid = adjoint_pressure_mid_RT';
p_RT_high = adjoint_pressure_high_RT';

dif_low_dt1 = p_kWave_low - p_RT_low;
dif_mid_dt1 = p_kWave_mid - p_RT_mid;
dif_high_dt1 = p_kWave_high - p_RT_high;

error_low_dt1 = sum(dif_low_dt1(:).^2)/sum(p_kWave_low(:).^2);
error_mid_dt1 = sum(dif_mid_dt1(:).^2)/sum(p_kWave_mid(:).^2);
error_high_dt1 = sum(dif_high_dt1(:).^2)/sum(p_kWave_high(:).^2);

% % kWave
% figure;
% surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, p_kWave_low + p_kWave_mid + p_kWave_high, 'EdgeColor', 'none');
% view(2);
% box on;
% axis tight;
% xlabel('x [mm]');
% ylabel('y [mm]');
% pbaspect([1 2 1]);
% set(gca,'FontSize',fontSize);
% set(gcf, 'pos', positionY);
% if(saveFigures)
% saveas(gcf, 'Example09_adjoint_kWave', 'png');
% saveas(gcf, 'Example09_adjoint_kWave.fig');
% end
% 
% % RT
% figure;
% surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, p_RT_low + p_RT_mid + p_RT_high, 'EdgeColor', 'none');
% view(2);
% box on;
% axis tight;
% xlabel('x [mm]');
% pbaspect([1 2 1]);
% colorbar();
% caxis([min(p_kWave_low(:)) max(p_kWave_low(:))]);
% set(gca,'FontSize',fontSize);
% set(gcf, 'pos', positionBar);
% 
% 
% % Difference
% figure;
% surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, dif_low_dt1 + dif_mid_dt1 + dif_high_dt1, 'EdgeColor', 'none');
% view(2);
% box on;
% axis tight;
% xlabel('x [mm]');
% pbaspect([1 2 1]);
% colorbar();
% caxis([-max(p_RT_low(:))/10 max(p_RT_low(:))/10]);
% set(gca,'FontSize',fontSize);
% set(gcf, 'pos', positionBar);
% 
% % Line plot
% line_kWave = p_kWave_low(:, 256);
% line_RT = p_RT_low(:, 256);
% figure;
% plot(line_kWave, 'Color', 'b', 'LineWidth', 2);
% hold on;
% plot(line_RT, 'Color', 'r');
% grid on;
% plot(line_RT-line_kWave, 'Color', 'm', 'LineWidth', 2);
% grid on;
% set(gca,'FontSize',fontSize);
% xlabel('pixel');
% ylabel('Error');

%==============================
% DT = 4e-9
%==============================
load adjoint_pressure_RT_4e-9;
p_RT_low = adjoint_pressure_low_RT';
p_RT_mid = adjoint_pressure_mid_RT';
p_RT_high = adjoint_pressure_high_RT';

dif_low_dt2 = p_kWave_low - p_RT_low;
dif_mid_dt2 = p_kWave_mid - p_RT_mid;
dif_high_dt2 = p_kWave_high - p_RT_high;

error_low_dt2 = sum(dif_low_dt2(:).^2)/sum(p_kWave_low(:).^2);
error_mid_dt2 = sum(dif_mid_dt2(:).^2)/sum(p_kWave_mid(:).^2);
error_high_dt2 = sum(dif_high_dt2(:).^2)/sum(p_kWave_high(:).^2);

%==============================
% DT = 2e-9
%==============================
load adjoint_pressure_RT_2e-9;
p_RT_low = adjoint_pressure_low_RT';
p_RT_mid = adjoint_pressure_mid_RT';
p_RT_high = adjoint_pressure_high_RT';
p_RT_all = adjoint_pressure_high_RT';

dif_low_dt3 = p_kWave_low - p_RT_low;
dif_mid_dt3 = p_kWave_mid - p_RT_mid;
dif_high_dt3 = p_kWave_high - p_RT_high;

error_low_dt3 = sum(dif_low_dt3(:).^2)/sum(p_kWave_low(:).^2);
error_mid_dt3 = sum(dif_mid_dt3(:).^2)/sum(p_kWave_mid(:).^2);
error_high_dt3 = sum(dif_high_dt3(:).^2)/sum(p_kWave_high(:).^2);

% kWave
figure;
surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, p_kWave_low + p_kWave_mid + p_kWave_high, 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('x [mm]');
ylabel('y [mm]');
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
if(saveFigures)
saveas(gcf, 'Example09_adjoint_kWave', 'png');
saveas(gcf, 'Example09_adjoint_kWave.fig');
end

% RT
figure;
surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, p_RT_low + p_RT_mid + p_RT_high, 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('x [mm]');
pbaspect([1 2 1]);
colorbar();
caxis([min(p_kWave_low(:)) max(p_kWave_low(:))]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example09_adjoint_RT', 'png');
saveas(gcf, 'Example09_adjoint_RT.fig');
end

% Difference
figure;
surf(scaleFactor*gridR.xAxis, scaleFactor*gridR.yAxis, dif_low_dt3 + dif_mid_dt3 + dif_high_dt3, 'EdgeColor', 'none');
view(2);
box on;
axis tight;
xlabel('x [mm]');
pbaspect([1 2 1]);
colorbar();
caxis([-max(p_RT_low(:))/10 max(p_RT_low(:))/10]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example09_error_2e-9', 'png');
saveas(gcf, 'Example09_error_2e-9.fig');
end


% Line plot
line_kWave = p_kWave_mid(:, 256);
line_RT = p_RT_mid(:, 256);
figure;
% Subplot 1
subplot(2, 1, 1);
plot(line_kWave, 'Color', 'b', 'LineWidth', 2);
hold on;
plot(line_RT, 'Color', 'r');
axis([1 1024 -1e-4 3e-4]);
legend('kWave', 'RT', 'location', 'northwest');
grid on;
set(gca,'FontSize',fontSize);
xlabel('pixel');
ylabel('Amplitude');
% Subplot 2
subplot(2, 1, 2);
plot(line_RT-line_kWave, 'Color', 'm', 'LineWidth', 2);
axis([1 1024 -5e-5 5e-5]);
grid on;
set(gca,'FontSize',fontSize);
set(gcf, 'pos', [700 700 500 400]);
xlabel('pixel');
ylabel('Error');
if(saveFigures)
saveas(gcf, 'Example09_LineError_2e-9', 'epsc');
saveas(gcf, 'Example09_LineError_2e-9.fig');
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
axis([0 10 1e-5 2e-1]);
ylabel('Error convergence');
xlabel('t [ns]');
semilogy(scaleFactorDelta*conv_t, conv_mid, 'Color', 'g', 'LineWidth', 2);
semilogy(scaleFactorDelta*conv_t, conv_high, 'Color', 'b', 'LineWidth', 2);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', [700 700 500 400]);
legend('Bottom ball', 'Middle ball', 'Top ball', 'Location', 'southeast');
if(saveFigures)
saveas(gcf, 'Example09_ErrorConvergence', 'epsc');
saveas(gcf, 'Example09_ErrorConvergence.fig');
end

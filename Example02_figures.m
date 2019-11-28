

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex02_simulation2D_homo;

%close all;
clear all;

drawForward = 1;
drawAdjoint = 1;
drawMix = 1;


%===============================================================================================================
%===============================================================================================================
%=====================                 FORWARD PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawForward)

load gridRT_impulse;
load signalRT;
load sensor_data_kWave;

%===============================================================================================================
% SINOGRAM
%===============================================================================================================
nSources = 764;
dcol_pos = @(x) [x(:, end) x(:, 1:end-1)];
dcol_neg = @(x) [x(:, 2:end) x(:, 1)];

positionYNoBar     = [700 700 550 630];
positionNoYBar     = [700 700 600 630];
positionYBar       = [700 700 620 630];
positionNoYNoBar   = [700 700 530 630];
%positionNoYBar = [700 700 610 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% SORT DATA
%========================================
Nx = 128;
Ny = 256;
% k-Wave 
unsorted_data = sensor_data;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
sensor_data = sort_data; 

% RT smooth
unsorted_data = signalRT;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalRT = sort_data; 

%========================================
% NORMS
%========================================
% k-Wave
maxKW = max(sensor_data(:));
normKW = maxKW; 
signalKW_norm = sensor_data/normKW;
% RT smooth
maxRT = max(signalRT(:));
normRT = maxRT;
signalRT_norm = signalRT/normRT;

%========================================
% FIGURES - FORWARD
%========================================
fontSize = 16;
% k-Wave
figure;
surf(1e6*kgrid.t_array, 1:nSources, signalKW_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
%colorbar();
caxis([-0.4 1]);
xlabel('t [\mus]');
ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionYNoBar);
saveas(gcf, 'Example02_kWave_f2D.fig');
saveas(gcf, 'Example02_kWave_f2D', 'png');


% RT smooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.4 1]);
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example02_RT_f2D.fig');
saveas(gcf, 'Example02_RT_f2D', 'png');

%========================================
% FIGURES - ERROR
%========================================
% Error - RT smoot
signal = dcol_neg(signalRT_norm);
error_RT = signal - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_RT, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.5 0.5]);
xlabel('t [\mus]');
ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionYNoBar);
saveas(gcf, 'Example02_error_RT_f2D.fig');
saveas(gcf, 'Example02_error_RT_f2D', 'png');

end

%===============================================================================================================
%===============================================================================================================
%=====================                 ADJOINT PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawAdjoint)

load adjoint_kWave.mat;
load adjoint_RT;
load gridRT_impulse;

%========================================
% Draw parameters
%========================================
axisGrid = [0 1e3*(Rgrid.Nx-1)*Rgrid.dx 0 1e3*(Rgrid.Ny-1)*Rgrid.dy];
position     = [700 700 320 630];
positionY    = [700 700 340 630];
positionBar  = [700 700 380 630];
positionYBar = [700 700 410 630];
fontSize = 15;
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% kWave Adjoint
%========================================
pixelAReverse = pressure_adjoint_kWave.p_final;
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
saveas(gcf, 'Example02_kWave_adjoint', 'png');
saveas(gcf, 'Example02_kWave_adjoint.fig');


%========================================
% RT 
%========================================
pixelAReverse = adjoint_RT;
maxPixelRT = max(real(pixelAReverse(:)));
pixelRT = max(0, real(pixelAReverse)/maxPixelRT);
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example02_RT_adjoint', 'png');
saveas(gcf, 'Example02_RT_adjoint.fig');


%========================================
% Error
%========================================
adjointKW_norm = max(0, pressure_adjoint_kWave.p_final/max(pressure_adjoint_kWave.p_final(:)));
adjointRT_norm = max(0, adjoint_RT/max(adjoint_RT(:)));
pixel = adjointKW_norm - adjointRT_norm;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixel', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example02_error_RT', 'png');
saveas(gcf, 'Example02_error_RT.fig');

end


%===============================================================================================================
%===============================================================================================================
%=====================             MIX ADJOINT PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawMix)
    
load adjointKWaveForward_RT;
load adjointRTForward_kWave;

%========================================
% Draw parameters
%========================================
axisGrid = [0 1e3*(Rgrid.Nx-1)*Rgrid.dx 0 1e3*(Rgrid.Ny-1)*Rgrid.dy];
position     = [700 700 320 630];
positionY    = [700 700 340 630];
positionBar  = [700 700 380 630];
positionYBar = [700 700 410 630];
fontSize = 15;
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% RT Adjoint - kWave data
%========================================
pixelAReverse = real(adjointKWaveForward_RT);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example02_RT_adjoint_kWave_data', 'png');
saveas(gcf, 'Example02_RT_adjoint_kWave_data.fig');

%========================================
% kWave Adjoint - RT data
%========================================
pixelAReverse = real(adjointRTForward_kWave.p_final);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
saveas(gcf, 'Example02_kWave_adjoint_RT_data', 'png');
saveas(gcf, 'Example02_kWave_adjoint_RT_data.fig');

%========================================
% Error mix kWave Forward - RT recon
%========================================
adjointKW = max(0, pressure_adjoint_kWave.p_final/max(pressure_adjoint_kWave.p_final(:)));
adjointRT = max(0, real(adjointKWaveForward_RT)/max(real(adjointKWaveForward_RT(:))));
pixelKWave = adjointKW - adjointRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example02_mix_error_RT', 'png');
saveas(gcf, 'Example02_mix_error_RT.fig');


end



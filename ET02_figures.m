

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex02_simulation2D_homo;

close all;
%clear all;

drawForward = 0;
drawAdjoint = 1;
drawMix = 1;
saveFigures = 0;


%===============================================================================================================
%===============================================================================================================
%=====================                 FORWARD PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawForward)

%load gridRT_impulse;
load signalRT;
load sensor_data_kWave;

%=============================
% Initial Pressure
%=============================
fontSize = 16;
position     = [700 700 300 600];
positionY    = [700 700 320 600];
positionBar  = [700 700 363 600];
positionYBar = [700 700 390 600];


load initial_pressure_nonsmooth;
X = 0:Rgrid.dx:(Rgrid.Nx-1)*Rgrid.dx;
Y = 0:Rgrid.dy:(Rgrid.Ny-1)*Rgrid.dy;
figure;
surf(1e3*X, 1e3*Y, initial_pressure_nonsmooth', 'EdgeColor', 'none');
axis tight;
box on;
pbaspect([1 2 1]);
xlabel('x [mm]');
ylabel('y [mm]');
view(2);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
if(saveFigures)
saveas(gcf, 'Example02_initialPressure_nonsmooth', 'png'); 
saveas(gcf, 'Example02_initialPressure_nonsmooth.fig'); 
end

load initial_pressure;
figure;
surf(1e3*X, 1e3*Y, initial_pressure', 'EdgeColor', 'none');
axis tight;
box on;
pbaspect([1 2 1]);
colorbar();
xlabel('x [mm]');
view(2);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_initialPressure', 'png'); 
saveas(gcf, 'Example02_initialPressure.fig'); 
end


%===============================================================================================================
% SINOGRAM
%===============================================================================================================
nSources = 764;
dcol_pos = @(x) [x(:, end) x(:, 1:end-1)];
dcol_neg = @(x) [x(:, 2:end) x(:, 1)];

position     = [700 700 530 630];
positionY    = [700 700 550 630];
positionBar  = [700 700 620 630];
positionYBar = [700 700 620 630];
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
normRT = maxKW;
signalRT_norm = signalRT/normRT;

%========================================
% FIGURES - FORWARD
%========================================
fontSize = 16;
% k-Wave
figure;
surf(1e6*kgrid.t_array, 1:nSources, signalKW_norm, 'EdgeColor', 'none');
axis([0 25 1 nSources]);
view(2);
box on;
%colorbar();
caxis([-.4 1]);
xlabel('t [\mus]');
ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
if(saveFigures)
saveas(gcf, 'Example02_kWave_f2D.fig');
saveas(gcf, 'Example02_kWave_f2D', 'png');
end

% RT smooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_norm, 'EdgeColor', 'none');
axis([0 25 1 nSources]);
view(2);
box on;
caxis([-.4 1]);
colorbar();
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_RT_f2D.fig');
saveas(gcf, 'Example02_RT_f2D', 'png');
end

%========================================
% FIGURES - ERROR
%========================================
% Error 
error_RT = signalRT_norm - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_RT, 'EdgeColor', 'none');
axis([0 25 1 nSources]);
view(2);
box on;
caxis([-.1 .1]);
colorbar();
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf,'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_error_RT_f2D.fig');
saveas(gcf, 'Example02_error_RT_f2D', 'png');
end

REE = sum(error_RT(:).^2)/sum(signalKW_norm(:).^2)

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
position     = [700 700 300 600];
positionY    = [700 700 340 600];
positionBar  = [700 700 363 600];
positionYBar = [700 700 390 600];

fontSize = 16;
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% kWave Adjoint
%========================================
pixelAReverse = pressure_adjoint_kWave.p_final;
normKW = max(pixelAReverse(:));
pixelKW = pixelAReverse/normKW;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKW', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
caxis([-.1 1]);
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
if(saveFigures)
saveas(gcf, 'Example02_kWave_adjoint', 'png');
saveas(gcf, 'Example02_kWave_adjoint.fig');
end

%========================================
% RT 
%========================================
pixelAReverse = adjoint_RT;
pixelRT = pixelAReverse/normKW;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.1 1]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_RT_adjoint', 'png');
saveas(gcf, 'Example02_RT_adjoint.fig');
end

%========================================
% Error
%========================================
pixel = pixelKW - pixelRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixel', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.1 .1]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_error_RT', 'png');
saveas(gcf, 'Example02_error_RT.fig');
end
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
position     = [700 700 300 600];
positionY    = [700 700 340 600];
positionBar  = [700 700 363 600];
positionYBar = [700 700 390 600];

fontSize = 16;
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% kWave Adjoint - RT data
%========================================
pixelAReverse = real(adjointRTForward_kWave.p_final);
%normKW = max(pixelAReverse(:));
pixelKW2 = pixelAReverse/normKW;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKW2', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
if(saveFigures)
saveas(gcf, 'Example02_kWave_adjoint_RT_data', 'png');
saveas(gcf, 'Example02_kWave_adjoint_RT_data.fig');
end
%========================================
% RT Adjoint - kWave data
%========================================
pixelAReverse = real(adjointKWaveForward_RT);
pixelRT = pixelAReverse/normKW;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_RT_adjoint_kWave_data', 'png');
saveas(gcf, 'Example02_RT_adjoint_kWave_data.fig');
end
%========================================
% Error mix kWave Forward - RT recon
%========================================
errorKW = pixelKW2 - pixelRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, errorKW', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.1 .1]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
pbaspect([1 2 1]);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
if(saveFigures)
saveas(gcf, 'Example02_mix_error_RT', 'png');
saveas(gcf, 'Example02_mix_error_RT.fig');
end

%e1 = pixelKW2 - pixelKW;
%REE = sum(e1(:).^2)/sum(pixelKW(:).^2)
e2 = pixelRT - pixelKW;
REE = sum(e2(:).^2)/sum(pixelKW(:).^2)

end



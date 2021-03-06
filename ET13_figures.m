cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex13_adjoint3D_het;

clear all;
close all;

saveFigures = 0;
load adjoint_pressure;
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
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;


%==================================================
% LOAD Adjoint - LOW
%==================================================
% Adjoint low
adjoint_kWave_low       = adjoint_pressure_low;
adjoint_kWave_low_slice = adjoint_pressure_low(:, :, floor(Nz/2));
normKWave_low           = sqrt(sum(adjoint_kWave_low(:).^2));
normKWave_low_slice     = sqrt(sum(adjoint_kWave_low_slice(:).^2));
normKWave_low_max       = max(adjoint_kWave_low(:));

% Adjoint RT low - dt1
adjoint_RT_low_matrix    = importdata('adjoint_pressure_RT_low_4e-8.dat', ' ', 0);
adjoint_RT_low_dt1       = matrix2cube(adjoint_RT_low_matrix, Nz);
adjoint_RT_low_dt1_slice = adjoint_RT_low_dt1(:, :, floor(Nz/2));
% Adjoint RT low - dt2
adjoint_RT_low_matrix    = importdata('adjoint_pressure_RT_low_2e-8.dat', ' ', 0);
adjoint_RT_low_dt2       = matrix2cube(adjoint_RT_low_matrix, Nz);
adjoint_RT_low_dt2_slice = adjoint_RT_low_dt2(:, :, floor(Nz/2));
% Adjoint RT low - dt3
adjoint_RT_low_matrix    = importdata('adjoint_pressure_RT_low_1e-8.dat', ' ', 0);
adjoint_RT_low_dt3       = matrix2cube(adjoint_RT_low_matrix, Nz);
adjoint_RT_low_dt3_slice = adjoint_RT_low_dt3(:, :, floor(Nz/2));

%==================================================
% LOAD Adjoint - MID
%==================================================
% Adjoint mid
adjoint_kWave_mid       = adjoint_pressure_mid;
adjoint_kWave_mid_slice = adjoint_pressure_mid(:, :, floor(Nz/2));
normKWave_mid           = sqrt(sum(adjoint_kWave_mid(:).^2));
normKWave_mid_slice     = sqrt(sum(adjoint_kWave_mid_slice(:).^2));
normKWave_mid_max       = max(adjoint_kWave_mid(:));

% Adjoint RT mid - dt1
adjoint_RT_mid_matrix    = importdata('adjoint_pressure_RT_mid_4e-8.dat', ' ', 0);
adjoint_RT_mid_dt1       = matrix2cube(adjoint_RT_mid_matrix, Nz);
adjoint_RT_mid_dt1_slice = adjoint_RT_mid_dt1(:, :, floor(Nz/2));
% Adjoint RT mid - dt2
adjoint_RT_mid_matrix    = importdata('adjoint_pressure_RT_mid_2e-8.dat', ' ', 0);
adjoint_RT_mid_dt2       = matrix2cube(adjoint_RT_mid_matrix, Nz);
adjoint_RT_mid_dt2_slice = adjoint_RT_mid_dt2(:, :, floor(Nz/2));
% Adjoint RT mid - dt3
adjoint_RT_mid_matrix    = importdata('adjoint_pressure_RT_mid_1e-8.dat', ' ', 0);
adjoint_RT_mid_dt3       = matrix2cube(adjoint_RT_mid_matrix, Nz);
adjoint_RT_mid_dt3_slice = adjoint_RT_mid_dt3(:, :, floor(Nz/2));

%==================================================
% LOAD Adjoint - HIGH
%==================================================
% Adjoint high
adjoint_kWave_high       = adjoint_pressure_high;
adjoint_kWave_high_slice = adjoint_pressure_high(:, :, floor(Nz/2));
normKWave_high           = sqrt(sum(adjoint_kWave_high(:).^2));
normKWave_high_slice     = sqrt(sum(adjoint_kWave_high_slice(:).^2));
normKWave_high_max       = max(adjoint_kWave_high(:));

% Adjoint RT high - dt1
adjoint_RT_high_matrix    = importdata('adjoint_pressure_RT_high_4e-8.dat', ' ', 0);
adjoint_RT_high_dt1       = matrix2cube(adjoint_RT_high_matrix, Nz);
adjoint_RT_high_dt1_slice = adjoint_RT_high_dt1(:, :, floor(Nz/2));
% Adjoint RT high - dt2
adjoint_RT_high_matrix    = importdata('adjoint_pressure_RT_high_2e-8.dat', ' ', 0);
adjoint_RT_high_dt2       = matrix2cube(adjoint_RT_high_matrix, Nz);
adjoint_RT_high_dt2_slice = adjoint_RT_high_dt2(:, :, floor(Nz/2));
% Adjoint RT high - dt3
adjoint_RT_high_matrix    = importdata('adjoint_pressure_RT_high_1e-8.dat', ' ', 0);
adjoint_RT_high_dt3       = matrix2cube(adjoint_RT_high_matrix, Nz);
adjoint_RT_high_dt3_slice = adjoint_RT_high_dt3(:, :, floor(Nz/2));


%==================================================
% ERROR COMPUTATION
%==================================================
RT_norms = 0;
% NORMS
if(RT_norms)
normRT_low_dt1        = sqrt(sum(adjoint_RT_low_dt1(:).^2));
normRT_low_dt1_slice  = sqrt(sum(adjoint_RT_low_dt1_slice(:).^2));
normRT_low_dt2        = sqrt(sum(adjoint_RT_low_dt2(:).^2));
normRT_low_dt2_slice  = sqrt(sum(adjoint_RT_low_dt2_slice(:).^2));
normRT_low_dt3        = sqrt(sum(adjoint_RT_low_dt3(:).^2));
normRT_low_dt3_slice  = sqrt(sum(adjoint_RT_low_dt3_slice(:).^2));
normRT_mid_dt1        = sqrt(sum(adjoint_RT_mid_dt1(:).^2));
normRT_mid_dt1_slice  = sqrt(sum(adjoint_RT_mid_dt1_slice(:).^2));
normRT_mid_dt2        = sqrt(sum(adjoint_RT_mid_dt2(:).^2));
normRT_mid_dt2_slice  = sqrt(sum(adjoint_RT_mid_dt2_slice(:).^2));
normRT_mid_dt3        = sqrt(sum(adjoint_RT_mid_dt3(:).^2));
normRT_mid_dt3_slice  = sqrt(sum(adjoint_RT_mid_dt3_slice(:).^2));
normRT_high_dt1       = sqrt(sum(adjoint_RT_high_dt1(:).^2));
normRT_high_dt1_slice = sqrt(sum(adjoint_RT_high_dt1_slice(:).^2));
normRT_high_dt2       = sqrt(sum(adjoint_RT_high_dt2(:).^2));
normRT_high_dt2_slice = sqrt(sum(adjoint_RT_high_dt2_slice(:).^2));
normRT_high_dt3       = sqrt(sum(adjoint_RT_high_dt3(:).^2));
normRT_high_dt3_slice = sqrt(sum(adjoint_RT_high_dt3_slice(:).^2));
else
normRT_low_dt1        = normKWave_low;
normRT_low_dt1_slice  = normKWave_low_slice;
normRT_low_dt2        = normKWave_low;
normRT_low_dt2_slice  = normKWave_low_slice;
normRT_low_dt3        = normKWave_low;
normRT_low_dt3_slice  = normKWave_low_slice;
normRT_mid_dt1        = normKWave_mid;
normRT_mid_dt1_slice  = normKWave_mid_slice;
normRT_mid_dt2        = normKWave_mid;
normRT_mid_dt2_slice  = normKWave_mid_slice;
normRT_mid_dt3        = normKWave_mid;
normRT_mid_dt3_slice  = normKWave_mid_slice;
normRT_high_dt1       = normKWave_high;
normRT_high_dt1_slice = normKWave_high_slice;
normRT_high_dt2       = normKWave_high;
normRT_high_dt2_slice = normKWave_high_slice;
normRT_high_dt3       = normKWave_high;
normRT_high_dt3_slice = normKWave_high_slice;
end
% Error Computation - Low
ea_dt1_low   = adjoint_kWave_low/normKWave_low-adjoint_RT_low_dt1/normRT_low_dt1;
ea_dt1_low_s = sum(ea_dt1_low(:).^2)/sum((adjoint_kWave_low(:)/normKWave_low).^2);
ea_dt2_low   = adjoint_kWave_low/normKWave_low-adjoint_RT_low_dt2/normRT_low_dt2;
ea_dt2_low_s = sum(ea_dt2_low(:).^2)/sum((adjoint_kWave_low(:)/normKWave_low).^2);
ea_dt3_low   = adjoint_kWave_low/normKWave_low-adjoint_RT_low_dt3/normRT_low_dt3;
ea_dt3_low_s = sum(ea_dt3_low(:).^2)/sum((adjoint_kWave_low(:)/normKWave_low).^2);
errorConvergence_low = [ea_dt3_low_s ea_dt2_low_s ea_dt1_low_s];
% Error Computation - Mid
ea_dt1_mid   = adjoint_kWave_mid/normKWave_mid-adjoint_RT_mid_dt1/normRT_mid_dt1;
ea_dt1_mid_s = sum(ea_dt1_mid(:).^2)/sum((adjoint_kWave_mid(:)/normKWave_mid).^2);
ea_dt2_mid   = adjoint_kWave_mid/normKWave_mid-adjoint_RT_mid_dt2/normRT_mid_dt2;
ea_dt2_mid_s = sum(ea_dt2_mid(:).^2)/sum((adjoint_kWave_mid(:)/normKWave_mid).^2);
ea_dt3_mid   = adjoint_kWave_mid/normKWave_mid-adjoint_RT_mid_dt3/normRT_mid_dt3;
ea_dt3_mid_s = sum(ea_dt3_mid(:).^2)/sum((adjoint_kWave_mid(:)/normKWave_mid).^2);
errorConvergence_mid = [ea_dt3_mid_s ea_dt2_mid_s ea_dt1_mid_s];
% Error Computation - High
ea_dt1_high   = adjoint_kWave_high/normKWave_high-adjoint_RT_high_dt1/normRT_high_dt1;
ea_dt1_high_s = sum(ea_dt1_high(:).^2)/sum((adjoint_kWave_high(:)/normKWave_high).^2);
ea_dt2_high   = adjoint_kWave_high/normKWave_high-adjoint_RT_high_dt2/normRT_high_dt2;
ea_dt2_high_s = sum(ea_dt2_high(:).^2)/sum((adjoint_kWave_high(:)/normKWave_high).^2);
ea_dt3_high   = adjoint_kWave_high/normKWave_high-adjoint_RT_high_dt3/normRT_high_dt3;
ea_dt3_high_s = sum(ea_dt3_high(:).^2)/sum((adjoint_kWave_high(:)/normKWave_high).^2);
errorConvergence_high = [ea_dt3_high_s ea_dt2_high_s ea_dt1_high_s];



figure;
x_axis_conv = 1e9*[1e-8 2e-8 4e-8];
semilogy(x_axis_conv, errorConvergence_high, 'Color', 'b', 'LineWidth', 2); 
hold on;
semilogy(x_axis_conv, errorConvergence_mid, 'Color', 'g', 'LineWidth', 2); 
semilogy(x_axis_conv, errorConvergence_low, 'Color', 'r', 'LineWidth', 2); 
legend('Top ball', 'Middle ball', 'Bottom ball', 'location', 'southeast');
axis([0 50 1e-3 2])
xlabel('\Deltat [ns]')
ylabel('REE')
grid on;
set(gca,'FontSize',fontSize);
if(saveFigures)
saveas(gcf, 'Example13_ErrorConvergence', 'epsc');
saveas(gcf, 'Example13_ErrorConvergence.fig');
end

%==================================================
% PLOT SLICE - 3 Balls
%==================================================
positionYBar = [700 700 375 600];
positionY    = [700 700 315 600];
positionBar  = [700 700 350 600];
% kWave
figure;
imagesc(1e3*y_axis, 1e3*x_axis, (adjoint_kWave_low_slice+adjoint_kWave_mid_slice+adjoint_kWave_high_slice));
view(2);
box on;
axis tight;
xlabel('y [mm]');
ylabel('x [mm]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
pbaspect([1 2 1])
%colorbar();
if(saveFigures)
saveas(gcf, 'Example13_adjoint_kWave', 'png');
saveas(gcf, 'Example13_adjoint_kWave.fig');
end
% RT
figure;
imagesc(1e3*y_axis, 1e3*x_axis, (adjoint_RT_low_dt3_slice+adjoint_RT_mid_dt3_slice+adjoint_RT_high_dt3_slice));
view(2);
box on;
axis tight;
xlabel('y [mm]');
colorbar();
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
pbaspect([1 2 1])
caxis([-4.5e-6 9e-6])
if(saveFigures)
saveas(gcf, 'Example13_adjoint_RT', 'png');
saveas(gcf, 'Example13_adjoint_RT.fig');
end

% Difference
figure;
imagesc(1e3*y_axis, 1e3*x_axis, (adjoint_kWave_low_slice+adjoint_kWave_mid_slice+adjoint_kWave_high_slice)-(adjoint_RT_low_dt3_slice+adjoint_RT_mid_dt3_slice+adjoint_RT_high_dt3_slice));
view(2);
box on;
axis tight;
xlabel('y [mm]');
colorbar();
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
pbaspect([1 2 1])
caxis([-4.5e-6 9e-6])
if(saveFigures)
saveas(gcf, 'Example13_adjoint_error', 'png');
saveas(gcf, 'Example13_adjoint_error.fig');
end

%%  % Adjoint
%%  figure;
%%  imagesc(adjoint_RT_low_dt3_slice/normRT_low_dt3_slice);
%%  colorbar();
%%  pbaspect([1 2 1])
%%  figure;
%%  imagesc(adjoint_RT_mid_dt3_slice/normRT_mid_dt3_slice);
%%  colorbar();
%%  pbaspect([1 2 1])
%%  figure;
%%  imagesc(adjoint_RT_high_dt3_slice/normRT_high_dt3_slice);
%%  colorbar();
%%  pbaspect([1 2 1])
%%  % kWave
%%  figure;
%%  imagesc(adjoint_kWave_low_slice/normKWave_low_slice);
%%  colorbar();
%%  pbaspect([1 2 1])
%%  figure;
%%  imagesc(adjoint_kWave_mid_slice/normKWave_mid_slice);
%%  colorbar();
%%  pbaspect([1 2 1])
%%  figure;
%%  imagesc(adjoint_kWave_high_slice/normKWave_high_slice);
%%  colorbar();
%%  pbaspect([1 2 1])


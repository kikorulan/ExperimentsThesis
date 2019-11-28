% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex01_caustics;

clear all;
close all;


% Select caustics to compute
caustic1 = 1;
caustic2 = 1;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 384;     % number of grid points in the x (row) direction
Ny = 128;     % number of grid points in the y (column) direction
Nz = 128;     % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
dz = 1e-4;        % grid point spacing in the y direction [m]
dimensions = [Nx Ny Nz; dx dy dz];
dlmwrite('dimensions.dat', dimensions, 'delimiter', ' ');

%=================================================================================================================
%=================================================================================================================
%==================                        =======================================================================
%==================    CAUSTIC TYPE 1      =======================================================================
%==================                        =======================================================================
%=================================================================================================================
%=================================================================================================================
if(caustic1)
%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
v1 = -.3;
radi = Ny/5;
kernelSize = 10;
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
dimZ = Nz + 2*kernelSize;
c = c0*ones(dimX, dimY, dimZ);
c = addSphere(c, floor(dimX/2),   floor(dimY/2), floor(dimZ/2), radi, c0*v1);

cSlice = c(:, :, floor(dimZ/2));
K = ones(kernelSize, kernelSize, kernelSize);
cConvSlice = convn(cSlice, K, 'same')/kernelSize/kernelSize/kernelSize;
cConvSlice = cConvSlice(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
cConv = convn(c, K, 'same')/kernelSize/kernelSize/kernelSize;
cConv = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);

% Build initial pressure
%c_matrix = cube2matrix(repmat(cConvSlice, [1 1 Nz]));
c_matrix = cube2matrix(cConv);
dlmwrite('sound_speed_caustic1.dat', c_matrix, 'delimiter', ' ');
dlmwrite('initial_pressure.dat', c_matrix, 'delimiter', ' ');


% SIMULATION
system('./RTsolver dimensions.dat sound_speed_caustic1.dat initial_pressure.dat sensor_caustic1.dat');

%=======================
% TRAJECTORIES
%=======================
% Import data
filenameData = 'output_data/Trajectory0_caustic1.dat';
trajectories = importdata(filenameData, ' ', 0);

% Read number of rays and steps
[nSteps nRays] = size(trajectories);
%nRays = floor(nRays/2);
xCoord = trajectories(:, 1:6:nRays);
yCoord = trajectories(:, 2:6:nRays);
zCoord = trajectories(:, 3:6:nRays);

nRays = floor(nRays/6);
% Plot the figure
figure;
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Nx*dx 0 Ny*dy 0 Nz*dz]);
hold on;
colours = winter(nRays);
for n = 1:3:nRays
    plot3(xCoord(:, n), yCoord(:, n), zCoord(:, n), 'Color', colours(n, :));
end
pbaspect([3 1 1]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
box on;

% Sphere
[X, Y, Z] = sphere;
xS = (radi*dx*X + Nx*dx/2);
yS = (radi*dx*Y + Ny*dy/2);
zS = (radi*dx*Z + Nz*dz/2);
surf(xS, yS, zS, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'c');

% Caustic
[X, Y, Z] = sphere;
xC = (5*dx*X + 0.023);
yC = (5*dx*Y + (Ny-1)*dy/2);
zC = (5*dx*Z + (Nz-1)*dz/2);
surf(xC, yC, zC, 'EdgeColor', 'none', 'FaceColor', 'r');
saveas(gcf, 'Example01_3D_caustic1_rays', 'png');
saveas(gcf, 'Example01_3D_caustic1_rays.fig');

%=======================
% SOUND SPEED
%=======================
figure;
imagesc(permute(cConv(floor(Nx/2), :, :), [2 3 1]));
title('X slice');
xlabel('z axis');
ylabel('y axis');
colorbar();
pbaspect([1 1 1]);
saveas(gcf, 'Example01_3D_caustic1_xSlice', 'png');
saveas(gcf, 'Example01_3D_caustic1_xSlice.fig');

figure;
imagesc(permute(cConv(:, floor(Ny/2), :), [3 1 2]));
title('Y slice');
xlabel('x axis');
ylabel('z axis');
colorbar();
pbaspect([3 1 1]);
saveas(gcf, 'Example01_3D_caustic1_ySlice', 'png');
saveas(gcf, 'Example01_3D_caustic1_ySlice.fig');

figure;
imagesc(cConv(:, :, floor(Nz/2))');
title('Z slice');
xlabel('x axis');
ylabel('y axis');
colorbar();
pbaspect([3 1 1]);
saveas(gcf, 'Example01_3D_caustic1_zSlice', 'png');
saveas(gcf, 'Example01_3D_caustic1_zSlice.fig');
end

%=================================================================================================================
%=================================================================================================================
%==================                        =======================================================================
%==================    CAUSTIC TYPE 2      =======================================================================
%==================                        =======================================================================
%=================================================================================================================
%=================================================================================================================
if(caustic2)
%==============================
% define the properties of the propagation medium
%==============================
c0 = 1500;
v1 = -0.3;
kernelSize = 20;
% Auxiliary matrices
dimY = Nx + 2*kernelSize;
dimX = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
vx = kernelSize + floor(Ny/2);
vy = floor(Nx/4);
focus = floor(16*Nx/60);
c = addParabola(c, vx, vy, focus, c0*v1);
%c(:, Ny/10:end) = c(:, Ny/10:end)*(1 + v1);
% Kernel convolution
K = ones(kernelSize);
%cConv = c;
cConv_slice = conv2(c, K, 'same')/kernelSize/kernelSize;
cConv_slice = cConv_slice(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
cConv_slice = permute(cConv_slice, [2 1]);
cConv = repmat(cConv_slice, [1 1 Nz]);
% Write Cube
c_matrix = cube2matrix(cConv);
dlmwrite('sound_speed_caustic2.dat', c_matrix, 'delimiter', ' ');


% SIMULATION
system('./RTsolver dimensions.dat sound_speed_caustic2.dat initial_pressure.dat sensor_caustic2.dat');

%=======================
% TRAJECTORIES
%=======================
% Import data
filenameData = 'output_data/Trajectory0_caustic2.dat';
trajectories = importdata(filenameData, ' ', 0);

% Read number of rays and steps
[nSteps nRays] = size(trajectories);
%nRays = floor(nRays/2);
xCoord = trajectories(:, 1:6:nRays);
yCoord = trajectories(:, 2:6:nRays);
zCoord = trajectories(:, 3:6:nRays);

nRays = floor(nRays/6);
% Plot the figure
figure;
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Nx*dx 0 Ny*dy 0 Nz*dz]);
hold on;
colours = winter(nRays);
for n = 1:2:nRays
    plot3(xCoord(:, n), yCoord(:, n), zCoord(:, n), 'Color', colours(n, :));
end
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
pbaspect([3 1 1]);
box on;

% Hyperbolic paraboloid
a = 1/4/(focus-vy)/dx;
y = 0:dy:Ny*dx;
x = a*(y-(vx-kernelSize)*dx).^2 + vy*dy;
z = y;
z(:) = Nz*dz/6;
xPar = [x; x];
yPar = [y; y];
zPar = [z; 4*z];
surf(xPar, yPar, zPar, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'c');

% Caustic
zC = 0.004:1e-4:0.0085;
xC = zC;
xC(:) = 0.0152;
yC = zC;
yC(:) = 0.0062;
plot3(xC, yC, zC, 'Color', 'r', 'LineWidth', 5, 'MarkerSize', 12);
saveas(gcf, 'Example01_3D_caustic2_rays', 'png');
saveas(gcf, 'Example01_3D_caustic2_rays.fig');

%=======================
% SOUND SPEED
%=======================
figure;
imagesc(permute(cConv(floor(Nx/2), :, :), [2 3 1]));
title('X slice');
xlabel('z axis');
ylabel('y axis');
colorbar();
pbaspect([1 1 1]);
saveas(gcf, 'Example01_3D_caustic2_xSlice', 'png');
saveas(gcf, 'Example01_3D_caustic2_xSlice.fig');

figure;
imagesc(permute(cConv(:, floor(Ny/2), :), [3 1 2]));
title('Y slice');
xlabel('x axis');
ylabel('z axis');
colorbar();
pbaspect([3 1 1]);
saveas(gcf, 'Example01_3D_caustic2_ySlice', 'png');
saveas(gcf, 'Example01_3D_caustic2_ySlice.fig');

figure;
imagesc(cConv(:, :, floor(Nz/2))');
title('Z slice');
xlabel('x axis');
ylabel('y axis');
colorbar();
pbaspect([3 1 1]);
saveas(gcf, 'Example01_3D_caustic2_zSlice', 'png');
saveas(gcf, 'Example01_3D_caustic2_zSlice.fig');
end


%=================================================================================================================
%=================================================================================================================
%==================                        =======================================================================
%==================    CAUSTIC TYPE 3      =======================================================================
%==================                        =======================================================================
%=================================================================================================================
%=================================================================================================================

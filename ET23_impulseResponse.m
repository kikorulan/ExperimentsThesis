%================================================================================
% ET23
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex23_impulseResponse;

clear all;
close all;

load sensor_data;
%==============================
% Parameters
%==============================
% Save results
saveResults = 1;
run_kWave   = 0;


% Plot parameters
position    = [500 500 300 400];
positionHor = [500 500 600 300];
axisNA = [-3e-5 3e-5];
% Parameters
factor = 50;
offset = 8e-9; % -2e-8
dt = 1.6e-7; % 1e-7
Npoints = 5;
t = -Npoints*dt:dt/factor:Npoints*dt;
epsilon = 1/dt*10; %1/dt/10
c = 1500;

%================================================================================
%=================== KWAVE IMPULSE RESPONSE  ====================================
%================================================================================
if(run_kWave)
% Create the computational grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
% Create the time array
xMax = 2*Nx*dx;
tMax = 2*xMax/c;

Filter = 1;
t_array = 1;

%==============================
% k-Wave simulation
%==============================
medium.sound_speed = c;
kgrid = makeGrid(Nx, dx, Ny, dy);
kgrid.t_array = 0:dt:tMax;
% Define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1;
% Define the acoustic parameters to record
sensor.record = {'p'};

% define a single source point
source.p0 = zeros(Nx, Ny);
source.p0(end - Nx/4, Ny/2) = 1;
% Input arguments
input_args = {'PMLInside', false};
% Run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save sensor_data sensor_data kgrid;
end
%================================================================================
%=================== NUMERICAL APPROXIMATION ====================================
%================================================================================
factorTime = 1e6*50;
%==============================
% Delta Function
%==============================
deltaG = @(x, a) a/sqrt(pi)*exp(-x.*x.*a^2);
deltaGP = @(x, a) -2*x.*a^3/sqrt(pi).*exp(-x.*x.*a^2);
deltaGPvector = deltaGP(t, epsilon);
% Plot of derivative of delta
figure;
plot(factorTime*t, deltaGPvector, 'LineWidth', 1.5);
grid on;
%box on;
xlabel('t [\mus]');
ylabel('amplitude');
xlim([-20 20]);
set(gcf, 'pos', position);
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET23_Impulse_Delta', 'epsc');
end
%==============================
% Greens Function 
%==============================
% Step function
H = @(x) 0.5*sign(x) + 0.5; 
G = @(c, t, offset) H(c*(t))./(c*(t+offset));
Gvector = G(c, t, offset);
Gvector(isinf(Gvector)) = 0;

% Plot G
figure;
plot(factorTime*t, Gvector, 'LineWidth', 1.5);
grid on;
%box on;
xlabel('t [\mus]');
ylabel('amplitude');
xlim([-20 20]);
set(gcf, 'pos', position);
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET23_Impulse_Green', 'epsc');
end
% Plot G derivative
Gt = conv(deltaGPvector, Gvector, 'same');
figure;
plot(factorTime*t, Gt, 'LineWidth', 1.5);
grid on;
%box on;
xlabel('t [\mus]');
ylabel('amplitude');
xlim([-20 20]);
set(gcf, 'pos', position);
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET23_Impulse_DeltaGreen', 'epsc');
end

maxG = max(Gt)
minG = min(Gt)
ratio = max(Gt)/min(Gt)

%================================================================================
%=================== COMPARISON BETWEEN THE TWO =================================
%================================================================================
xMax = 2*kgrid.Nx*kgrid.dx;
tMax = 2*xMax/c;
factorTime = 1e6;
% Normalise signals
Gt_NA = Gt/max(Gt);
Gt_KW = spline(kgrid.t_array, sensor_data.p, 0:dt:tMax);
Gt_KW = Gt_KW/max(Gt_KW);
Npoints_subset = Npoints*factor;%Npoints/2;

% Build comparable signals: Numerical approximation
delayNA = find(Gt_NA == 1)-1;
delayKW = find(Gt_KW == 1);
Gt_NA_subset = Gt_NA(delayNA-Npoints_subset/2:delayNA+Npoints_subset);
Gt_KW_subset = Gt_KW(delayKW-Npoints_subset/2:delayKW+Npoints_subset);
time = (-Npoints_subset/2)*dt:dt:Npoints_subset*dt;

% Plot signals
figure;
hold on;
plot(factorTime*time, Gt_NA_subset, 'r', 'LineWidth', 2);
plot(factorTime*time, Gt_KW_subset, 'b');
grid on; box on;
xlabel('t [\mus]');
ylabel('amplitude');
axis([-20 30 -.6 1.2]);
legend('Numerical approx.', 'k-Wave');
set(gcf, 'position', positionHor);
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET23_Impulse_kWaveComparison', 'epsc');
end
% Error
errorGt = Gt_NA_subset - Gt_KW_subset;
figure;
plot(factorTime*time, errorGt, 'g', 'LineWidth', 2);
grid on;
xlabel('t [\mus]');
ylabel('amplitude');
axis([-20 30 -.15 .15]);
set(gcf, 'position', positionHor);
set(gca, 'FontSize', 15);
if(saveResults)
saveas(gcf, 'ET23_Impulse_error', 'epsc');
end

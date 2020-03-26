%================================================================================
% ET21
% Single GB to obtain error at the boundary of domain and compare with kWave
% Inspired in Formula (35) of "Fast Multiscale Gaussian Wavepacket Transforms" by Qian and Ying
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex21_singleGB;

close all;
clear all;

load sensor_data_singleGB_kWave;

saveResults   = 0;
approximateGB = 1;
runSingleGB   = 1;
runMultipleGB = 0;
%====================================================================================================
%=======                     ========================================================================
%=======    DOMAIN           ========================================================================
%=======                     ========================================================================
%====================================================================================================

% Define Domain
N = 128;
dx = 1e-4;
x_axis = 0:dx:(N-1)*dx;
y_axis = x_axis;
[Y, X] = meshgrid(x_axis, y_axis);

%====================================================================================================
%=======                                 ============================================================
%=======    GAUSSIAN BEAM (PHI APPROX)   ============================================================
%=======                                 ============================================================
%====================================================================================================
% Define GB function
phi_var = @(f, Ll, sl, index_x, index_y, X, Y) ...
         (sqrt(pi/N/Ll)*sl)^2 * ...
         exp(2*pi*i*((X - index_x/Ll)*f(1) + (Y - index_y/Ll)*f(2))) .* ...
         exp(-sl*sl*pi*pi*((X - index_x/Ll).*(X - index_x/Ll) + (Y - index_y/Ll).*(Y - index_y/Ll)));

% Parametrize gaussian beam function
f0 = 1.5e3;
f = [0; f0];
Ll = N*pi;
index_x = N*N*pi/2*dx;
index_y = N*N*pi/2*dx;
sl = N;

% Plot Gaussian beam
GB_phi = phi_var(f, Ll, sl, index_x, index_y, X, Y);
GB_phi = GB_phi/max(real(GB_phi(:)));
singleGB = GB_phi;
save singleGB singleGB;

%====================================================================================================
%=======                                 ============================================================
%=======    GAUSSIAN BEAM (RAY TRACING)  ============================================================
%=======                                 ============================================================
%====================================================================================================
if(approximateGB)
% Define functions
amplitude_D = @(alpha, gamma, c0, t) 0.5./sqrt(1 + i*alpha*(gamma-t*c0*c0)/c0);
phi_D = @(alpha, gamma, delta, c0, t) 1/c0/c0*(c0*(gamma - t*c0*c0) + i*alpha*((gamma - t*c0*c0).*(gamma - t*c0*c0) + delta.*delta./(1 + i*alpha*t*c0)));
GB_D = @(alpha, gamma, delta, c0, w, t) amplitude_D(alpha, gamma, c0, t).*exp(i*2*pi*w*phi_D(alpha, gamma, delta, c0, t));

% Parametrize
c0 = 1500;
x0 = [index_x/Ll; index_y/Ll];
w =  f0;
p0 = f/w/c0;
normP0 = norm(p0);
sigma = sl;
matrix = [X(:)'; Y(:)'];
coord = GB.findCoordinates(matrix, x0, p0);
gamma = reshape(coord(1, :), [N, N]);
delta = reshape(coord(2, :), [N, N]);
alpha = pi/2*sigma*sigma/w;

% Plot Gaussian beam
t = 0;%-1e-6;
GB_r = GB_D(alpha, gamma, delta, c0, w, t);
figure;
surf(x_axis, y_axis, real(2*GB_r), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();


% Plot Gaussian beam
t = 0;
GB_r = GB_D(alpha, gamma, delta, c0, w, t);
figure;
surf(x_axis, y_axis, real(GB_phi-2*GB_r), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();

% Generate difference as function of alpha
t = 0;
alpha_ini = 1;
alpha_end = 2;
alpha_n   = 20;
vec = (alpha_ini:(alpha_end-alpha_ini)/(alpha_n-1):alpha_end);
alpha   = sigma*sigma/w*vec;
dif_vec = zeros(size(alpha));
for ii = 1:alpha_n
    GB_a = GB_D(alpha(ii), gamma, delta, c0, w, t);
    val = sum((real(GB_phi(:)) - real(2*GB_a(:))).^2);
    dif_vec(ii) = val;
end
figure;
plot(vec, dif_vec);

singleGB_approx = GB_r;
save singleGB_approx singleGB_approx;
end
%====================================================================================================
%=======                                      =======================================================
%=======    PRESSURE TIME SERIES (SINGLE GB)  =======================================================
%=======                                      =======================================================
%====================================================================================================
if(runSingleGB)
% Sensor data
xS = [x_axis; 0*x_axis];
dt    = 2e-8;
tSync = 1e-8;
tMax  = 1e-5;
t_array = tSync:dt:(tMax+tSync);
GBF = zeros(N, length(t_array));

% Parametrize
c0 = 1500;
x0 = [index_x/Ll; index_y/Ll];
w =  f0;
p0 = f/w/c0;
normP0 = norm(p0);
sigma = sl;
alpha = pi/2*sigma*sigma/w;

% Loop over 
for ii = 1:N
    GBF(ii, :) = GB.gb(xS(:, ii), x0, p0, -t_array, alpha, w);
end

% Plot GBF
figure;
surf(t_array, x_axis, real(GBF), 'EdgeColor', 'none');
view(2);
axis tight

% Plot differences with kWave
figure;
surf(t_array, x_axis, real(sensor_data), 'EdgeColor', 'none');
view(2);
axis tight

figure;
surf(t_array, x_axis, real(sensor_data-GBF), 'EdgeColor', 'none');
view(2);
axis tight

sensor_data_singleGB_approx = GBF;
save sensor_data_singleGB_approx t_array sensor_data_singleGB_approx;
%%  % Plot 5 sensors
%%  for ii = 0:9
%%      n_sensor = 64 + 5*ii;
%%      figure;
%%      hold on;
%%      plot(real(GBF(n_sensor, :)));
%%      plot(sensor_data(n_sensor, :));
%%      title(['Sensor ', int2str(n_sensor)]);
%%      legend('GB', 'KW');
%%  end

%==================================================
%=======    EVALUATE ALPHA               ==========
%==================================================
%%  normGB = max(real(GBF(:)));
%%  normKW = max(real(sensor_data(:))); 
%%  
%%  alpha_ini = 1;
%%  alpha_end = 2;
%%  alpha_n   = 20;
%%  vec = (alpha_ini:(alpha_end-alpha_ini)/(alpha_n-1):alpha_end);
%%  alpha = sigma*sigma/w*vec;
%%  GBF_alpha = zeros(N, length(t_array), alpha_n);
%%  % Loop 
%%  for jj = 1:alpha_n 
%%      for ii = 1:N
%%          GBF_alpha(ii, :, jj) = gb(xS(:, ii), x0, p0, -t_array, alpha(jj), w);
%%      end
%%  end
%%  % Plot error
%%  error_vec = zeros(size(vec));
%%  for jj = 1:alpha_n 
%%      GBF_j = GBF_alpha(:, :, jj);
%%      error_alpha(jj) = sum(real(sensor_data(:) - GBF_j(:)).^2);
%%  end
%%  figure;
%%  plot(vec, error_alpha);
%%  title('Error as function of alpha');

%==================================================
%=======    EVALUATE SYNCHRONIZATION     ==========
%==================================================
%%  load sensor_data_singleGB;
%%  normGB = max(real(GBF(:)));
%%  normKW = max(real(sensor_data(:))); 
%%  
%%  
%%  t_ini = -2e-8;
%%  t_end = 2e-8;
%%  t_n   = 50;
%%  t_vec = (t_ini:(t_end-t_ini)/(t_n-1):t_end);
%%  GBF_t = zeros(N, length(t_array), t_n);
%%  alpha = pi/2*sigma*sigma/w;
%%  % Loop 
%%  for jj = 1:t_n
%%      for ii = 1:N
%%          t_array_aux = t_vec(jj) + (0:dt:tMax);
%%          GBF_t(ii, :, jj) = gb(xS(:, ii), x0, p0, -t_array_aux, alpha, w);
%%      end
%%  end
%%  % Plot error
%%  error_vec = zeros(size(vec));
%%  for jj = 1:t_n 
%%      GBF_j = GBF_t(:, :, jj);
%%      error_t(jj) = sum(real(sensor_data(:) - GBF_j(:)).^2);
%%  end
%%  figure;
%%  plot(t_vec, error_t);
%%  title('Error as function of t sync');
end

%====================================================================================================
%=======                                        =====================================================
%=======    PRESSURE TIME SERIES (MULTIPLE GB)  =====================================================
%=======                                        =====================================================
%====================================================================================================
if(runMultipleGB)
%========================================
% FFT (MANUAL)
%========================================
Nf = N;
fU0_manual = zeros(Nf, Nf);
x_vec = [Y(:) X(:)];
f_vec = zeros(Nf*Nf, 2);
fWidth = 5e3;
for jj = 1:Nf
    disp(jj)
    for ii = 1:Nf
        f_vec(ii + (jj-1)*Nf, :) = [-fWidth+(ii-1)*2*fWidth/Nf,  -fWidth+(jj-1)*2*fWidth/Nf];
        fU0_manual(ii, jj) = 1/N*sum(exp(-2*i*pi*sum(bsxfun(@times, x_vec, f_vec(ii + (jj-1)*Nf, :)), 2)).*singleGB(:));
    end
end
f_axis = f_vec(1:Nf, 1);
figure;
surf(f_axis, f_axis, real(fU0_manual), 'EdgeColor', 'none');
title('Real fft u0');
colorbar();
view(2);
figure;
surf(f_axis, f_axis, imag(fU0_manual), 'EdgeColor', 'none');
title('Imag fft u0');
colorbar();
view(2);

%========================================
% OBTAIN C COEFFICIENTS
%========================================
[valMax, indexM] = max(abs(fU0_manual(:)));
[indexC, indexR] = ind2sub([Nf, Nf], indexM);
threshold = 1e-2;
fU0_manual(abs(fU0_manual) < threshold) = 0;

% Create GB
Wl_1 = 2*fWidth/Nf;
Ll_1 = 2*Wl_1*2;
sl_1 = Wl_1;

% N elements
Ngb = sum(sum(abs(fU0_manual) > 0));
Nc = 5;
phi_GB = zeros(N, N);
c_vec = zeros(Nc, Nc, Ngb);
fU0_manualC = fU0_manual;
for kk = 1:Ngb
    [valMax, indexM] = max(abs(fU0_manualC(:)));
    % Obtain C coefficients
    epsilon = f_vec(indexM, :);
    for ii = 1:Nc
        for jj = 1:Nc
            c_vec(ii, jj, kk) = 1/Ll_1*exp(2*pi*i*((ii-1)*epsilon(2) + (jj-1)*epsilon(1))/Ll_1)*fU0_manual(indexM);
        end
    end
    % Obtain GB
    for ii = 1:Nc
        for jj = 1:Nc
            phi_GB = phi_GB + c_vec(ii, jj, kk)*(sqrt(pi/N/Ll_1)*sl_1).^2*exp(2*pi*i*((X-(ii-1)/Ll_1)*f_vec(indexM, 2) + (Y-(jj-1)/Ll_1)*f_vec(indexM, 1))).*exp(-sl_1^2*pi*pi*((X-(ii-1)/Ll_1).^2 + (Y-(jj-1)/Ll_1).^2));
        end
    end
    % Eliminate coefficient
    fU0_manualC(indexM) = 0;
end
figure;
surf(x_axis, y_axis, real(phi_GB), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
if(saveResults)
saveas(gcf, 'ET21_MultGB', 'epsc')
end

figure;
surf(x_axis, y_axis, real(phi_GB) - real(singleGB), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
if(saveResults)
saveas(gcf, 'ET21_MultGB_error', 'epsc')
end

%========================================
% COMPUTE GAUSSIAN BEAM
%========================================
% Sensor data
xS = [x_axis; 0*x_axis];
dt    = 2e-8;
tSync = 1e-8;
tMax  = 1e-5;
t_array = tSync:dt:(tMax+tSync);
GBF = zeros(N, length(t_array));

c0 = 1500;
% Loop over GB
for mm = 1:N
    fU0_manualC = fU0_manual;
    for kk = 1:Ngb
        [valMax, indexM] = max(abs(fU0_manualC(:)));
        % Obtain C coefficients
        epsilon = f_vec(indexM, :);
        for ii = 1:Nc
            for jj = 1:Nc
                % Parametrize
                x0 = [(ii-1)/Ll_1; (jj-1)/Ll_1];
                w =  norm(epsilon);
                p0 = fliplr(epsilon)'/w/c0;
                normP0 = norm(p0);
                sigma = sl_1;
                alpha = pi/2*sigma*sigma/w;
                GBF(mm, :) = GBF(mm, :) + c_vec(ii, jj, kk)*(GB.gb(xS(:, mm), x0, p0, -t_array, alpha, w) + GB.gb(xS(:, mm), x0, p0, t_array, alpha, w));
            end
        end
        % Eliminate coefficient
        fU0_manualC(indexM) = 0;
    end
end

% Norm
normGB = max(real(GBF(:)));
normKW = max(real(sensor_data(:)));
% Plot GBF
figure;
surf(t_array, x_axis, real(GBF)/2, 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
if(saveResults)
saveas(gcf, 'ET21_MultiGB_sinogram', 'epsc');
end

% Plot differences with kWave
figure;
surf(t_array, x_axis, real(sensor_data), 'EdgeColor', 'none');
view(2);
axis tight

figure;
surf(t_array, x_axis, real(sensor_data-GBF/2), 'EdgeColor', 'none');
view(2);
axis tight
colorbar();
if(saveResults)
saveas(gcf, 'ET21_MultiGB_errorSinogram', 'epsc');
end

% Plot 5 sensors
for ii = 1:4
    n_sensor = 64 + 20*(ii-1);
    figure;
    hold on;
    plot(real(GBF(n_sensor, :))/normGB, 'LineWidth', 2);
    plot(sensor_data(n_sensor, :)/normKW, 'LineWidth', 1.5, 'Color', 'r');
    title(['Sensor ', int2str(n_sensor)]);
    legend('GB', 'KW');
    grid on;
    if (saveResults)
    saveas(gcf, ['ET21_MultiGB_sensor', int2str(n_sensor)], 'epsc');
    end
end

%==================================================
%=======    COMPUTE ERROR BY SENSOR      ==========
%==================================================
error_sensor = zeros(size(x_axis));
energy_sensor = zeros(size(x_axis));
for ii = 1:N
    error_sensor(ii) = sum(real(sensor_data(ii, :)/normKW - GBF(ii, :)/normGB).^2);
    energy_sensor(ii) = sum(real(GBF(ii, :)/normGB).^2);
end
% Plot error and energy
figure;
semilogy(x_axis, error_sensor, 'LineWidth', 2);
hold on;
semilogy(x_axis, error_sensor./max(energy_sensor, 1e-20), 'LineWidth', 2);
semilogy(x_axis, energy_sensor, 'LineWidth', 2);
legend('Error sensor', 'Error sensor (norm energy)', 'Energy sensor');
title('Error by sensor (multiple GB)');
grid on;
if(saveResults)
saveas(gcf, 'ET21_MultiGB_errorSensor', 'epsc')
end

end

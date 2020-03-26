%================================================================================
% ET22
% Obtain time solution from Gaussian Beams
% Formula (35) of "Fast Multiscale Gaussian Wavepacket Transforms" by Qian and Ying
% Run using MATLAB 2016 or newer
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex22_vesselsGB;
close all;
clear all;


load vessels_GBdecomposition; 
load sensor_data_kWave;

%==================================================
%=======                     ======================
%=======    DOMAIN           ======================
%=======                     ======================
%==================================================
% Define Domain
Nx = gb1.Nx;
dx = gb1.dx;
x_axis = gb1.x_axis;
y_axis = gb1.y_axis;
[Y, X] = meshgrid(x_axis, y_axis);
L = gb1.L;
bl = gb1.bl;
Ll = gb1.Ll;
sl = gb1.sl;
nBoxes = gb1.nBoxes;
box = gb1.box;
epsilon = gb1.epsilon;
cCoeff = gb1.cCoeff;
extra_c = gb1.extra_c;

%===============================================================================================
%=======                                     ===================================================
%=======    SUM OF GAUSSIAN BEAMS (SENSORS)  ===================================================
%=======                                     ===================================================
%===============================================================================================
 
% Sensor data
dt    = 2e-8;
tSync = 3e-8;
tMax  = 1e-5;
t_array = tSync:dt:(tMax+tSync);
GBF = zeros(Nx, length(t_array));
c0 = 1500;
% Loop over 
for ii = 1:L
    for index_x = 1:bl(ii)+extra_c
        disp([int2str(ii), int2str(index_x)])
        for index_y = 1:bl(ii)+extra_c
            for index_box_x = 1:bl(ii)
                for index_box_y = 1:bl(ii)
                    coeff = cCoeff{ii}(index_x + (bl(ii) + extra_c)*(index_y-1), index_box_x, index_box_y);
                    if (abs(coeff) > 1e-5)
                        for index_sensor = 1:Nx
                            x_sensor = [dx*(index_sensor-1); 0];
                            %x0 = [(index_x-ceil(extra_c/2))/Ll(ii); gb1.factorXY*(index_y-ceil(extra_c/2))/Ll(ii)];
                            %x0 = [(index_x-ceil(gb1.extra_c/2)); gb1.factorXY*(index_y-ceil(gb1.extra_c/2))]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
                            x0 = [(index_x-0.5); gb1.factorXY*(index_y-0.5)]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
                            w =  sqrt(sum(epsilon{ii}(index_box_x, index_box_y, :).^2));
                            p0 = permute(epsilon{ii}(index_box_x, index_box_y, :), [3 2 1])/w/c0;
                            normP0 = norm(p0);
                            sigma = sl(ii);
                            alpha = pi/2*sigma*sigma/w;
                            GBF(index_sensor, :) = GBF(index_sensor, :) + coeff*GB.gb(x_sensor, x0, p0, t_array, alpha, w);
                            GBF(index_sensor, :) = GBF(index_sensor, :) + coeff*GB.gb(x_sensor, x0, p0, -t_array, alpha, w);
                        end
                    end
                end
            end
        end
    end
end


%%  % Choose 3 Sensors
%%  sensor_index = [10, 30, 60];
%%  GBF3 = zeros(3, length(t_array));
%%  % Loop over 
%%  for ii = 1:L
%%      for index_x = 1:bl(ii)+extra_c
%%          disp([int2str(ii), int2str(index_x)])
%%          for index_y = 1:bl(ii)+extra_c
%%              for index_box_x = 1:bl(ii)
%%                  for index_box_y = 1:bl(ii)
%%                      coeff = cCoeff{ii}(index_x + (bl(ii) + extra_c)*(index_y-1), index_box_x, index_box_y);
%%                      if (abs(coeff) > 0)
%%                          x1 = [dx*(sensor_index(1)-1); 0];
%%                          x2 = [dx*(sensor_index(2)-1); 0];
%%                          x3 = [dx*(sensor_index(3)-1); 0];
%%                          x0 = [(index_x-ceil(extra_c/2))/Ll(ii); gb1.factorXY*(index_y-ceil(extra_c/2))/Ll(ii)];
%%                          w =  sqrt(sum(epsilon{ii}(index_box_x, index_box_y, :).^2));
%%                          p0 = permute(epsilon{ii}(index_box_x, index_box_y, :), [3 2 1])/w/c0;
%%                          normP0 = norm(p0);
%%                          sigma = sl(ii);
%%                          alpha = pi/2*sigma*sigma/w;
%%                          GBF3(1, :) = GBF(1, :) + coeff*GB.gb(x1, x0, p0, t_array, alpha, w);
%%                          GBF3(1, :) = GBF(1, :) + coeff*GB.gb(x1, x0, p0, t_array, alpha, w);
%%                          GBF3(2, :) = GBF(2, :) + coeff*GB.gb(x2, x0, p0, t_array, alpha, w);
%%                          GBF3(2, :) = GBF(2, :) + coeff*GB.gb(x2, x0, p0, t_array, alpha, w);
%%                          GBF3(3, :) = GBF(3, :) + coeff*GB.gb(x3, x0, p0, t_array, alpha, w);
%%                          GBF3(3, :) = GBF(3, :) + coeff*GB.gb(x3, x0, p0, -t_array, alpha, w);
%%                      end
%%                  end
%%              end
%%          end
%%      end
%%  end


%==================================================
% GAUSSIAN BEAM
%==================================================
pos = [0 0 1000 800];
normKW = max(sensor_data(:)); %sqrt(sum(sensor_data(:).^2));
normGB = max(real(GBF(:))); %sqrt(sum(real(GBF(:)).^2));
% Full set of sensors
figure;
imagesc(real(GBF))

figure;
imagesc(sensor_data);

% Difference
figure;
imagesc(real(GBF)/normGB - sensor_data/normKW);


%%  % Offset
%%  offset_pos = @(x, index) [x(:, index+1:end) zeros(size(x, 1), index)];
%%  for ii = 1:4
%%      figure;
%%      imagesc((offset_pos(real(GBF), ii))/normGB - sensor_data/normKW);
%%      colorbar()
%%  end
%%  
%%  figure;
%%  imagesc(real(GBF)/normGB - sensor_data/normKW);
%%  
%%  % Offset
%%  offset_neg = @(x, index) [zeros(size(x, 1), index) x(:, 1:end-index) ];
%%  for ii = 1:4
%%      figure;
%%      imagesc((offset_neg(real(GBF), ii))/normGB - sensor_data/normKW);
%%      colorbar()
%%  end


error_sensor = zeros(size(x_axis));
energy_sensor = zeros(size(x_axis));
for ii = 1:Nx
    error_sensor(ii) = sum(real(sensor_data(ii, :)/normKW - GBF(ii, :)/normGB).^2);
    energy_sensor(ii) = sum(real(sensor_data(ii, :)/normKW).^2);
end
% Plot error and energy
figure;
semilogy(x_axis, error_sensor, 'LineWidth', 2);
hold on;
semilogy(x_axis, error_sensor./max(energy_sensor, 1e-20), 'LineWidth', 2);
semilogy(x_axis, energy_sensor, 'LineWidth', 2);
legend('Error', 'Error/Energy', 'Energy');
title('Error by sensor (multiple GB)');
grid on;

%==================================================
% SAVE
%==================================================
sensor_data_GB = GBF;
t_array_GB = t_array;
save sensor_data_GB sensor_data_GB t_array_GB;

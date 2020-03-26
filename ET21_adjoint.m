%================================================================================
% ET21
% From the Gaussian Beam, compute Adjoint
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex21_singleGB;

close all;
clear all;

load adjoint_singleGB_kWave;
load sensor_data_singleGB_GBparameters;

adjoint_kWave = adjoint_singleGB_kWave.p_final;
%==================================================
%=======                     ======================
%=======    DOMAIN           ======================
%=======                     ======================
%==================================================
% Define Domain
Nx = 128;
Ny = 128;
dx = 1e-4;
dy = 1e-4;
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;
[Y, X] = meshgrid(x_axis, y_axis);
x_vec = [Y(:)'; X(:)';];
L = gb1.L;
bl = gb1.bl;
Ll = gb1.Ll;
sl = gb1.sl;
nBoxes = gb1.nBoxes;
box = gb1.box;
epsilon = gb1.epsilon;
cCoeff = gb1.cCoeff;
extra_c = gb1.extra_c;
t_len = size(gb1.u0, 2);

%===============================================================================================
%=======                                     ===================================================
%=======    SUM OF GAUSSIAN BEAMS (SENSORS)  ===================================================
%=======                                     ===================================================
%===============================================================================================
 
% Sensor data
dt = 2e-8;
GBF = zeros(1, Nx*Ny);
c0 = 1500;
% Loop over 
%x0 = [(index_x-ceil(extra_c/2))/Ll(ii); (index_y-ceil(extra_c/2))/Ll(ii)];

%%  % GB 1
%%  x0 = [(Nx-1)*dx*3/8; 0];
%%  Nt = 750*2/8;
%%  tp = Nt*dt;
%%  f = permute(epsilon{2}(5, 2, :), [3 2 1]);
%%  w =  sqrt(sum(f.^2));
%%  p0 = f/w/c0;
%%  sigma = sl(2);
%%  alpha = pi/2*sigma*sigma/w;
%%  GBF = GBF + GB.gb(x_vec, x0, p0, tp, alpha, w);
%%  GBF = GBF + GB.gb(x_vec, x0, p0, -tp, alpha, w);
%%  f = permute(epsilon{2}(4, 7, :), [3 2 1]);
%%  f(2) = f(2);
%%  w =  sqrt(sum(f.^2));
%%  p0 = f/w/c0;
%%  sigma = sl(2);
%%  alpha = pi/2*sigma*sigma/w;
%%  GBF = GBF + GB.gb(x_vec, x0, p0, tp, alpha, w);
%%  GBF = GBF + GB.gb(x_vec, x0, p0, -tp, alpha, w);
%%  
%%  % GB 2
%%  f = permute(epsilon{2}(4, 1, :), [3 2 1]);
%%  f(2) = f(2);
%%  w =  sqrt(sum(f.^2));
%%  p0 = f/w/c0;
%%  sigma = sl(2);
%%  alpha = pi/2*sigma*sigma/w;
%%  GBF = GBF + GB.gb(x_vec, x0, p0, tp, alpha, w);
%%  GBF = GBF + GB.gb(x_vec, x0, p0, -tp, alpha, w);


GBF = zeros(1, Nx*Ny);
% Loop over 
for ii = 1:L
    for index_x = 1:bl(ii)+extra_c
        disp([int2str(ii), int2str(index_x)])
        for index_y = 1:bl(ii)+extra_c
            for index_box_x = 1:bl(ii)
                for index_box_y = 1:bl(ii)
                    coeff = cCoeff{ii}(index_x + (bl(ii) + extra_c)*(index_y-1), index_box_x, index_box_y);
                    if (abs(coeff) > 0.001)
                        a = [index_x index_y bl(ii) index_box_x index_box_y coeff];
                        %x0 = [(Nx-1)*dx*(index_x-1)/(bl(ii)); 0];
                        %x0 = [(index_x-ceil(gb1.extra_c/2)); 0]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
                        %Nt = 750*(index_y-1)/bl(ii);
                        x0 = [(index_x-0.5); 0]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
                        Nt = (t_len)*(index_y-0.5)/(bl(ii)+extra_c);
                        tp = (Nt)*dt;
                        f = permute(epsilon{ii}(index_box_x, index_box_y, :), [3 2 1]);
                        w =  sqrt(sum(f.^2));
                        p0 = f/w/c0;
                        sigma = sl(ii);
                        alpha = pi/2*sigma*sigma/w;
                        GBF = GBF + coeff*GB.gb(x_vec, x0, p0, tp, alpha, w);
                        GBF = GBF + coeff*GB.gb(x_vec, x0, p0, -tp, alpha, w);
                    end
                end
            end
        end
    end
end



%%  alpha_ini = -4;
%%  alpha_end = 4;
%%  alpha_n   = 20;
%%  vec = (alpha_ini:(alpha_end-alpha_ini)/(alpha_n-1):alpha_end);
%%  error_vec = zeros(size(vec));
%%  for jj = 1:alpha_n
%%      disp(vec(jj))
%%      GBF = zeros(1, Nx*Ny);
%%      % Loop over 
%%      for ii = 1:L
%%          for index_x = 1:bl(ii)+extra_c
%%              disp([int2str(ii), int2str(index_x)])
%%              for index_y = 1:bl(ii)+extra_c
%%                  for index_box_x = 1:bl(ii)
%%                      for index_box_y = 1:bl(ii)
%%                          coeff = cCoeff{ii}(index_x + (bl(ii) + extra_c)*(index_y-1), index_box_x, index_box_y);
%%                          if (abs(coeff) > 0.001)
%%                              a = [index_x index_y bl(ii) index_box_x index_box_y coeff];
%%                              %x0 = [(Nx-1)*dx*(index_x-1)/(bl(ii)); 0];
%%                              %x0 = [(index_x-ceil(gb1.extra_c/2)); 0]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
%%                              %Nt = 750*(index_y-1)/bl(ii);
%%                              x0 = [(index_x-0.5); 0]/(gb1.bl(ii) + gb1.extra_c)*(gb1.Nx-1)*gb1.dx;
%%                              Nt = (t_len)*(index_y-.5)/(bl(ii)+extra_c);
%%                              tp = (Nt+vec(jj))*dt;
%%                              f = permute(epsilon{ii}(index_box_x, index_box_y, :), [3 2 1]);
%%                              w =  sqrt(sum(f.^2));
%%                              p0 = f/w/c0;
%%                              sigma = sl(ii);
%%                              alpha = pi/2*sigma*sigma/w;
%%                              GBF = GBF + coeff*GB.gb(x_vec, x0, p0, tp, alpha, w);
%%                              GBF = GBF + coeff*GB.gb(x_vec, x0, p0, -tp, alpha, w);
%%                          end
%%                      end
%%                  end
%%              end
%%          end
%%      end
%%      GBF_reshape = reshape(GBF, [Ny, Nx]);
%%      GBF_reshape = GBF_reshape';
%%      normGB = max(real(GBF_reshape(:)));
%%      normKW = max(real(adjoint_kWave(:)));
%%      error = real(adjoint_kWave)/normKW-real(GBF_reshape)/normGB;
%%      error_alpha(jj) = sum((error(:)).^2);
%%  end
%%  figure;
%%  plot(vec, error_alpha);

%=========================================================================
% PLOT
%=========================================================================
GBF_reshape = reshape(GBF, [Ny, Nx]);
GBF_reshape = GBF_reshape';

normGB = max(real(GBF_reshape(:)));
normKW = max(real(adjoint_kWave(:)));

figure;
imagesc(real(GBF_reshape));

figure;
imagesc(real(adjoint_kWave));

figure;
imagesc(real(adjoint_kWave)/normKW-real(GBF_reshape)/normGB);
colorbar();

%=========================================================================
% SAVE
%=========================================================================
adjoint_GB = real(GBF_reshape);
save adjoint_singleGB_GB adjoint_GB;

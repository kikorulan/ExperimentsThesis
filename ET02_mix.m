%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex02_simulation2D_homo;

close all;
%clear all;

load signalRT.mat
load adjoint_kWave.mat
load sensor_data_kWave.mat

%================================================================================
% kWave reconstruction using RT data
%================================================================================

sensor_data_RT = signalRT;
% run the time reversal reconstruction
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% Build sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p_mask(end, :) = 1;
source_adjoint.p_mask(:, 1) = 1;
source_adjoint.p_mask(:, end) = 1;
source_adjoint.p = fliplr(sensor_data_RT);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

%================================================================================
% RT reconstruction using kWave data
%================================================================================
%========================================
% Compute filters
%========================================
Rgrid.inverse_filter(100);

%========================================
% Compute reverse signal
%========================================
nSources = 764;
for n = 1:nSources
    source(n).setForwardSignal(sensor_data(n, :));
end
for n = 1:nSources
    disp(n)
    Rgrid.inverse_beam(source(n));
end
%Rgrid.computeAdjointParallel(source);

pixelAReverseSensors = Rgrid.inverse_signal(source);


adjointKWaveForward_RT = Rgrid.pixelAReverse;
adjointRTForward_kWave = p0_recon_adjoint;
save adjointKWaveForward_RT adjointKWaveForward_RT;
save adjointRTForward_kWave adjointRTForward_kWave;

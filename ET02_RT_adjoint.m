%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/Ex02_simulation2D_homo;

%close all;
%clear all;

load signalRT.mat
load adjoint_kWave.mat;
% Measure computational time
tic;
start_time = clock;
 
 
%========================================
% Compute filters
%========================================
Rgrid.inverse_filter(100);

%========================================
% Compute reverse signal
%========================================
%Rgrid.inverse_beam(1);

%aReverse = Rgrid.inverse_beam_adjoint();
nSources = 764;
for n = 1:nSources
    source(n).setForwardSignal(signalRT(n, :));
end
for n = 1:nSources
    disp(n)
    Rgrid.inverse_beam(source(n));
end
%Rgrid.computeAdjointParallel(source);

pixelAReverseSensors = Rgrid.inverse_signal(source);

adjoint_RT = Rgrid.pixelAReverse;
save adjoint_RT adjoint_RT;

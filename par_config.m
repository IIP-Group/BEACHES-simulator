% -----------------------------------------------------
% -- Simple MIMO simulator with estimated LoS channel
% -- 2018 (c) studer@cornell.edu, Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
function par = par_config(sim_scenario, include_ANM_NOMP)

fullpath = mfilename('fullpath');
simulator_path = strrep(fullpath, mfilename(), '');
addpath(genpath(simulator_path));
par.simulator_path = simulator_path;

if include_ANM_NOMP
    nomp = exist('extractSpectrum', 'file');
    anm = exist('ast_denoise', 'file');
    if ~(nomp == 2 && anm == 2)
        error('In order to simulate NOMP and/or ANM you must add their codes to the path. See README.md for more info.');
    end
end
% --------------------  parameters -----------------------------
switch(sim_scenario)
    case 'a'
        par.runId = 1; % simulation ID (used to reproduce results)    % 100: 256*16 Los - 101: 128*8 Los - 200: 256*16 nLoS - 201: 128*8 nLos - 300 128*8 Rayleigh
        par.B = 128; % number of BS antennas
        par.U = 8; % number of UEs
        par.channel = 'QuadMMLoS'; % 'Rayleigh', 'LoS', 'QuadMMLoS', 'QuadFreeSpace', 'QuadMMnLoS'
        par.p_fa = 0.1;         % false alarm rate for NOMP, optimized for each scenario by simulation a range of values of p_fa and picking the one that minimizes MSE.
        par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE'};
        par.SNRdB_list_L = {[-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15]}; % list of SNR [dB] values to be simulated
        par.n_channel_trials_L = {[100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000]};  % NOTE: overall number of trials per SNR point is (par.n_channel_trials(snr_index)*par.trials_per_channel*par.n_parallel_workers)
        if include_ANM_NOMP
            par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE', 'NOMP', 'ANM'};
            par.SNRdB_list_L = {[-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15]};
            par.n_channel_trials_L = {[100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000]};  
        end
    case 'b'
        par.runId = 2; 
        par.B = 256; 
        par.U = 16; 
        par.channel = 'QuadMMLoS'; 
        par.p_fa = 0.1;  
        par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE'};
        par.SNRdB_list_L = {[-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15]}; % list of SNR [dB] values to be simulated
        par.n_channel_trials_L = {[100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000]};  % NOTE: overall number of trials per SNR point is (par.n_channel_trials(snr_index)*par.trials_per_channel*par.n_parallel_workers)
        if include_ANM_NOMP
            par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE', 'NOMP', 'ANM'};
            par.SNRdB_list_L = {[-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15], [-10 -5 0 5 10 15]};
            par.n_channel_trials_L = {[100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 1000]};  
        end
        
    case 'c'
        par.runId = 3; 
        par.B = 128; 
        par.U = 8; 
        par.channel = 'QuadMMnLoS'; 
        par.p_fa = 0.6;
        par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE'};
        par.SNRdB_list_L = {[-10 -5 0 5 8], [-10 -5 0 5 10 12], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10]}; 
        par.n_channel_trials_L = {[100 100 100 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 ], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000]};  
        if include_ANM_NOMP
            par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE', 'NOMP', 'ANM'};
            par.SNRdB_list_L = {[-10 -5 0 5 8], [-10 -5 0 5 10 12], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10]}; 
            par.n_channel_trials_L = {[100 100 100 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 ], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000]};  
        end
        
    case 'd'
        par.runId = 4; 
        par.B = 256; 
        par.U = 16; 
        par.channel = 'QuadMMnLoS'; 
        par.p_fa = 0.8;
        par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE'};
        par.SNRdB_list_L = {[-10 -5 0 5 8], [-10 -5 0 5 10 12], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10]}; % list of SNR [dB] values to be simulated
        par.n_channel_trials_L = {[100 100 100 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 ], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000]};  % NOTE: overall number of trials per SNR point is (par.n_channel_trials(snr_index)*par.trials_per_channel*par.n_parallel_workers)
        if include_ANM_NOMP
            par.denoiser_list = {'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE', 'NOMP', 'ANM'};
            par.SNRdB_list_L = {[-10 -5 0 5 8], [-10 -5 0 5 10 12], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10], [-10 -5 0 5 10]}; 
            par.n_channel_trials_L = {[100 100 100 1000 1000], [100 100 100 1000 1000 1000], [100 100 100 1000 1000 ], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000], [100 100 100 1000 1000]};  
        end
        
    otherwise
        error('Simulation scenario not defined in par_config function.');
end
 
par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
% list of channel denoisers to be simulated. Choose from: % 'Perfect CSI', 'ML', 'BEACHES', 'BEACHES (hw)', 'BEACHES (fp)', 'exact MSE', 'ANM', 'NOMP'

if length(par.n_channel_trials_L) ~= length(par.SNRdB_list_L) || length(par.denoiser_list) ~= length(par.SNRdB_list_L)
    error('Number of elements in the cell arrays par.denoiser_list and par.SNRdB_list_L and par.n_channel_trials_L must be the same!')
end

for k = 1:length(par.n_channel_trials_L)
   if length(par.n_channel_trials_L{k}) ~=  length(par.SNRdB_list_L{k})
      error('Number of elements in each vector element of the cell arrays par.n_channel_trials_L and par.SNRdB_list_L must match!') 
   end
end


par.trials_per_channel = 20; % number of noise realizations within each channel realization (or one channel coherence interval)
par.n_parallel_workers = 4;  % Number of parallel workers, whose error rate results are fused to make better averaging. 
                             % enables use of matlab parfor which results in considerable simulation speedups, depending on the machine.
                             % It should not exceed 8! Please choose either 4 or 8 for performance or 1 when debugging.

% ######## Less needed to modify:
par.MinAngleSep = 1;    % Minumum angle separation between users in degrees in the LoS and nonLoS model 
par.dmin = 10;          % Minimum distance of users from the BS, based on mmMAGIC model
par.dmax = 110;         % Maximum distance of users from the BS, based on mmMAGIC model
par.ch_normalization = 'per-user';  % 'overall', 'per-user', 'no-normalization'
%par.n_channel_sets_available = 15;

par.MAX_MEM = 0.5e9;    % 0.5 Gbytes for load and store of mat files
par.N_channel_realizations = 1000;  % Number of channel realization per request from function 'generate_channels()'
if any(max(cell2mat(par.n_channel_trials_L)) > par.N_channel_realizations)
    error(['Requested number of channel realizations per SNR point should not exceed ' num2str(par.N_channel_realizations) '.\n'])
end

par.SAVE_CHANNEL = 1;
if strcmp(par.channel, 'Rayleigh')
    par.SAVE_CHANNEL = 0;
end


%% Initializations : DO NOT MODIFY
% -------------------------------------------------------------------
% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
        
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');


end
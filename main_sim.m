% -----------------------------------------------------
% -- main script that launches parallel error_rate_simulator tasks
% -- AUGUST 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
clearvars;
% Set the parameters in the par_config()
sim_scenario = 'b';
BER_figID = 1;      % figure ID for the BER plot
MSE_figID = 2;      % figure ID for the MSE plot
INCLUDE_ANM_NOMP = false;   % whether or not to include ANM and NOMP in simulations. 
                            % Please first download the ANM and NOMP matlab
                            % codes to ./denoisers directory to enable their
                            % simulation. See README.md for more info.
par = par_config(sim_scenario, INCLUDE_ANM_NOMP);
% merge subsets of channel realizations into one mat file containing 1000
% realizations. They are initially broken into small files to meet GitHub
% max file size limit.
merge_mat_files;
t0 = tic;
for den_idx = 1:length(par.denoiser_list)
    
    par.denoiser = par.denoiser_list{den_idx};
    par.SNRdB_list = par.SNRdB_list_L{den_idx};
    par.n_channel_trials = par.n_channel_trials_L{den_idx};
    
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('>> starting the simulation for %s channel estimation\n', par.denoiser);
    
    
    % parfor loop for multiple concurrent workers
    parfor worker_id = 1:par.n_parallel_workers
        res_worker(worker_id) = ber_mse_simulator(par, worker_id-1);
    end
    
    res.MSE = zeros(length(par.SNRdB_list),1); % vector error rate
    res.VER = zeros(length(par.SNRdB_list),1); % vector error rate
    res.SER = zeros(length(par.SNRdB_list),1); % symbol error rate
    res.BER = zeros(length(par.SNRdB_list),1); % bit error rate
    res.num_trials = zeros(length(par.SNRdB_list),1);
    for worker_id = 1:par.n_parallel_workers
        res.MSE = res.MSE + res_worker(worker_id).MSE;
        res.VER = res.VER + res_worker(worker_id).VER;
        res.SER = res.SER + res_worker(worker_id).SER;
        res.BER = res.BER + res_worker(worker_id).BER;
        res.num_trials = res.num_trials + res_worker(worker_id).num_trials;
    end
    sim_name{den_idx} = ['ERR-' num2str(par.runId), '_' num2str(par.B) 'x' num2str(par.U) '_' par.channel '_' par.mod '_' par.denoiser];
    legend_str{den_idx} = par.denoiser;
    save([par.simulator_path 'results/sim_res/' sim_name{den_idx}],'par','res');
    

    fprintf('>> finished the simulation for %s channel estimation in %g seconds. Current time: %s\n', par.denoiser, toc(t0)-t0, datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
    t0 = tic;
    
end


fprintf('************* Simulation finished at time: %s *************\n', datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'))

%% Plot BER and MSE results
sim_name2 = sim_name;
legend_str2 = legend_str;
sim_name2(~isempty(strfind(sim_name2, 'Perfect CSI'))) = [];
legend_str2(ismember(legend_str2,'Perfect CSI')) = [];

plot_BER(BER_figID, par.runId, sim_name, legend_str);


plot_MSE(MSE_figID, par.runId, sim_name2, legend_str2);
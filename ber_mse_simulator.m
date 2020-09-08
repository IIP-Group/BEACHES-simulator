% -----------------------------------------------------
% -- main simulator for 1-bit detection in massive MU-MIMO-OFDM systems
% -- MARCH 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------

function res = ber_mse_simulator(par, worker_id)

% use runId random seed (enables reproducibility)
rng(1000*par.runId + worker_id, 'twister');

% -- start simulation

% initialize result arrays 
res.MSE = zeros(length(par.SNRdB_list),1); % channel denoising MSE
res.VER = zeros(length(par.SNRdB_list),1); % vector error rate
res.SER = zeros(length(par.SNRdB_list),1); % symbol error rate
res.BER = zeros(length(par.SNRdB_list),1); % bit error rate
res.num_trials = par.n_channel_trials.'*par.trials_per_channel;
res.NOMP_K = zeros(length(par.SNRdB_list),2);

%% trials loop
fprintf(['Denoiser ' par.denoiser ' - Worker ' num2str(worker_id) ' started ...\n']);
H_load = generate_channels(par, worker_id);
for snr_idx = 1:length(par.SNRdB_list)
    for ch_trial = 1:par.n_channel_trials(snr_idx)
        if mod(ch_trial,10) == 1
            fprintf('Denoiser %s - Worker %d: SNR = %d dB, %d-th channel trial: - time: %s \n', par.denoiser, worker_id, par.SNRdB_list(snr_idx), ch_trial, datetime('now','TimeZone','local','Format','d-MMM HH:mm'));
        end
        H = H_load(:,:,ch_trial);
        N0 = (norm(H, 'fro')^2/par.B)*par.Es*10^(-par.SNRdB_list(snr_idx)/10);
        
        for trial = 1:par.trials_per_channel

            % generate transmit symbols
            bits = randi([0 1],par.U,par.Q);
            idx = bi2de(bits,'left-msb')+1;
            s = par.symbols(idx).';
            awgn_noise = (randn(par.B,1)+1i*randn(par.B,1));
            channel_noise = (randn(par.B,par.U)+1i*randn(par.B,par.U));
            
            switch par.denoiser
                case 'Perfect CSI'
                    Hest = H;
                case 'ML' 
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = Hnoisy;
                case 'BEACHES (hw)'
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = BEACHES_hw(par,Hnoisy, N0_est);
                case 'exact MSE'
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = exactMSE(par,Hnoisy, H);
                case 'BEACHES'
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = BEACHES(par,Hnoisy, H, N0_est, par.denoiser,0);
                case 'BEACHES (fp)'
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = BEACHES_fp(par,Hnoisy, H, N0_est);
                case {'ANM', 'AN_deb_N0', 'AN_deb_N0d2', 'AN_deb', 'AN_est'}
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    Hest = AN_denoiser(par, Hnoisy, N0_est, par.denoiser);
                case 'NOMP'
                    N0_est = N0/par.U/par.Es;
                    Hnoisy = H + channel_noise*sqrt(N0_est/2);
                    [Hest, num_freq] = NOMP_denoiser(par,Hnoisy, N0_est);
                    res.NOMP_K(snr_idx,1) = res.NOMP_K(snr_idx,1) + sum(num_freq)/nnz(num_freq);    % avg number of K
                    res.NOMP_K(snr_idx,2) = res.NOMP_K(snr_idx,1) + 1 - nnz(num_freq)/length(num_freq);    % avg number of failures in NOMP
                otherwise
                    error('par.denoiser not defiend')
            end
            
            n = sqrt(N0/2)*awgn_noise;
            y = H*s + n;
            
            %snr(snr_idx, trial_num) = (norm(H*s)/norm(n))^2;
            % Perform LMMSE detection with the esimated channel matrix

            [idxhat,bithat] = LMMSE(par,Hest,y,N0);

            
            err = (idx~=idxhat);
            res.VER(snr_idx) = res.VER(snr_idx) + sum(any(err));
            res.SER(snr_idx) = res.SER(snr_idx) + sum(sum(err))/par.U;
            res.BER(snr_idx) = res.BER(snr_idx) + sum(sum(bits ~= bithat))/(par.U*par.Q);
            res.MSE(snr_idx) = res.MSE(snr_idx) + mean(abs(H(:) - Hest(:)).^2);
            
        end
    end
end


end
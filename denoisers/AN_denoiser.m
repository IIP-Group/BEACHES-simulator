% -----------------------------------------------------
% -- Denoisers with Hard thresholding function
% -- 2018 (c) sm2675@cornell.edu, studer@cornell.edu
% -----------------------------------------------------

function [Hdenoised, num_iterations] = AN_denoiser(par, Hn, N0, denoiser)

Hdenoised = zeros(size(Hn));
num_iterations = zeros(par.U, 1);

for uu=1:par.U
    
    hnoinsy = Hn(:,uu);
    
    error_flag = true;
    N_trial = 0;
    max_num_trials = 100;
    while(error_flag)
        try
            N_trial = N_trial + 1;
            if N_trial > max_num_trials
                error(['ANM function failed even after ' num2str(max_num_trials) ' trials! Shutting down!']);
            end
            switch (denoiser)
                case 'AN_deb'
                    ast_out = ast_denoise(hnoinsy);
                    Hdenoised(:,uu) = ast_out.debiased;
                    num_iterations(uu) = ast_out.it_count;
                case 'AN_est'
                    ast_out = ast_denoise(hnoinsy);
                    Hdenoised(:,uu) = ast_out.estimate;
                    num_iterations(uu) = ast_out.it_count;
                case 'AN_deb_N0d2'
                    noise_std = sqrt(N0/2);
                    ast_out = ast_denoise(hnoinsy, 'noise_std',noise_std);
                    Hdenoised(:,uu) = ast_out.debiased;
                    num_iterations(uu) = ast_out.it_count;
                case 'AN_deb_N0'
                    noise_std = sqrt(N0);
                    ast_out = ast_denoise(hnoinsy, 'noise_std',noise_std);
                    Hdenoised(:,uu) = ast_out.debiased;
                    num_iterations(uu) = ast_out.it_count;
                case 'ANM'
                    noise_std = sqrt(N0);
                    ast_out = ast_denoise(hnoinsy, 'noise_std',noise_std);
                    Hdenoised(:,uu) = ast_out.debiased;
                    num_iterations(uu) = ast_out.it_count;
            end
            error_flag = false;
            %disp('Generated channel successfully!');
        catch
            warning('AN_denoiser.m: ANM failed in this trial. Trying again ...');
        end
    end
    
    
end
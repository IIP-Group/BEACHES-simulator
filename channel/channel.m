% -----------------------------------------------------
% -- LoS MIMO channel with some minumum angle seperation between users
% -- 2018 (c) studer@cornell.edu, sm2675@cornell.edu
% -- Last update: 10/14/2018
% -----------------------------------------------------

function H = channel(par)

error_flag = true;
N_trial = 0;
max_num_trials = 100;
while(error_flag)
    try
        [x_UE, y_UE] = users_locations(par);
        N_trial = N_trial + 1;
        if N_trial > max_num_trials
            error(['Channel generation failed after ' num2str(max_num_trials) ' trials! Shutting down!']);
        end
        switch (par.channel)
            case 'Rayleigh'
                H = sqrt(0.5)*(randn(par.B,par.U)+1i*randn(par.B,par.U));
            case 'FreeSpace'
                [~,H] = ch_LoS(par, x_UE,y_UE);
            case 'QuadMMLoS'
                H = ch_mmMAGIC_UMi_LOS(par,x_UE,y_UE);
            case 'QuadMMnLoS'
                H = ch_mmMAGIC_UMi_NLOS(par,x_UE,y_UE);
            otherwise
                error('par.channel not defined')
        end
        error_flag = false;
        %disp('Generated channel successfully!');
    catch
        warning('Channel generation failed. Trying again ...');
    end
end


end
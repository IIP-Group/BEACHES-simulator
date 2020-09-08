% -----------------------------------------------------
% -- Generate channels with Quadriga channel model and store them for 
% -- future simulations; or if the channels with the given paramters are
% -- are already generated, retrieve them.
% -- 2018 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -- Last update: 12/05/2018
% -----------------------------------------------------
function H = generate_channels(par, rq_num)
%% generate or retrieve channel realizations

if par.MAX_MEM < par.B*par.U*par.N_channel_realizations*16
    error('requested number of channel realizations is too large!')
end

% First check to see if there is a mat file containing channel realizations
% with the parameters given in par. If there is, retrieve it.
gen_channel_dir = [par.simulator_path 'channel/generated_channels/'];
filesAndFolders = dir(gen_channel_dir);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory                    
N_Files = length(filesInDir);
scenario = [num2str(par.B) '_' num2str(par.U) '_' num2str(par.channel) '_' num2str(par.MinAngleSep) '_' num2str(par.dmin) '_' num2str(par.dmax)];

found_flag = false;

for i = 1:N_Files
    filename = filesInDir(i).name;
    if ~isempty(strfind(filename, [scenario 'S' num2str(rq_num)]))
        H_load0 = load(filename);
        H_load = H_load0.H;
        found_flag = true;
        break;
    end
end


if found_flag && size(H_load,3) >= par.N_channel_realizations
    H = H_load;
    disp('Channel realizations loaded successfully!');
elseif found_flag && size(H_load,3) < par.N_channel_realizations
    disp(['There is not enough channel realizations in the stored file. Generating ' num2str(par.N_channel_realizations - size(H_load,3)) ' more realizations...']);
    H = zeros(par.B, par.U, par.N_channel_realizations);
    N3 = size(H_load,3);
    H(:,:,(1:N3)) = H_load;
    for chc_idx = 1:par.N_channel_realizations - N3
        H(:,:,N3 + chc_idx) = channel(par);
    end
    if par.SAVE_CHANNEL
        save([gen_channel_dir 'H_stored_', scenario 'S' num2str(rq_num) '.mat'], 'H');
    end
else
    disp(['Channel realizations with the given scenario not found. Generating ' num2str(par.N_channel_realizations) ' realizations...']);
    H = zeros(par.B, par.U, par.N_channel_realizations);
    for ch_idx = 1:par.N_channel_realizations
        H(:,:,ch_idx) = channel(par);
    end
    if par.SAVE_CHANNEL
    save([gen_channel_dir 'H_stored_', scenario 'S' num2str(rq_num) '.mat'], 'H');
    end
end


%% Normalize the channel entries

for i = 1:size(H, 3)
    Hstore_i = H(:,:,i);

    switch (par.ch_normalization)
    case 'overall'
        H(:,:,i) = Hstore_i/sqrt(mean(abs(Hstore_i(:)).^2));
    case 'per-user'
        for u = 1:par.U
            H(:,u,i) = Hstore_i(:,u)/sqrt(mean(abs(Hstore_i(:,u)).^2));
        end
    otherwise
        error('par.ch_normalization not defined')
    end
    
    if any(isnan(H(:,:,i)))
        error(['The ' num2str(i) '-th channel realization has NaN entries!'])
    end

end
    

end
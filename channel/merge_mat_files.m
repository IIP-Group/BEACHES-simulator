% -----------------------------------------------------
% -- Merge mat files storing channel realizations into one large file
% -- They were broken into small files only to be able to upload them to
% -- GitHub
% -----------------------------------------------------

fullpath = mfilename('fullpath');
file_path = strrep(fullpath, mfilename(), '');
REMOVE_ORIGINAL_SUB_FILES = true;

H = zeros(256,16,1000);

gen_channel_dir = [file_path '/generated_channels/'];


done = false;

while (~done)
    
    filesAndFolders = dir(gen_channel_dir);     % Returns all the files and folders in the directory
    filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
    N_Files = length(filesInDir);
    
    for i = 1:N_Files
        filename = filesInDir(i).name;
        if ~isempty(strfind(filename, '_sub'))
            header = regexprep(filename,'_sub\w*.mat','');
            done = false;
            break;
        end
        done = true;
    end
    
    if done
        return;
    end
    
    k = 1;
    for i = 1:N_Files
        filename = filesInDir(i).name;
        if ~isempty(strfind(filename, header)) && ~isempty(strfind(filename, '_sub'))
            load(filename);
            H(:,:,(k-1)*250+1:k*250) = Hsub;
            k = k + 1;
            delete(sprintf('%s',[file_path '/generated_channels/' filename]))
        end
        
    end
    
    save([file_path '/generated_channels/' header], 'H');
    
end

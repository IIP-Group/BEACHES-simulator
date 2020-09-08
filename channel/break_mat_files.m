% -----------------------------------------------------
% -- Break mat files storing channel realizations into smaller
% -- subsets, so that the file sizes are reduced to be uploaded to GitHub.
% -----------------------------------------------------

fullpath = mfilename('fullpath');
file_path = strrep(fullpath, mfilename(), '');
REMOVE_ORIGINAL_FILES = true;

file_names = { 'H_stored_256_16_QuadMMLoS_1_10_110S0',...
    'H_stored_256_16_QuadMMLoS_1_10_110S1',...
    'H_stored_256_16_QuadMMLoS_1_10_110S2',...
    'H_stored_256_16_QuadMMLoS_1_10_110S3',...
    'H_stored_256_16_QuadMMnLoS_1_10_110S0',...
    'H_stored_256_16_QuadMMnLoS_1_10_110S1',...
    'H_stored_256_16_QuadMMnLoS_1_10_110S2',...
    'H_stored_256_16_QuadMMnLoS_1_10_110S3',...
    };

for i = 1:length(file_names)
    load([file_names{i}]);
    Hsub = H(:,:,1:250);
    save([file_path '/generated_channels/' file_names{i} '_sub1'], 'Hsub')
    Hsub = H(:,:,251:500);
    save([file_path '/generated_channels/' file_names{i} '_sub2'], 'Hsub')
    Hsub = H(:,:,501:750);
    save([file_path '/generated_channels/' file_names{i} '_sub3'], 'Hsub')
    Hsub = H(:,:,751:1000);
    save([file_path '/generated_channels/' file_names{i} '_sub4'], 'Hsub')
end

if REMOVE_ORIGINAL_FILES
    for i = 1:length(file_names)
        delete(sprintf('%s.mat',[file_path '/generated_channels/' file_names{i}]))
    end
end

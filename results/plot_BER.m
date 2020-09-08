% -----------------------------------------------------
% -- plotter function
% -- 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
function plot_BER(figID, runID, curves, legends)

fullpath = mfilename('fullpath');
my_path = strrep(fullpath, mfilename(), '');
marker_style = {'bx-','ks-','rv-','go:','m*--','b>--','kx-', 'go--', 'rd-', 'cx:'};
% Check to see if the simulation results for the specified runID exists
gen_channel_dir = [my_path 'sim_res/'];
filesAndFolders = dir(gen_channel_dir);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
N_Files = length(filesInDir);
header = ['ERR-' num2str(runID) '_'];
n_curves = numel(curves);

for curve_idx = 1:n_curves
    
    found_flag = false;
    
    for i = 1:N_Files
        filename = filesInDir(i).name;
        if ~isempty(strfind(filename, header)) && ~isempty(strfind(filename, [curves{curve_idx} '.mat']))
            load(filename);
            found_flag = true;
            break;
        end
    end
    
    if ~found_flag
        error('The requested simulation result not found!');
    end

    % normalize results
    res.BER = res.BER./res.num_trials;

    % -- show results (generates fairly nice Matlab plot)
    
    
    figure(figID)   
    if curve_idx==1
        semilogy(par.SNRdB_list,res.BER,marker_style{curve_idx},'LineWidth',2, 'DisplayName', legends{curve_idx});
        hold on
    else
        semilogy(par.SNRdB_list,res.BER,marker_style{curve_idx},'LineWidth',2, 'DisplayName', legends{curve_idx});
    end
    
  
end

figure(figID)
hold off
grid on
xlabel('signal-to-noise ratio [dB]','FontSize',12)
ylabel('uncoded bit error rate','FontSize',12)
ylim([1e-4 1])
xlim([-10 15])

plot_title = ['B = ' num2str(par.B) ', U = ' num2str(par.U) ', ' par.mod ', ' par.channel];
title(plot_title);
legend('Location','northeast')
legend('show', 'FontSize',12);
set(gca,'FontSize',12)
savefig([my_path 'figures/' 'BER_' plot_title])

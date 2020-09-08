% ------------------------------------------------------------------------
% -- BEACHES denoiser with hardware optimizations that do not affect its 
% -- performance. The optimization is that only the entries of the sorted
% -- input noisy channel are candidates of the threshold value, without
% -- finding the optimal threshold value between two consecutive sorted
% -- entries. Please refer to the paper for more details.
% -- 2018 (c) rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu
% ------------------------------------------------------------------------

function Hdenoised = BEACHES_hw(par,Hn, N0)

hdenoised = zeros(size(Hn));
SURE = zeros(par.B,1);
tau_max = inf;
for uu=1:par.U
    hnoisy = fft(Hn(:,uu))/sqrt(par.B);
    N = length(hnoisy);
    % find threshold for shrinkage function in NlogN time (sort + linear scan)
    
    sorth = sort(abs(hnoisy),'ascend');
    cumsum2 = 0;
    cumsuminv = sum(1./abs(hnoisy));
    tau_opt = inf;
    suremin = inf;
    for bb = 1:N
        tau = sorth(bb);
        SURE(bb) = cumsum2 + (N-bb+1)*tau^2 + N*N0 - 2*N0*(bb-1)-tau*N0*cumsuminv;
        cumsum2 = cumsum2 + sorth(bb).^2;
        cumsuminv = cumsuminv - 1/sorth(bb);
        if SURE(bb)<suremin
            suremin = SURE(bb);
            tau_opt = tau;
        end
    end
    tau_opt = min(tau_opt, tau_max);

    
    hdenoised(:,uu) = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    
end

Hdenoised = ifft(hdenoised)*sqrt(par.B);

end
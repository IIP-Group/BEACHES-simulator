% ------------------------------------------------------------------------
% -- SOFT thresholding with MSE-optimal threshold
% -- 2018 (c) rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu
% ------------------------------------------------------------------------

function Hdenoised = exactMSE(par,Hn, H)

hdenoised = zeros(size(Hn));
MSE = zeros(par.B,1);

for uu=1:par.U
    hnoisy = fft(Hn(:,uu))/sqrt(par.B);
    h = fft(H(:,uu))/sqrt(par.B);
    N = length(hnoisy);
    
    sorth = sort(abs(hnoisy),'ascend');
    for i = 1:N
        tau = sorth(i);
        hdenoised_t = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau,0);
        MSE(i) = norm(h-hdenoised_t)^2;
    end
    [~, tau_index] = min(MSE);
    tau_opt = sorth(tau_index);
    
    hdenoised(:,uu) = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    
end

Hdenoised = ifft(hdenoised)*sqrt(par.B);

end
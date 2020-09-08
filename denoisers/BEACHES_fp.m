% -----------------------------------------------------
% -- Fast SURE denoiser
% -- 2018 (c) rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu
% ******************** How does it work? ***********************
% Uses cart2pol and pol2cart instead of z/|z|.
% Plus more detailed fixed point simulation of variables.
% **********************************************************************
function Hdenoised = BEACHES_fp(par,Hn, H, N0_O)

H_W = 16;    % Total number of bits for the entries of Hn
H_F = 8;    % number of fractional bits for the entries of Hn
hf_W = 10;
hf_F = 8;
Hn = fixP(Hn, H_W, H_F);
hdenoised = zeros(size(Hn));
SURE = zeros(par.B,1);
MSE = zeros(par.B,1);
N0 = fixP(N0_O/(par.B),16,15);

for uu=1:par.U
    
    hnoisy1 = fft(Hn(:,uu))/(par.B);
    %hnoisyMax(uu) = max(max(real(hnoisy1)), max(imag(hnoisy1)));
    %hnoisyMin(uu) = min(min(real(hnoisy1)), min(imag(hnoisy1)));
    hnoisy = fixP(hnoisy1, hf_W, hf_F);
    [anghn, abshn] = cart2pol(real(hnoisy), imag(hnoisy));
    anghnoisy = fixP(anghn, hf_W, hf_F-1);  % 3.7 format
    abshnoisy = fixP(abshn, hf_W, hf_F);    % 2.8 format
    
    h = fft(H(:,uu))/sqrt(par.B);
    N = length(hnoisy);
    
    % find threshold for shrinkage function in NlogN time (sort + linear scan)
    
    sorth = sort(abshnoisy,'ascend');
    sum2 = 0;
    % -- reciprocal unit:
    recip = fixP(1./sorth, 12,2);
    % ----------- Solution for numerical problem: -----------
    recip(sorth < 2^(-hf_F)) = 0;
    suminv = fixP(sum(recip),20,2);
    suremin = inf;
    for bb=1:length(sorth)
        tau = sorth(bb);
        tau2 = fixP(sorth(bb)^2, 20,16);
        sum2 = fixP(fixP(sum2,26,16) + tau2,26,16);
        T1pT2 = fixP(fixP(sum2, 18,8) + fixP((N-bb+1)*tau2, 18,8),18,8);
        T3pT4 = fixP(fixP(2*N0*(bb-1),17,8) + fixP(tau*fixP(N0*suminv,24,8),24,8),24,8);
        SURE(bb) = fixP(T1pT2 - T3pT4, 24,8);
        %SURE(bb) = fixP((fixP(sum2, 24,8) + fixP((N-bb+1)*tau2, 24,8) - 2*N0*(bb-1)-fixP(tau*fixP(N0*suminv,24,8),24,8)), 24,8);
        %sum2 = fixP(fixP(sum2,26,16) + tau2,26,16);
        suminv = suminv - recip(bb);
        if SURE(bb)<suremin
            suremin = SURE(bb);
            tau_opt = tau;
        end
    end
    
    
    
    [x,y] = pol2cart(anghnoisy, max(abshnoisy-tau_opt,0));
    hndenoisedtemp = x + 1i*y;
    %hndenoisedtemp(abshnoisy < 2^(-hf_F)) = 0;
    hdenoised(:,uu) = fixP(hndenoisedtemp, hf_W+1, hf_F+1);
    
end

Hdenoised = ifft(hdenoised)*(par.B);
Hdenoised = fixP(Hdenoised, H_W, H_F);

end




function A_FP = fixP(A, W, F)

A1 = floor(A*2^F)/2^F;
A_FP = sign(A).*min(abs(A1), 2^(W-F-1));

end
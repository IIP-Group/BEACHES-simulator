% =========================================================================
% Newtonized OMP. The extractSpectrum function and examples downloaded from:
% https://bitbucket.org/wcslspectralestimation/continuous-frequency-estimation/src/NOMP
% 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% =========================================================================

function [Hdenoised, num_freq] = NOMP_denoiser(par,Hest, N0_est)

N = par.B;
sigma = sqrt(N0_est);
num_freq = zeros(par.U,1);
S = eye(N);
tau = sigma^2 * (log(N) - log( log(1/(1-par.p_fa))));

Hdenoised = zeros(size(Hest));
for u = 1:par.U
    [omega_est, gain_est, ~] = extractSpectrum(Hest(:,u), S, tau);
    if isempty(omega_est)
        hdenoised = Hest(:,u);
    else
        num_freq(u) = length(omega_est);
        hdenoised = exp(1i* (0:(N-1)).' * omega_est.')/sqrt(N) * gain_est;
    end
    Hdenoised(:,u) = hdenoised;
end

end
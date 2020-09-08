% ------------------------------------------------------------------------
% -- optimal BEACHES denoiser. Refer to the BEACHES paper for more details. 
% -- 2018 (c) rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu
% ------------------------------------------------------------------------

function Hdenoised = BEACHES(par,Hn, H, E0, denoiser,show)

hdenoised = zeros(size(Hn));
SURE = zeros(par.B,1);
MSE = zeros(par.B,1);
e = 0.01;
CORDIC_SCALE = 1.0; 
%tau_max = 4*sqrt(-E0*log(e)); % Fixed threshold upper bound
tau_max = inf;
for uu=1:par.U
    hnoisy = CORDIC_SCALE*fft(Hn(:,uu))/sqrt(par.B);
    h = fft(H(:,uu))/sqrt(par.B);
    N = length(hnoisy);
    switch (denoiser)
        
        case 'BEACHES' % complex SURE
            % find threshold for shrinkage function in NlogN time (sort + linear scan)
            
            sorth = sort(abs(hnoisy),'ascend');
            cumsum2 = 0;
            cumsuminv = sum(1./abs(hnoisy));
            tau_opt = inf;
            suremin = inf;
            for bb = 1:N
                tau = sorth(bb);
                SURE(bb) = cumsum2 + (N-bb+1)*tau^2 + N*E0 - 2*E0*(bb-1)-tau*E0*cumsuminv;
                cumsum2 = cumsum2 + sorth(bb).^2;
                cumsuminv = cumsuminv - 1/sorth(bb);
                if SURE(bb)<suremin
                    suremin = SURE(bb);
                    tau_opt = tau;
                end
            end
            tau_opt = min(tau_opt, tau_max);
            
        case 'BEACHESopt' % complex SURE
            % find optimal threshold for shrinkage function
            
            sorth = sort(abs(hnoisy),'ascend');
            cumsum = 0;
            cumsuminv = sum(1./abs(hnoisy));
            tau_opt = inf;
            suremin = inf;
            tau_low = 0;
            if(show)
                tau_opt_k = 0;
                k_opt = 1;
            end
            for bb = 1:N
                tau_high = sorth(bb);
                tau = max(tau_low,min(tau_high,E0/(2*(N-bb+1))*cumsuminv));
                tau_low = sorth(bb);
                SURE(bb) = cumsum + (N-bb+1)*tau^2 + N*E0 - 2*E0*(bb-1)-tau*E0*cumsuminv;              
                cumsum = cumsum + sorth(bb).^2;
                cumsuminv = cumsuminv - 1/sorth(bb);
                if SURE(bb)<suremin
                    suremin = SURE(bb);
                    tau_opt = tau;
                end
                if((show) && (tau~=tau_high)) % save the values for display when BEACHES != BEACHESopt
                    tau_opt_k = tau;
                    k_opt = bb;
                end            
            end
            if((show)&&(tau_opt == tau_opt_k)) % this is just for displaying (when show = 1)
                if(k_opt==1)
                    fprintf(' k_opt = %3.0f \t \t tau_low = %6.4f \t tau_opt = %6.4f \t tau_high = %6.4f\n', k_opt, 0, tau_opt, sorth(k_opt));
                    fprintf('             \t SURE(tau_low) = %6.4f \t SURE(tau_opt) = %6.4f \t SURE(tau_high) = %6.4f\n',...
                        N*E0,...
                        sum(sorth(1:k_opt).^2)+(N-k_opt+1)*tau_opt^2 + N*E0 - 2*E0*(k_opt-1)-tau_opt*E0*sum(1./sorth(k_opt+1:N)),...
                        sum(sorth(1:k_opt).^2)+(N-k_opt+1)*sorth(k_opt)^2 + N*E0 - 2*E0*(k_opt-1)-sorth(k_opt)*E0*sum(1./sorth(k_opt+1:N)));
                else
                    fprintf(' k_opt = %3.0f \t \t  tau_low  = %6.4f \t \t  tau_opt  = %6.4f \t \t  tau_high  = %6.4f\n', k_opt, sorth(k_opt-1), tau_opt, sorth(k_opt)); 
                    fprintf('             \t SURE(tau_low) = %6.4f \t SURE(tau_opt) = %6.4f \t SURE(tau_high) = %6.4f\n',...
                        sum(sorth(1:k_opt-1).^2)+(N-k_opt)*sorth(k_opt-1)^2 + N*E0 - 2*E0*(k_opt-1-1)-sorth(k_opt-1)*E0*sum(1./sorth(k_opt:N)),...
                        sum(sorth(1:k_opt).^2)+(N-k_opt+1)*tau_opt^2 + N*E0 - 2*E0*(k_opt-1)-tau_opt*E0*sum(1./sorth(k_opt+1:N)),...
                        sum(sorth(1:k_opt).^2)+(N-k_opt+1)*sorth(k_opt)^2 + N*E0 - 2*E0*(k_opt-1)-sorth(k_opt)*E0*sum(1./sorth(k_opt+1:N)));
                end
            end  
            
        case 'MSE-SOFT'
            sorth = sort(abs(hnoisy),'ascend');
            for i = 1:N
                tau = sorth(i);
                hdenoised_t = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau,0);
                MSE(i) = norm(h-hdenoised_t)^2;
            end
            [~, tau_index] = min(MSE);
            tau_opt = sorth(tau_index);
            
        otherwise
            error('par.denoiser not defined');
    end

    
    
    hdenoised(:,uu) = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    
end

Hdenoised = ifft(hdenoised/CORDIC_SCALE)*sqrt(par.B);

end
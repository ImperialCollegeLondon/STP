function [ arma_noise_var ] = estimate_ARMA_noise_var( phi,theta,sig_var )

%estimate_ARMA_noise_var estimates the noise variance of a ARMA process
%with predescribed statistics and variance
options = optimoptions('fsolve');
options.Display='off';
arma_noise_var = fsolve(@(x) obj(x,phi,theta,sig_var),0.01,options);


    function err = obj(xx,phi,theta,sig_var)

        [~, sigma2_y] = acf(phi,theta,1,xx);

        err = sig_var - sigma2_y;
        
    end


end


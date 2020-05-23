function [ ea , dDe ] = ComputeVapPressure2( esat , Ta , Rswtm1 , Rswtm2 , dDetm1 , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe )
%	Inputs:
%   esat - vapor pressure at saturation [Pa]
%   Ta - temperature [deg. C]
%   Rswtm1 , Rswtm2 - global shortwave radiation [W m^-2] at t-1 and t-2
%   dDetm1 - stochastic component of vapor pressure deficit [Pa] at t-1
%	a0...a3 - regression coefficients of deterministic component vapor pressure
%	deficit [-]
%	dDem - average of vapor pressure deficit devaitions
%	rhodDe - lag-1 autocorrelation of the process
%	sigmadDe - standard deviation of the process
%   Outputs:
%   ea - air ambient vapor pressure [Pa]
%   dDe - stochastic component of vapor pressure deficit [Pa]
%% Compute vapor pressure
% Dedet=a0+a1.*(Ta.^3)+a2.*Rswtm1+a3.*Rswtm2; % Slower
Dedet=a0+a1.*(Ta.*Ta.*Ta)+a2.*Rswtm1+a3.*Rswtm2; 
epsDe=normrnd(0,1);
dDe=dDem+rhodDe.*(dDetm1-dDem)+epsDe.*sigmadDe.*sqrt(1-rhodDe.^2);
dDe=mean2(dDe);
De=Dedet+dDe;
ea=single(esat-De); 
ea(ea<0)=1;
ea(ea>esat)=esat(ea>esat);
end
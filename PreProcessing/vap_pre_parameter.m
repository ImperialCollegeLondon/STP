function [ a , dDem , rhodDe , sigmadDe ] = vap_pre_parameter( ea , esat , Ta , Rsw )
% TEMPERATURE PARAMETER ESTIMATION - Monthly basis
%	Inputs:
%   ea - air ambient vapor pressure [Pa]
%   esat - vapor pressure at saturation [Pa]
%   Ta - temperature [deg. C]
%   Rsw - global shortwave radiation [W m^-2]
%	Outputs:
%	a - regression coefficients of deterministic component vapor pressure
%	deficit [-]
%	dDem - average of vapor pressure deficit devaitions
%	rhodDe - lag-1 autocorrelation of the process
%	sigmadDe - standard deviation of the process
%% Initilzing
n=length(ea);
ea=reshape(ea,1,n);
esat=reshape(esat,1,n);
Rsw=reshape(Rsw,1,n);
Ta=reshape(Ta,1,n);
Ta3=sign(Ta).*abs(Ta).^3;
Y=esat(3:end)-ea(3:end);
%% Coefficients
XX=[ones(1,length(Y)); Ta3(3:end); Rsw(2:end-1); Rsw(1:end-2);];
a=regress(Y',XX');
%% Parameters
Des=a'*XX;
Des(Des<0)=0;
dDe=Y-Des;
dDe=dDe(not(isnan(dDe)));
dDem=mean(dDe);
R=autocorr(dDe,10);
rhodDe=R(2);
sigmadDe=std(dDe);
end
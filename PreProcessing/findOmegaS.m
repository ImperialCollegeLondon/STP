function [ omegaS ] = findOmegaS( betaS , tet , xiS )
%findOmegaS Find the slope in the direction of the wind omegaS.
%Following Liston and Elder (2006).
%   Inputs:
%   betaS - maximum slope for each pixel [deg]
%   tet - wind direction [deg]
%   xiS - terrain slope azimuth [deg from north]
%   Output:
%   omegaS - the slope in the direction of the wind
%% Find slope direction
omegaS=betaS.*cosd(180+tet-xiS);
%% Scale to [-0.5 0.5]
tmp=reshape(omegaS,[],1);
tmp(isnan(tmp))=[];
[muhat,sigmahat] = normfit(tmp);
tmp=normcdf(omegaS,muhat,sigmahat);
omegaS=-0.5+tmp;
omegaS(omegaS<-0.5)=-0.5;
omegaS(omegaS>0.5)=0.5;
end
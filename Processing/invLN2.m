function [ rField ] = invLN2( qField , war , imf , CV )
%invLN A function that convert the quantile field to a homogenoues
%precipitation field using a inverse LogNormal distribution.
%   Inputs:
%   qField - normalized quantile field
%   war - wet area ration [-]
%   imf - areal mean intensity [mm h^-1]
%   CV - coefficient of variation for the precipitation
%   Outputs:
%   rField - rain field [mm h^-1]
%% Initilazing
miu=bsxfun(@times,log(imf./(war.*sqrt(CV^2+1))),ones(size(qField)));
s=sqrt(log(CV^2+1));
%% Invrese LogNormal - rainfall assingments
rField=logninv(qField,miu,s);
end
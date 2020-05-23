function [ cos_ST ] = SolarIlluminationAngle( SLO , Alt , zeta_S , zeta_T )
%SolarIlluminationAngle Compute the Solar Illumination Angle
%   Inputs:
%   SLO - slope [-]
%   Alt - solar altitude [rad]
%   zeta_S - azimute angle [rad]
%   zeta_T - local aspect (clock wise direction from the north) [angular degree from North]
%   Output:
%   cos_ST - local solar illumination angle [rad]
%% Initilazing
SLO=atand(SLO);
SLO=SLO.*pi./180;
zeta_T=zeta_T.*pi./180;
cos_ST=cos(SLO).*sin(Alt)+sin(SLO).*cos(Alt).*cos(zeta_S-zeta_T);
end
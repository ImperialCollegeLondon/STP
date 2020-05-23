function [ EB , ED , EB1 , EB2 , ED1 , ED2 , PARB , PARD ] = ComputeRadiationForcings( h_S , E0 , Ro , LWP0 , N , domainDTM , Tdew , beta_A , alpha_A , omega_A1 , omega_A2 , uo , un , rho_g , cos_ST , S , SvF , Ct , corrAOD , zAOD )
%	Inputs:
%   h_S - solar altidute [rad]
%   E0 - correction factor
%   Ro - Solar constant [W/m^2]
%   LWP0 - liquid water path [g/m^2]
%	N - cloud fraction [-]
%   domainDTM - elevation map [m]
%   Tdew - dew point temperature [°C] 
%   beta_A - Angstrom turbidity parameter [-]
%   alpha_A - Angstrom turbidity parameter [-]
%   omega_A1 - aerosol single-scattering albedo band 1 [0.74-0.94]
%   omega_A2 - aerosol single-scattering albedo band 2 [0.74-0.94]
%	uo - ozone amount in vertical column [0.22-0.35] [cm]
%	un - total nitrogen dioxide amount [0.0001-0.046] [cm]
%	rho_g - spatial average regional albedo [0.05-0.3]
%   cos_ST - local solar illumination angle [rad]
%   S - shadow effect [boolean, 0 is shaded]
%   SvF - sky view factor [-]
%   Ct - terrain factor [-]
%   corrAOD - elevation correction for the aerosol optical depth (after
%   Ingold, 2001)
%   zAOD - elevation of the aerosol optical depth station [m]
%   Outputs:
%	EB - sky beam irradiance without terrain effects [W/m^2] 
%	ED - sky total diffuse irradiance without terrain effects [W/m^2] 
%   EB1 - sky beam irradiance VIS band[0.29 um - 0.70um ] [W/m^2] 
%	EB2 - sky beam irradiance NIR band [0.70 um - 4.0 um ] [W/m^2] 
%   ED1 - sky total diffuse flux at the ground  VIS band [0.29 um -0.70um ] [W/m^2] 
%	ED2 - sky total diffuse flux at the ground NIR band  [0.70 um-4.0 um ] [W/m^2] 
%% Constants
Sop=Ro.*E0; % Actual Solar constant [W/m^2]
%% Partition Energy two bands Gueymard (2004; 2008)
So1=Sop*0.4651;  % Extraterrestrial radiation VIS band [0.29 um - 0.70 um ] [W m^-2]
So2=Sop*0.5195; % Extraterrestrial radiation NIR band [0.70 um - 4.0 um ] [W m^-2]
%% Cloudiness
LWP=LWP0*N; % Liquid Water Path while cloudiness  [g/m^2]
if LWP<1.1
    Ntmp=0;
else
    Ntmp=N;
end
[ beta_A ] = AODexp( beta_A , domainDTM , zAOD , corrAOD );
if Ntmp==0 % Completely clear sky
    [EB1,EB2,ED1,ED2,~,~,~,~,Mb,Mg] = SetClearSkyRadiation(h_S,domainDTM,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
elseif Ntmp>0 % Overcast Condition
    [Eb1,Eb2,~,~,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,domainDTM,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    [EB1,EB2,ED1,ED2] = SetCloudySkyRadiation(LWP,h_S,Eb1,Eb2,Edp1,Edp2,N,rho_s1,rho_s2,rho_g);
end
%% Correction slope
EB1=S.*cos_ST.*EB1;
EB2=S.*cos_ST.*EB2;
%% Correcting the sky view factor
ED1=ED1.*SvF;
ED2=ED2.*SvF;
%% Reflected radiation
ER1=Ct.*rho_g.*(EB1.*cos_ST+(1-SvF).*ED1);
ER2=Ct.*rho_g.*(EB2.*cos_ST+(1-SvF).*ED2);
%% Summing the bands
EB = EB1 + EB2 + ER1; % beam irradiance [W/m^2]
ED = ED1 + ED2 + ER2; % total diffuse flux at the ground [W/m^2]
%% PAR Radiation estimation
PARB = EB1.*Mb;
PARD = Mg.*(EB1+ED1)-PARB;
end
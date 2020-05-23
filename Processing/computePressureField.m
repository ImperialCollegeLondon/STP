function [ P ] = computePressureField( P0 , Z0 , Ta , domainDTM )
%computePressureField Compute the pressure field based on a reference
%location
%   Inputs:
%   P0 - pressure at refernce point [hPa]
%   Z0 - elevation at reference point [m]
%   Ta - Temperature field [Deg Celcius]
%   domainDTM - elevation field [m]
%   Output:
%   P - Pressure field [hPa]
%% Initilazing
g=9.80665; % Gravitational acceleration constant [m/s]
M=0.0289644; % molar mass of Earth's air [kg/mol]
R=8.31432; % Universal gas constant [N*m/mol*deg.K]
%% Compute pressure field
TaK=Ta+273; % [Deg. Kelvin]
T0=TaK(domainDTM==Z0);
T0=T0(1);
L0=(TaK-T0)./(domainDTM-Z0);
P=P0.*(1+(L0./T0).*(domainDTM-Z0)).^((-g.*M)./(R.*L0));
P(isnan(L0))=P0;
end
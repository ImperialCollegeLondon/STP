function [ u , v ] = p2z( tet , W )
%p2z Polar to zonal winds
%   Inputs:
%   tet - wind direction [deg]
%   W - wind velocity [m s^-1]
%   Output:
%   u , v - zonal wind component [m s^-1]
%% Code
u=W.*sind(tet);
v=W.*cosd(tet);
end
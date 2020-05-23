function [ W , tet ] = z2p( u , v )
%z2p Zonal to polar winds
%   Inputs:
%   u , v - zonal wind component [m s^-1]
%   Output:
%   tet - wind direction [deg]
%   W - wind velocity [m s^-1]
%% Code
W=sqrt(u.^2+v.^2);
tet=zeros(size(W));
tet(u==0&v<0)=180;
tet(u>0&v==0)=90;
tet(u<0&v==0)=270;
tet(u>0&v>0)=atand(u(u>0&v>0)./v(u>0&v>0));
tet(u<0&v<0)=180+atand(u(u<0&v<0)./v(u<0&v<0));
tet(u>0&v<0)=180+atand(u(u>0&v<0)./v(u>0&v<0));
tet(u<0&v>0)=360+atand(u(u<0&v>0)./v(u<0&v>0));
end
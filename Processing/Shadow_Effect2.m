function [ ShF ] = Shadow_Effect2( DTM , S_zen , zeta_S , HZr , Z )
%	Inputs:
%	DTM - matrix digital elevation model [m]
%	S_zen - solar height [rad]
%	zeta_S - solar azimuth [rad]
%	HZr - Horizon angle [rad] from N
%	Z - Azimuth [angular degree] from N
%	Output:
%	ShF - Shadow Effect 1 --> no shadow 0 --> shadow
%% Initilazing
Z= Z*pi/180; %%% [rad]
Z=[Z 2*pi]; %%[rad]
ShF=zeros(size(DTM));
%% Vec
vec_m=(abs(Z-zeta_S));
[~,p1]=min(vec_m);
vec_m(p1)=Inf;
[~,p2]=min(vec_m);
p3=p1;
if p3==length(Z)
    p3=1;
end
p4=p2;
if p4==length(Z)
    p4=1;
end
%% Find ShF
HZz1=HZr(:,:,p3);
HZz2=HZr(:,:,p4);
HZi=HZz1+(HZz2-HZz1).*((zeta_S - Z(p1))./(Z(p2)-Z(p1)));
ShF(HZi>S_zen)=1;
ShF(isnan(DTM))=NaN;
end
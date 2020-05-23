function [ LR , dLR , T0 , T0amp , Th1 , Th2 ] = computeLR( T , e1 , e2 , m )
%computeLR Function that compute the lapse rate for two given ground stations
%   Inputs:
%   T - Temperature [deg Celsius]
%   e1 and e2 - elevations of the ground station [m a.s.l.]
%   m - month to analyze
%   Outputs:
%   LR - lapse rate [deg. C per km, usually negative with increasing hight]
%   dLR - date of observation
%   T0 - temperature [deg. C] at 0 elevation a.s.l.
%   T0amp - amplitude between max and min daily temperature [deg. C]
%% Compute LR
Th1=T(month(T(:,3))==m,1);
Th2=T(month(T(:,3))==m,2);
diffE=(e2-e1)/1000;
diffT=Th2-Th1;
LR=diffT/diffE;
%% Daily dates
dLR=T(month(T(:,3))==m,3);
tmp=datevec(dLR(:,1));
dLR=datenum(tmp(:,1),tmp(:,2),tmp(:,3),0,0,0);
%% Temperature [deg. C] at 0 elevation
T0=Th1(:)+(e1/1000).*-LR(:);
for i=1:size(Th1,1)
    aCorrMax(i,1)=max(T0(dLR==dLR(i)));
    aCorrMin(i,1)=min(T0(dLR==dLR(i)));
end
T0amp=aCorrMax-aCorrMin;
%% Dates
dLR=T(month(T(:,3))==m,3);
end
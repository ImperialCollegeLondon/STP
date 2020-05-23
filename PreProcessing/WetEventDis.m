function [ cloudDis ] = WetEventDis( clt , T , dt )
%WetEventDis A function that calculates the distance (in time) from the closest
%wet event.
%   Inputs:
%   clt - cloud time series [fraction]
%   T - a time series indicating a wet [1] or dry [0] event
%   dt - time interval [min]
%   Outputs:
%   cloudDis - time series of distance of each dry time step from the
%   closest wet event [min]
%% Initilazing
cloudDis=zeros(size(clt));
%% Find cloudiness for distance from a wet event
parfor i=1:size(clt,1)
    if T(i)==0
        cloudDis(i)=min(abs(i-find(T>=1)))*dt; % Time [min] to / from closest rain event
    end
end
end
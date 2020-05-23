function [ radarIsEvent ] = findRadarEvents( radarWAR , thresholdWAR , thresholdTime )
%findRadarEvents Define wet event base on the radar data
%   Inputs:
%   radarWAR - wer area ratio [-]
%   thresholdWAR - wet area threshold to consider as wet event [-]
%   threshold time - minimum time of wet event lentgh [min]
%   Output:
%   radarIsEvent - time series of the wet events indexed by ID
%% Initilazing
radarIsEvent=zeros(size(radarWAR));
WAR(radarWAR>=thresholdWAR,1)=1;
for i=6:size(WAR,1)-6 % Correcting time series of radar data for wet events hiatus
    if sum(WAR(i-5:i-1))>=1 && sum(WAR(i+1:i+5))>=1 && WAR(i)==0
        WAR(i)=1;
    end
end
es=[]; ee=[]; c=1;
%% Finding wet events
if WAR(1)==1
    es(size(es,1)+1,1)=1;
end
for i=2:size(WAR,1)
    if WAR(i)==1 && WAR(i-1)==0
        es(size(es,1)+1,1)=i;
        radarIsEvent(i,1)=c;
    end
    if WAR(i)==1 && WAR(i-1)==1
        radarIsEvent(i,1)=c;
    end
    if WAR(i)==0 && WAR(i-1)==1
        ee(size(ee,1)+1,1)=i-1;
        c=c+1;
    end
    if i==size(WAR,1) && WAR(i)==1 && WAR(i-1)==1
        ee(size(ee,1)+1,1)=i-1;
    end
end
%% Checking for event length treshold
checkTime=(ee-es)<thresholdTime/5;
for i=1:size(checkTime,1)
    if checkTime(i)==1
        radarIsEvent(radarIsEvent==i)=0;
    end
end
end
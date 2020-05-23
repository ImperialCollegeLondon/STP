function [ gaugeIsEvent ] = findGaugeEvents( gauge , thresholdG , thresholdTime )
%findGaugeEvents Define wet event base on the gauge data
%   Inputs:
%   gauge - rainfall [mm] per 10-min
%   thresholdG - rainfall threshold to consider [mm]
%   threshold time - minimum time of wet event lentgh [min]
%   Output:
%   gaugeIsEvent - time series of the wet events indexed by ID
%% Initilazing
gaugeIsEvent=zeros(size(gauge,1),1);
G(gauge>=thresholdG,1)=1;
for i=3:size(G,1)-5 % Correcting time series of gauge data for wet events hiatus, i.e. up to 20 minutes dry period btw wet events consider also wet
    if sum(G(i-1))>=1 && sum(G(i+1:i+5))>=1 && G(i)==0
        G(i)=1;
    end
end
es=[]; ee=[]; c=1;
%% Finding wet events
if G(1)==1
    es(size(es,1)+1,1)=1;
end
for i=2:size(G,1)
    if G(i)==1 && G(i-1)==0
        es(size(es,1)+1,1)=i;
        gaugeIsEvent(i,1)=c;
    end
    if G(i)==1 && G(i-1)==1
        gaugeIsEvent(i,1)=c;
    end
    if G(i)==0 && G(i-1)==1
        ee(size(ee,1)+1,1)=i-1;
        c=c+1;
    end
    if i==size(G,1) && G(i)==1 && G(i-1)==1
        ee(size(ee,1)+1,1)=i-1;
    end
end
%% Checking for event length treshold
checkTime=(ee-es)<thresholdTime/10;
for i=1:size(checkTime,1)
    if checkTime(i)==1
        gaugeIsEvent(gaugeIsEvent==i)=0;
    end
end
end
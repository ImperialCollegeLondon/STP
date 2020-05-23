function [ tsMatrix ] = tsTable( tsWetDry , tsMonth )
%tsTable Preaper a table fotmat of the dry / wet periods time series
%% Generating the time series table
tsMatrix(1,1)=tsWetDry(1); % Wet or dry condition
tsMatrix(1,2)=tsMonth(1); % Month of the start of the event
tsMatrix(1,3)=1; % Start of event
c=2;
for i=2:size(tsWetDry,2)
    if tsWetDry(i)~=tsWetDry(i-1)
        tsMatrix(c,1)=tsWetDry(i);
        tsMatrix(c,2)=tsMonth(i);
        tsMatrix(c,3)=i;
        tsMatrix(c-1,4)=i-1;
        c=c+1;
    end
end
tsMatrix(end,4)=size(tsWetDry,2);
end
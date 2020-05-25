
clear;clc

load('H:\CODE_MATLAB\calEvents.mat','calEvents','speed','EventSUMMARY','ObsSUMMARY');
load('H:\CODE_MATLAB\calEvents_Lagrangian.mat','calEvents_Lagrangian');

%%

colm = flip(pink(80),1);
colm(1:20,:) = [];
evi = 1;
for calEvent = calEvents_Lagrangian'
    
    [rain] = originalData(calEvent,'double');
    timestep = 12*2; % 2 hours 
    [autoCorr,lags,dur] = getACF(rain,timestep);
    
    plot(lags,nanmean(autoCorr,1),'color',colm(evi,:));
    
    hold on;
    drawnow
    evi = evi+1;
end
xlabel('Lag');
ylabel('acf in lagrangian coor');
savePlot(['Figs\Lagrangian\TemporalEvol'],...
        'wholepage',true,'onlyPng',true,'needreply','N');
%%
evi = 0;
Rains = [];
for calEvent = calEvents_Lagrangian'
    [R] = originalData(calEvent,'double');
    [RTime] = getTime(calEvent);
    Rmon = RTime(1).Month;
    % remove those events with missing data
    if any(isnan(R(:)))
        % # to do #
        % correct NaN val into ...
    else
        evi = evi+1;
        Rains{evi,1} = R;
    end
end
[ ARMA, mean_spatial_correlation ] = ARMA_estimate( Rains, 1, true );
%%
try
    [ ARMA, mean_spatial_correlation ] = ARMA_estimate( Rains, 1, true );
    pause(1)
    savePlot(['Figs\Lagrangian\event',num2str(evi),'-ARMAfit'],...
        'wholepage',true,'onlyPng',true,'needreply','N');
    close(gcf);
catch me
    fprintf('>>> Check Event No.%03d\n', evi);
end

save('Birmingham\ARMA.mat','ARMA','Rains')





%%
for i=0:4
    for j=0:4
        Mdl=arima(i,0,j);
        [~,~,LLF] = estimate(Mdl,G);
        NumActive = i+j+1;
        AIC(i+1,j+1) = aicbic(LLF,NumActive);
    end
end

% Choose the model to use
[i,j]=find(min(min(AIC))==AIC);
Mdl=arima(i-1,0,j-1);



function [autoCorr,lags,dur] = getACF(rain,timestep)
% dur: unit-[h]
autoCorr = [];
[autoCorr,lags,dur] = deal([]);
dur = size(rain,3)/12;

rain = reshape(rain,[],size(rain,3));
for i = 1:size(rain,1)
    [autoCorr(i,:),lags,~] = autocorr(rain(i,:),timestep);
end

end


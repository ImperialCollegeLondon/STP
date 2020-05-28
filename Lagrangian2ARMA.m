
clear;clc
load('H:\CODE_MATLAB\calEvents.mat','calEvents','speed','EventSUMMARY','ObsSUMMARY');
load('H:\CODE_MATLAB\calEvents_Lagrangian.mat','calEvents_Lagrangian');
%% Obtain R field.
evi = 0;
[Rains,Rmon] = deal([]);
for calEvent = calEvents_Lagrangian'
    [R] = originalData(calEvent,'double');
    [RTime] = getTime(calEvent);
    % remove those events with missing data
    if any(isnan(R(:)))
        % # to be confirmed #
        % correct NaN val into ...
        if nanmean(isnan(R(:)))<0.1
            R(isnan(R(:))) = rand(nansum(isnan(R(:))),1)*1*nanmean(R(:));
            evi = evi+1;
            Rains{evi,1} = R;
            Rmon(evi,1) = RTime(1).Month;
        end
    else
        evi = evi+1;
        Rains{evi,1} = R;
        Rmon(evi,1) = RTime(1).Month;
    end
end

%% Estimate ARMA coef. (p [1,4], q [0,4])
[ ARMA, mean_spatial_correlation ] = ARMA_estimate( Rains, 1, true );
pause(1)
savePlot(['Figs\Lagrangian\event',num2str(evi),'-ARMAfit'],...
    'wholepage',true,'onlyPng',true,'needreply','N');
close(gcf);
save('Birmingham\ARMA.mat','ARMA','mean_spatial_correlation','-v7.3')
save('Birmingham\ARMA_Rains.mat','Rains','-v7.3')
%% Plot Fitting result
figure;
colm = flip(pink(14),1);

for evi = 1:size(mean_spatial_correlation,1)
    autoCorr = mean_spatial_correlation(evi,:);
    lags = 0:length(autoCorr)-1;
    plot(lags,autoCorr,'color',colm(Rmon(evi),:));
    hold on;
    drawnow
end
xlabel('Lag');
ylabel('acf in lagrangian coor');
savePlot(['Figs\Lagrangian\TemporalEvol'],'onlyPng',true,'needreply','N');

% for i=0:4
%     for j=0:4
%         Mdl=arima(i,0,j);
%         [~,~,LLF] = estimate(Mdl,G);
%         NumActive = i+j+1;
%         AIC(i+1,j+1) = aicbic(LLF,NumActive);
%     end
% end
% 
% % Choose the model to use
% [i,j]=find(min(min(AIC))==AIC);
% Mdl=arima(i-1,0,j-1);


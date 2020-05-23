function [ thresholdI ] = I_Threshold_estimation( )
%%
%IMF_Threshold_estimation This is an example for estimating the IMF
%threshold needed to simulate propely the number of days exceeding 1 mm of rainfall.
%% Load data
load('CV.mat') % Rainfall coefficient of variation [-], monthly data
load('ARMA.mat') % ARMA coefficients
load('expoSpatialCorrelation.mat') % Exponential coefficient of the spatial correlation
load('NHRO.mat') % Non homogenic rainfall ocurrence; in this example nu,ber of days above a certain threshold, can be also hours or minutes
load('Precipitation.mat');
Precipitation.data=data;
Precipitation.WAR=WAR;
Precipitation.IMF=IMF;
load('triNormTrans.mat')
load('Cloud.mat')
Cloud=load('Cloud.mat','ar');
Cloud.cloudDis=cloudDis;
Cloud.Dry=Dry;
Cloud.Wet=Wet;
load('obsOCC.mat')
%% Running a 30 years simulation for the domain
PreProc_30y_run_imf
%% Setting the objective target - Number of days exceeding 1 mm daily rainfall
objFun=[];
for m=1:12
    objFun(m,:)=nanmean(NHRO(m).occurence);
end
objFun=reshape(sum(objFun),13,13);
%% Estimating the IMF threshold
thresholdI=zeros(13,13);
for idx=1:13*13
    fun = @(x) minFun(x,RFIELD,objFun,idx);
    thresholdI(idx)=fminbnd(fun,0,3); % Search between threshold of 0 to 3 [mm]
end
%% Plot - comparing observed and simulated wet days (above 1 mm threshold)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(objFun)
axis('square')
colorbar
caxis([145 166])
title('Averaged rainfall days > 1 mm [Observed]')
[ simNdays ] = countDays2( thresholdI , RFIELD );
subplot(1,2,2)
imagesc(simNdays)
axis('square')
colorbar
caxis([145 166])
title('Averaged rainfall days > 1 mm [Simulated]')
savefig('Averaged_rainfall_days_1_mm.fig');
end
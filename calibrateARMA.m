%
% calibrate ARMA coef.
%

% method:
% all time series of variables are kept same as observed.
% only advection are randomly sampled.
clear;clc
close all

dx = 1;%2
dt = 5;
domainSize = 110;%13
WAR_threshold = 0.02;%0.1;
YEAR = 2007;
stormi = 1:10;
% ARMA_test = struct('ar',0.98,'ma',[]);
load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham\ARMA.mat')
ARMA_test = ARMA(1);

times = tic;
rng('shuffle');
addpath(genpath(cd))

load('Birmingham\MeanArealStats.mat')
load('Files\radarIsEvent.mat','radarIsEvent');
load('Files\obsOCC.mat','obsOCC');
load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham\CV.mat')
Precipitation = load('Files\Precipitation.mat');

% observed Storm arrival process
tsWetDry = radarIsEvent(:,1)>=1;
tsMonth = datetime(datevec(radarIsEvent(:,2))).Month;
tsYear = datetime(datevec(radarIsEvent(:,2))).Year;
[ tsMatrix ] = tsTable( double(tsWetDry(tsYear == YEAR)') , tsMonth(tsYear == YEAR)' );
tsMatrix(tsMatrix(:,1)==0,:) = [];% only simulate wet period to save computational time
tsMatrix = tsMatrix(stormi,:);
% observed WAR, IMF
obsWAR = MeanArealStats.WAR(tsYear==YEAR,1);
obsIMF = MeanArealStats.IMF(tsYear==YEAR,1);
obsIMF = fillmissing(obsIMF, 'previous');
obsWAR = fillmissing(obsWAR, 'previous');
% [simU,simV] = getSimUV(tsMatrix);
[simU,simV] = getObsUV(tsMatrix,MeanArealStats);

%% Generating rain fields %(only for a certain amount of storms)
simRain=NaN(105120,domainSize,domainSize);
QField = [];
progressbar('Generating rain fields') % Init single bar
for m=1:12
    tic
    if ~isempty(min(tsMatrix(tsMatrix(:,2)==m,3)))
    mmin(m,1)=min(tsMatrix(tsMatrix(:,2)==m,3));
    mmax(m,1)=max(tsMatrix(tsMatrix(:,2)==m,4));
    [ QField(mmin(m,1):mmax(m,1),:,:) ] = quantileFieldGen( [domainSize*2 domainSize*2] , ARMA_test , 5 , ...
        Precipitation.data(m).SpatialAlpha , [mmin(m,1) mmax(m,1)] , simU , simV , dx , dt );
    %% Non-homogeneous probability of precipitation occurrence
    for i=mmin(m,1):mmax(m,1)
        [ simRain(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , obsWAR(i) , obsOCC{m} );
    end
    %% Assigning rainfall intensity for the Gaussian field
    for i=1:size(tsMatrix,1)
        M=tsMatrix(i,2);
        if tsMatrix(i,1)==1 && m==M
            [ simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) ] = invLN2( simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) , ...
                obsWAR(tsMatrix(i,3):tsMatrix(i,4)) , obsIMF(tsMatrix(i,3):tsMatrix(i,4)) , CV(m));
        end
    end
    end
    progressbar(m/12) % Update progress bar
    toc
end
timee = toc(times);
%%
fprintf('Total time required: %3.1f seconds\n',timee);
load(['H:\CODE_MATLAB\',sprintf('PRS_CLEEHILL_%04d.mat',YEAR)],'DATA');
%%
simStorm = [];
h = figure;
for i=1:size(tsMatrix,1)
    stormInd = tsMatrix(i,3):tsMatrix(i,4);
    simStorm{i} = squeeze(simRain(stormInd,:,:));
    simStorm{i} = reshape(simStorm{i},size(simStorm{i},1),[]);
    simStorm{i} = fillmissing(simStorm{i}, 'previous');
    obsStorm{i} = permute(squeeze(originalData(extractOnePeriod(DATA,stormInd))),[3,1,2]);
    obsStorm{i} = reshape(obsStorm{i},size(simStorm{i},1),[]);
    subplot(2,1,1)
    plotACF(obsStorm{i});hold on;plotACF(simStorm{i});
    title('Autocorrelation');
    xlabel('Lags');ylabel('Acorr')
    ylim([-0.1,1])
    subplot(2,1,2)
    plot(obsIMF(tsMatrix(i,3):tsMatrix(i,4)));hold on;
    plot(nanmean(reshape(simStorm{i},size(simStorm{i},1),[]),2))
    title('Areal Mean Precipitation');
    xlabel('Time steps[5min]');ylabel('IMF[mm/h]')
    clf(h)
end



function acf = plotACF(rain)
acf = [];
for sitei = 1:size(rain,2)
    [acf(sitei,:),lags,bounds] = autocorr(rain(:,sitei),20);
end
plot(lags,nanmean(acf,1));
end


%% AUXILLARY function
function [obsU,obsV] = getObsUV(tsMatrix,MeanArealStats)
% obsU = MeanArealStats.U500(tsYear==YEAR,1);
% obsV = MeanArealStats.V500(tsYear==YEAR,1);
% obsU = fillmissing(obsU, 'previous');
% obsV = fillmissing(obsV, 'previous');
for i=1:size(tsMatrix,1)
    obsU(tsMatrix(i,3):tsMatrix(i,4),1)=MeanArealStats.U500(tsMatrix(i,3):tsMatrix(i,4));
    obsV(tsMatrix(i,3):tsMatrix(i,4),1)=MeanArealStats.V500(tsMatrix(i,3):tsMatrix(i,4));
end
end
function [simU,simV] = getSimUV(tsMatrix)
Advection = load('Files\Advection.mat');
%% Simulating advection
tic
[Sd,Sw]=deal(cell(1,size(tsMatrix,1)));
simU=NaN(tsMatrix(end),1);
simV=NaN(tsMatrix(end),1);
if size(tsMatrix,1)>10
    parfor i=1:size(tsMatrix,1)
        % Sd{i}=simulate(Advection.Dry(tsMatrix(i,2)).mdl,20000);
        Sw{i}=simulate(Advection.Wet(tsMatrix(i,2)).mdl,11000);
    end
else
    for i=1:size(tsMatrix,1)
        % Sd{i}=simulate(Advection.Dry(tsMatrix(i,2)).mdl,20000);
        Sw{i}=simulate(Advection.Wet(tsMatrix(i,2)).mdl,11000);
    end
end
toc

progressbar('Generating advection components') % Init single bar
for i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if i==1
        switch tsMatrix(i,1)
            case 0
                %
            case 1
                Ss=simulate(Advection.Wet(tsMatrix(i,2)).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1);
                YS=Ss(end,:);
                Ss=normcdf(Ss);
                Ss=round(Ss.*1000)./1000;
                Ss(Ss==1)=0.999;
                Ss(Ss==0)=0.001;
                Ss(:,1)=norminv(Ss(:,1),Advection.Wet(m).unorm(1),Advection.Wet(m).unorm(2));
                Ss(:,2)=norminv(Ss(:,2),Advection.Wet(m).vnorm(1),Advection.Wet(m).vnorm(2));
        end
        simU(tsMatrix(i,3):tsMatrix(i,4),1)=Ss(:,1);
        simV(tsMatrix(i,3):tsMatrix(i,4),1)=Ss(:,2);
    end
    if i>1
        switch tsMatrix(i,1)
            case 0
                %
            case 1
                E=find(min((Sw{i}(1:10000,1)-YS(1)).^2+(Sw{i}(1:10000,2)-YS(2)).^2)==(Sw{i}(1:10000,1)-YS(1)).^2+(Sw{i}(1:10000,2)-YS(2)).^2);
                E=E(1);
                S=Sw{i}(E:E+tsMatrix(i,4)-tsMatrix(i,3),:);
                YS=S(end,:);
                S=normcdf(S);
                S(:,1)=norminv(S(:,1),Advection.Wet(m).unorm(1),Advection.Wet(m).unorm(2));
                S(:,2)=norminv(S(:,2),Advection.Wet(m).vnorm(1),Advection.Wet(m).vnorm(2));
        end
        simU(tsMatrix(i,3):tsMatrix(i,4),1)=S(:,1);
        simV(tsMatrix(i,3):tsMatrix(i,4),1)=S(:,2);
    end
    progressbar(i/size(tsMatrix,1)) % Update progress bar
end
% simU=single(simU);
% simV=single(simV);
clear S Ss Sd Sw
end


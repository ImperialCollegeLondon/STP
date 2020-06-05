% STREAP Main code
tic
%% Initilazing
rng('shuffle');
addpath(genpath(cd))
%% PreProcessing
%% Load data
filePath = 'Birmingham/';%'Engelberger/';%
load([filePath,'MeanArealStats.mat']) % WAR [-], (CAR [-]) and IMF [mm h^-1] data
% U500, V500: [m s^-1]
% NO MeanArealStats.CAR;
load([filePath,'CV.mat']) % Rainfall coefficient of variation [-], monthly data
load([filePath,'ARMA.mat']) % ARMA coefficients
%
% load([filePath,'Gauge_10_min_precip.mat'], 'MERGED') % Boolean rain / no-rain (10-min interval), composite of 4 MeteoSwiss rain-gauge data
%
load([filePath,'expoSpatialCorrelation.mat']) % Exponential coefficient of the spatial correlation
load([filePath,'NHRO.mat']) % Non homogenic rainfall ocurrence; in this example number of days above a certain threshold, can be also hours or minutes
%
% load([filePath,'gpFitSO.mat']) % Generalized Pareto parameters for the observed and simulated daily data
% load([filePath,'thresholdI.mat']) % Simulated daily rain intensity threshold [mm] (from which daily rainfall > 1 mm)
%
dx = 1;%2
dt = 5;
domainSize = 110;%13
WAR_threshold = 0.1;%0.02;
%% Rainfall event analysis - gauged based
% isEventdt = 10;%10min
% [ gaugeIsEvent ] = findGaugeEvents( MERGED(:,1) , 0.1 , isEventdt ); 
% gaugeIsEvent(:,2) = MERGED(:,2);
% save('Files\gaugeIsEvent.mat','gaugeIsEvent','-v7.3');
% Boolean vector defining an event for each 10-min; in this case, 
% there is no meaning for the rainfall threshold (set arbitrary to 0.1 [mm]). 
% Similar file exists for radar analysis: find RadarEvents
[ radarIsEvent ] = findRadarEvents( MeanArealStats.WAR , WAR_threshold , 60 );
radarIsEvent(:,2) = MeanArealStats.Time;
save('Files\radarIsEvent.mat','radarIsEvent','-v7.3');
isEventdt = 5;
%% Storm arrival process
% Fitting distributions for the wet and dry periods
for m=1:12
    Precipitation.data(m).isWetEventG=radarIsEvent(month(radarIsEvent(:,2))==m);
    Precipitation.data(m).isWetEventGyear=year(radarIsEvent(month(radarIsEvent(:,2))==m,2));
    Precipitation.data(m).isWetEventG(Precipitation.data(m).isWetEventG>1)=1;
    [ Precipitation.data(m).dryG , Precipitation.data(m).wetG , ...
        Precipitation.data(m).dryFitG , Precipitation.data(m).wetFitG ] = ...
        WetDryG( Precipitation.data(m).isWetEventG , Precipitation.data(m).isWetEventGyear , 1 , isEventdt );
    if m==1
        export_fig Figs\tif\Wet_Dry_Fit_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\WetDryFitG',num2str(m),'.fig']); close(gcf);
end

%% Rainfall spatial correlation
% Finding alpha parameter for the FFT filter that match the spatial correlation exponential
% coefficient
% 13*2-1 % 110*2-1
dim = domainSize*2-1;%25;
parfor m=1:12
    [ tmp{m} ] = fftAutoCorr( [dim dim] , [1 40] , dx , expoSpatialCorr(m));
end
for m=1:12
    Precipitation.data(m).SpatialAlpha=tmp{m}(2);
end
clear tmp

%% Fit IMF distribution
IMF=MeanArealStats.IMF(MeanArealStats.IMF>0&MeanArealStats.WAR>=WAR_threshold); 
% Lower threshold of WAR>=0.1 and IMF>0 is set for the IMF statistics
T=MeanArealStats.Time(MeanArealStats.IMF>0&MeanArealStats.WAR>=WAR_threshold);
for m=1:12
    Precipitation.IMF(m).IMF=IMF(month(T)==m);
    [ ~ , ~ , ~ , ~ , MGGP , ~ , ~ ] = fitIMF( Precipitation.IMF(m).IMF(:) , 1 ); % For this example the Mixed Gamma and aGeneralized Pareto distribution is applied; other distributions can be applied using this code (Gamma, Generalized Pareto , Mixed Exponential and LogNormal
    Precipitation.IMF(m).case=5; % Out of 5 possible choices
    Precipitation.IMF(m).name='Mixed Gamma and Generalized Pareto';
    Precipitation.IMF(m).fit=MGGP;
    if m==1
        export_fig Figs\tif\Fit_IMF_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\fitIMF',num2str(m),'.fig']); close(gcf);
end
clear MGGP IMF

%% Fit WAR distribution
WAR=MeanArealStats.WAR(MeanArealStats.IMF>0&MeanArealStats.WAR>=WAR_threshold); 
% Lower threshold of WAR>=0.1 and IMF>0 is set for the IMF statistics
for m=1:12
    Precipitation.WAR(m).WAR=WAR(month(T)==m);
    [ ~ , ~ , John , ~ , ~] = fitWAR( Precipitation.WAR(m).WAR , 1 ); % For this example the Johnson distribution is applied; Beta or LogNormal distributions can be also applied using this code
    Precipitation.WAR(m).case=3;
    Precipitation.WAR(m).name='Johnson';
    Precipitation.WAR(m).fit=John;
    if m==1
        export_fig Figs\tif\Fit_WAR_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\fitWAR',num2str(m),'.fig']); close(gcf);
end
clear WAR John

%% Find temporal correlation between IMF and WAR
for m=1:12
    Precipitation.data(m).rho=corr([Precipitation.IMF(m).IMF Precipitation.WAR(m).WAR]);
end
save('Files\Precipitation','-struct','Precipitation','-v7.3')

IsEvent=zeros(size(MeanArealStats.IMF));
IsEvent(MeanArealStats.IMF>0&MeanArealStats.WAR>=WAR_threshold)=1; % Same thresholds as for the WAR and CAR assigned before

%% Estimating Varfima parameters
[ IsEventID ] = isEventID( IsEvent );
[ Norm , Copula , MaternCov ] = biNormTransform( IsEventID , MeanArealStats.IMF , MeanArealStats.WAR , MeanArealStats.Time , Precipitation);
save('Files\biNormTrans.mat','Norm','Copula','MaternCov','-v7.3')

%% VAR parameter estimation - Advection
for m=1:12
    % wet:
    U=MeanArealStats.U500(IsEvent==1&month(MeanArealStats.Time)==m);
    V=MeanArealStats.V500(IsEvent==1&month(MeanArealStats.Time)==m);
    [ Advection.Wet(m).mdl , Advection.Wet(m).unorm , Advection.Wet(m).vnorm , h ] = UVadvectionVAR( U , V , 12 , true );
    Advection.Wet(m).stdU=std(U);
    Advection.Wet(m).stdV=std(V);
    if m==1
        figure(h(2))
        export_fig Figs\tif\Fit_500hPa_wind_wet_example.tif -tif -transparent -r300
    end
    savefig(h(1),['Figs\Matlab\PACFWet',num2str(m),'.fig']);
    savefig(h(2),['Figs\Matlab\uvDistWet',num2str(m),'.fig']);
    close all
    % dry:
    U=MeanArealStats.U500(IsEvent==0&month(MeanArealStats.Time)==m);
    V=MeanArealStats.V500(IsEvent==0&month(MeanArealStats.Time)==m);
    [ Advection.Dry(m).mdl , Advection.Dry(m).unorm , Advection.Dry(m).vnorm , h ] = UVadvectionVAR( U , V , 12 , true );
    Advection.Dry(m).stdU=std(U);
    Advection.Dry(m).stdV=std(V);
    if m==1
        figure(h(2))
        export_fig Figs\tif\Fit_500hPa_wind_dry_example.tif -tif -transparent -r300
    end
    savefig(h(1),['Figs\Matlab\PACFDry',num2str(m),'.fig']);
    savefig(h(2),['Figs\Matlab\uvDistDry',num2str(m),'.fig']);
    close all
end
clear U V
save('Files\Advection.mat','-struct','Advection','-v7.3')

%% Estimating the non-homogeneous probability of precipitation occurrence fields (normalized)
% This is unknown parameter for the Engerberger catchment as the weather radar data
% lack the quality for a good estimation of the 5-min non-homogeneous
% probability of precipitation occurrence. The non-homogeneous probability
% of precipitation occurrence is therefore assume to follow the annual
% rainfall filter. ALTERNATIVLY - one can choose to apply the NHRO daily
% filter (as presented in AWE-2d paper)
load([filePath,'obsAnnualRain.mat'])%13*13 grids
ri=zeros(domainSize,domainSize);
for y=1:length(ObsAnnualRain)
    ri=ri+ObsAnnualRain{y};
end
ri=ri./y;
tmp=reshape(ri./mean2(ri),[1 domainSize*domainSize]);
for m=1:12
    obsOCC{m}=tmp;
end
clear tmp y m
save('Files\obsOCC.mat','obsOCC','-v7.3');

%% Processing
load('Files\radarIsEvent.mat','radarIsEvent');
load('Files\obsOCC.mat','obsOCC');
Advection = load('Files\Advection.mat');
load('Files\biNormTrans.mat','Norm','Copula','MaternCov')
Precipitation = load('Files\Precipitation.mat');
load([filePath, 'CV.mat'])
I = 1; % Number of realizations to simulate (the processing stage can be looped)
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
if isEventdt == 5    
    [ tsMonth , tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , dt/(24*60) );
else
    % for 10min gauge data % deleted by Yt
    [ tsMonth , tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , isEventdt/(24*60) );
    for i=1:length(tsWetDry) % For this example the [dt]-min temporal resolution is interpolated from the [dt*2]-min resolution
        tsWetDry5Min(i*2-1)=tsWetDry(i);
        tsWetDry5Min(i*2)=tsWetDry(i);
        tsMonth5Min(i*2-1)=tsMonth(i);
        tsMonth5Min(i*2)=tsMonth(i);
    end
    tsWetDry=tsWetDry5Min;
    tsMonth=tsMonth5Min;
    clear tsMonth5Min tsWetDry5Min
end
[ tsMatrix ] = tsTable( tsWetDry , tsMonth );
%% Varfima - simulating WAR, IMF
simWAR=zeros(tsMatrix(end),1);
simIMF=zeros(tsMatrix(end),1);
simCAR=zeros(tsMatrix(end),1);
parfor i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if tsMatrix(i,1)==1
        [ tsWARFIMA{i} ] = MultiBiVARFIMA( MaternCov , Precipitation , Copula , m , tsMatrix(i,4)-tsMatrix(i,3)+1 , 5 );
    end
end

for i=1:size(tsMatrix,1)
    if tsMatrix(i,1)==1
        E=abs(0-tsWARFIMA{i}.WARreal(1,:))+abs(0-tsWARFIMA{i}.WARreal(end,:));
        E=find(min(E)==E);
        if length(E)>1
            E=E(1);
        end
        simWAR(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.WARreal(:,E);
        simIMF(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.IMFreal(:,E);
        % simCAR(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.CARreal(:,E);
    end
end
clear tsWARFIMA m
%% Simulating advection
simU=zeros(tsMatrix(end),1);
simV=zeros(tsMatrix(end),1);
parfor i=1:size(tsMatrix,1)
    Sd{i}=simulate(Advection.Dry(tsMatrix(i,2)).mdl,20000);
    Sw{i}=simulate(Advection.Wet(tsMatrix(i,2)).mdl,20000);
end
poolobj = gcp('nocreate');
delete(poolobj);

progressbar('Generating advection components') % Init single bar
for i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if i==1
        switch tsMatrix(i,1)
            case 0
                Ss=simulate(Advection.Dry(tsMatrix(i,2)).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1);
                YS=Ss(end,:);
                Ss=normcdf(Ss);
                Ss=round(Ss.*1000)./1000;
                Ss(Ss==1)=0.999;
                Ss(Ss==0)=0.001;
                Ss(:,1)=norminv(Ss(:,1),Advection.Dry(m).unorm(1),Advection.Dry(m).unorm(2));
                Ss(:,2)=norminv(Ss(:,2),Advection.Dry(m).vnorm(1),Advection.Dry(m).vnorm(2));
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
                E=find(min((Sd{i}(1:10000,1)-YS(1)).^2+(Sd{i}(1:10000,2)-YS(2)).^2)==(Sd{i}(1:10000,1)-YS(1)).^2+(Sd{i}(1:10000,2)-YS(2)).^2);
                E=E(1);
                S=Sd{i}(E:E+tsMatrix(i,4)-tsMatrix(i,3),:);
                YS=S(end,:);
                S=normcdf(S);
                S(:,1)=norminv(S(:,1),Advection.Dry(m).unorm(1),Advection.Dry(m).unorm(2));
                S(:,2)=norminv(S(:,2),Advection.Dry(m).vnorm(1),Advection.Dry(m).vnorm(2));
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
simU=single(simU);
simV=single(simV);
clear S Ss Sd Sw
%% Generating rain fields
ARMA22 = ARMA(2);%ARMA(8);
% ARMA22 = ARMA;
simRain=zeros(105120,domainSize,domainSize);
QField = [];
progressbar('Generating rain fields') % Init single bar
tic
for m=1:12
    tic
    mmin(m,1)=min(tsMatrix(tsMatrix(:,2)==m,3));
    mmax(m,1)=max(tsMatrix(tsMatrix(:,2)==m,4));
    [ QField(mmin(m,1):mmax(m,1),:,:) ] = quantileFieldGen( [domainSize*2 domainSize*2] , ARMA22 , 5 , ...
        Precipitation.data(m).SpatialAlpha , [mmin(m,1) mmax(m,1)] , simU , simV , dx , dt );
    % Non-homogeneous probability of precipitation occurrence
    % for i=mmin(m,1):mmax(m,1)
    %     [ simRain(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simWAR(i) , obsOCC{m});%%%obsOCC{m} );
    % end
    % Homogeneous probability of precipitation occurrence
    for i=mmin(m,1):mmax(m,1)
         [ simRain(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simWAR(i) , obsOCC{m}*0+1);%%%obsOCC{m} );
    end
    %% Assigning rainfall intensity for the Gaussian field
    for i=1:size(tsMatrix,1)
        M=tsMatrix(i,2);
        if tsMatrix(i,1)==1 && m==M
            [ simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) ] = invLN2( simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) , ...
                simWAR(tsMatrix(i,3):tsMatrix(i,4)) , simIMF(tsMatrix(i,3):tsMatrix(i,4)) , CV(m));
        end
    end
    toc
    progressbar(m/12) % Update progress bar
end
simRain(isnan(simRain))=0;
simRain=single(simRain);
toc
% %% Non-homogeneous spatial rainfall accumulation
% % Simulation is divided to daily time steps to match the correction records
% t=(datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,23,55,0))'; % One year dates based on the year 2001
% tmp=squeeze(sum(reshape(simRain,[24*60/dt 365 domainSize domainSize]))).*dt/60;
% for i=1:domainSize*domainSize
%     for d=1:365
%         if tmp(d,i)<thresholdI(i)
%             k=tmp(d,i)/thresholdI(i);
%             simRain(1+24*60/dt*(d-1):24*60/dt*d,i)=simRain(1+24*60/dt*(d-1):24*60/dt*d,i).*k; % Below the 1 mm threshold is consider as "drizzle"
%             tmp(d,i)=0;
%         end
%     end
% end
% for d=1:365
%     m=month(t(d));
%     tmp2=squeeze(tmp(d,:,:));
%     gpCDF=gpcdf(tmp2,gpFitS{1}{m},gpFitS{2}{m},thresholdI);
%     gpCDF(gpCDF>0.999)=0.999;
%     gpINV=gpinv(gpCDF,gpFitO{1}{m},gpFitO{2}{m},1);
%     k=gpINV./tmp2; % Daily rainfall correction factor
%     k(isinf(k))=1;
%     simRain(1+24*60/dt*(d-1):24*60/dt*d,:,:)=simRain(...
%         1+24*60/dt*(d-1):24*60/dt*d,:,:).*repmat(reshape(k,[1 domainSize domainSize]),24*60/dt,1,1);
% end
% simRain(isnan(simRain))=0;
% simRain=single(simRain);

%% Saving output
save('Files\Output_Pr_Clt_500UV.mat','simRain','simIMF','simWAR','simU','simV','-v7.3'); % Matlab
% clear simRain

toc

pause(10);
run('H:\CODE_MATLAB\AWE-GEN-2D\Video_Map_STREAP.m')


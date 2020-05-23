% AWE-GEB-2D Main code
%% Initilazing
rng('shuffle');
addpath(genpath(cd))
%% PreProcessing
%% Load data
load('MeanArealStats.mat') % WAR [-], CAR [-] and IMF [mm h^-1] data
load('CV.mat') % Rainfall coefficient of variation [-], monthly data
load('ARMA.mat') % ARMA coefficients
load('Gauge_10_min_precip.mat', 'MERGED') % Boolean rain / no-rain (10-min interval), composite of 4 MeteoSwiss rain-gauge data
load('expoSpatialCorrelation.mat') % Exponential coefficient of the spatial correlation
load('NHRO.mat') % Non homogenic rainfall ocurrence; in this example number of days above a certain threshold, can be also hours or minutes
load('gpFitSO.mat') % Generalized Pareto parameters for the observed and simulated daily data
load('thresholdI.mat') % Simulated daily rain intensity threshold [mm] (from which daily rainfall > 1 mm)
load('LRdata.mat') % Temperature data for two stations to calculate lapse rate
load('CloudCover.mat') % Hourly cloud cover data for the temperature analysis
load('Ta_Stations.mat') % Temperature stations located within the domain
load('Rsw_Stations.mat') % Short wave radiation from stations
load('Ps_Stations.mat') % Pressure at ground level from stations
load('AOD.mat') % AERONET aerosol optical depth data
domainDTM=rasterread('engelberger_DTM.txt'); % 250 x 250 pixels [100 x 100 m] DTM of the Engelberger catchment
domainDTM=flipud(domainDTM);
%% Rainfall event analysis - gauged based
[ gaugeIsEvent ] = findGaugeEvents( MERGED(:,1) , 0.1 , 10 ); % Boolean vector defining an event for each 10-min; in this case, there is no meaning for the rainfall threshold (set arbitrary to 0.1 [mm]). Similar file exists for radar analysis: findRadarEvents
gaugeIsEvent(:,2)=MERGED(:,2);
save('Files\gaugeIsEvent.mat','gaugeIsEvent','-v7.3');
%% Storm arrival process
% Fitting distributions for the wet and dry periods
for m=1:12
    Precipitation.data(m).isWetEventG=gaugeIsEvent(month(gaugeIsEvent(:,2))==m);
    Precipitation.data(m).isWetEventGyear=year(gaugeIsEvent(month(gaugeIsEvent(:,2))==m,2));
    Precipitation.data(m).isWetEventG(Precipitation.data(m).isWetEventG>1)=1;
    [ Precipitation.data(m).dryG , Precipitation.data(m).wetG , Precipitation.data(m).dryFitG , Precipitation.data(m).wetFitG ] = WetDryG( Precipitation.data(m).isWetEventG , Precipitation.data(m).isWetEventGyear , 1 , 10 );
    if m==1
        export_fig Figs\tif\Wet_Dry_Fit_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\WetDryFitG',num2str(m),'.fig']); close(gcf);
end
%% Rainfall spatial correlation
% Finding alpha parameter for the FFT filter that match the spatial correlation exponential
% coefficient
parfor m=1:12
    [ tmp{m} ] = fftAutoCorr( [25 25] , [1 20] , 2 , expoSpatialCorr(m));
end
for m=1:12
    Precipitation.data(m).SpatialAlpha=tmp{m}(2);
end
clear tmp
%% Fit IMF distribution
IMF=MeanArealStats.IMF(MeanArealStats.IMF>0&MeanArealStats.WAR>=0.1); % Lower threshold of WAR>=0.1 and IMF>0 is set for the IMF statistics
T=MeanArealStats.Time(MeanArealStats.IMF>0&MeanArealStats.WAR>=0.1);
for m=1:12
    Precipitation.IMF(m).IMF=IMF(month(T)==m);
    [ ~ , ~ , ~ , ~ , MGGP , ~ , ~ ] = fitIMF( Precipitation.IMF(m).IMF , 1 ); % For this example the Mixed Gamma and aGeneralized Pareto distribution is applied; other distributions can be applied using this code (Gamma, Generalized Pareto , Mixed Exponential and LogNormal
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
WAR=MeanArealStats.WAR(MeanArealStats.IMF>0&MeanArealStats.WAR>=0.1); % Lower threshold of WAR>=0.1 and IMF>0 is set for the IMF statistics
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
%% Fit CAR distribution
% Seperate distributions for wet and dry periods
CAR=MeanArealStats.CAR(MeanArealStats.IMF>0&MeanArealStats.WAR>=0.1); % Lower threshold of WAR>=0.1 and IMF>0 is set for the IMF statistics
CARd=MeanArealStats.CAR(MeanArealStats.IMF==0|MeanArealStats.WAR<0.1); % Data for dry period
Td=MeanArealStats.Time(MeanArealStats.IMF==0|MeanArealStats.WAR<0.1);
for m=1:12
    Cloud.Wet(m).CAR=CAR(month(T)==m);
    Cloud.Dry(m).CAR=CARd(month(Td)==m);
    [ ~ , ~ , ~ , John , ~ , aicc ] = fitCAR( Cloud.Wet(m).CAR , 1 );
    Cloud.Wet(m).case=2; % For this example the Johnson distribution is applied; LogNormal, Normal and Weibull are also available
    Cloud.Wet(m).name='Johnson';
    Cloud.Wet(m).fit=John;
    if m==1
        export_fig Figs\tif\Fit_CAR_wet_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\fitCARw',num2str(m),'.fig']); close(gcf);
    [ ~ , ~ , ~ , John , ~ , aicc ] = fitCAR( Cloud.Dry(m).CAR , 1 );
    Cloud.Dry(m).case=2;
    Cloud.Dry(m).name='Johnson';
    Cloud.Dry(m).fit=John;
    if m==1
        export_fig Figs\tif\Fit_CAR_dry_example.tif -tif -transparent -r300
    end
    savefig(gcf,['Figs\Matlab\fitCARd',num2str(m),'.fig']); close(gcf);
end
clear CAR CARd T Td John
%% Find temporal correlation between CAR, IMF and WAR
for m=1:12
    Precipitation.data(m).rho=corr([Cloud.Wet(m).CAR Precipitation.IMF(m).IMF Precipitation.WAR(m).WAR]);
end
save('Files\Precipitation','-struct','Precipitation','-v7.3')
%% CAR autoregressive estimation for inter-storm periods
IsEvent=zeros(size(MeanArealStats.CAR));
IsEvent(MeanArealStats.IMF>0&MeanArealStats.WAR>=0.1)=1; % Same thresholds as for the WAR and CAR assigned before
[ Cloud.cloudDis ] = WetEventDis( MeanArealStats.CAR , IsEvent , 5 );
progressbar('CAR autoregressive estimation') % Init single bar
for m=1:12
    Cloud.ar(m).tThreshold=max(Cloud.cloudDis(month(MeanArealStats.Time)==m));
    [ Cloud.ar(m).disMiu , Cloud.ar(m).disStd ] = CloudWetEventDis( MeanArealStats.CAR(month(MeanArealStats.Time)==m) , Cloud.cloudDis(month(MeanArealStats.Time)==m) , 5 , Cloud.ar(m).tThreshold , 0 );
    [ Cloud.ar(m).mdl ] = cloudAR( MeanArealStats.CAR(month(MeanArealStats.Time)==m) , Cloud.cloudDis(month(MeanArealStats.Time)==m) , Cloud.ar(m).disMiu , Cloud.ar(m).disStd , Cloud.ar(m).tThreshold , Cloud.Dry(m).fit , 60 , 0 );
    progressbar(m/12) % Update progress bar
end
save('Files\Cloud','-struct','Cloud','-v7.3')
%% Estimating Varfima parameters
[ IsEventID ] = isEventID( IsEvent );
[ Norm , Copula , MaternCov ] = triNormTransform( IsEventID , MeanArealStats.IMF , MeanArealStats.WAR , MeanArealStats.CAR , MeanArealStats.Time , Precipitation , Cloud );
save('Files\triNormTrans.mat','Norm','Copula','MaternCov','-v7.3')
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
load('obsAnnualRain.mat')
ri=zeros(13,13);
for y=1:32
    ri=ri+ObsAnnualRain{y};
end
ri=ri./y;
tmp=reshape(ri./mean2(ri),[1 13*13]);
for m=1:12
    obsOCC{m}=tmp;
end
clear tmp y m
save('Files\obsOCC.mat','obsOCC','-v7.3');
%% Computing the flow accomulation and topographic curves - Temperature
plotDTM(domainDTM) % Plot the domain elevation [m], grid size: 100 m x 100 m
export_fig Figs\tif\DTM.tif -tif -transparent -r300
savefig(gcf,'Figs\Matlab\Engelberger_DTM.fig'); close(gcf);
% Accumulation:
AccThreshold=500;
[ AccMatrix ] = catchmentAccumulation( domainDTM , AccThreshold );
figure(1)
export_fig Figs\tif\FlowAccIdx.tif -tif -transparent -r300
figure(2)
export_fig Figs\tif\FlowAccDTM.tif -tif -transparent -r300
savefig(gcf,'Figs\Matlab\FlowAccDTM.fig'); close(gcf);
savefig(gcf,'Figs\Matlab\FlowAccIdx.fig'); close(gcf);
% Curves:
[ ~ , planc ] = catchmentCurv( domainDTM );
export_fig Figs\tif\Curvature.tif -tif -transparent -r300
savefig(gcf,'Figs\Matlab\Curvature.fig'); close(gcf);
% Combine curvature and flow accumulation indexes to define the elevation
% field threshold for the temperature inversion
% Compute and plot Z1 and Z2
CurvThreshold=[-10 10];
CurveAccDist=3;
[ Z1 ] = defineZ1( domainDTM , AccMatrix , AccThreshold , planc , CurvThreshold , CurveAccDist );
export_fig Figs\tif\Z1Channels.tif -tif -transparent -r300
savefig(gcf,'Figs\Matlab\Z1Channels.fig'); close(gcf);
Z2=domainDTM-Z1;
plotZ1Z2( Z1 , Z2 )
figure(1)
export_fig Figs\tif\Z1.tif -tif -transparent -r300
figure(2)
export_fig Figs\tif\Z2.tif -tif -transparent -r300
savefig(gcf,'Figs\Matlab\Z2.fig'); close(gcf);
savefig(gcf,'Figs\Matlab\Z1.fig'); close(gcf);
save('Files\Zlevels.mat','Z1','Z2');
%% Temperature lapse rate inversion
% Estimating hourly lapse rate from mountain and valley stations to account
% for the thermal inversion. A copula is calculated between the daily
% temperature amplitude and the lapse rate
StationEle(1)=454; % Station elevation a.s.l. [m]
StationEle(2)=2106;
progressbar('Calculate monthly lapse rate including inversion') % Init single bar
for m=1:12
    [ LR(m).data , LR(m).date , LR(m).T0 , LR(m).T0amp , LR(m).Th1 , LR(m).Th2 ] = computeLR( LRdata , StationEle(1) , StationEle(2) , m );
    progressbar(m/12) % Update progress bar
end
for m=1:12
    H=hour(LR(m).date);
    for h=0:23
        X=LR(m).T0amp(H==h);
        Y=LR(m).data(H==h);
        X(isnan(Y))=[];
        Y(isnan(Y))=[];
        Y(isnan(X))=[];
        X(isnan(X))=[];
        invLRcopula(m).h(h+1).corr=corr(X,Y);
        invLRcopula(m).h(h+1).miu1=mean(X); % Daily amplitude
        invLRcopula(m).h(h+1).sigma1=std(X);
        invLRcopula(m).h(h+1).miu2=mean(Y); % Lapse rate
        invLRcopula(m).h(h+1).sigma2=std(Y);
    end
end
save('Files\LR.mat','LR','invLRcopula','-v7.3');
%% Estimating mean areal air temperature (2-m) parameters
% For a reference field at 0 masl
dTbar=zeros(12,24); dTsigma=zeros(12,24);
progressbar('Calculate monthly mean areal air temperature (2-m) parameters') % Init single bar
for m=1:12
    [ Tb{m} , dTbar(m,:) , dTrho(m) , dTsigma(m,:) , Ti(m) , R2(m) ] = temperature_parameter( LR(m).T0 , CloudCover.data(month(CloudCover.time)==m) , LR(m).date , 1 , 8.4 , 46.8 );
    progressbar(m/12) % Update progress bar
end
save('Files\MeanArealTempParameters.mat','Tb','dTbar','dTrho','dTsigma','Ti','R2');
%% Estmating temperature AR(1)
mdl=arima(1,0,0);
progressbar('Estmating temperature AR(1) parameters') % Init single bar
for m=1:12
    EstMdl = estimate(mdl,LR(m).data);
    ARcoeffT(m,1)=EstMdl.AR{1};
    progressbar(m/12) % Update progress bar
end
save('Files\ARcoeffT.mat','ARcoeffT');
%% Shading effect
t=(datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,23,0,0))'; % Hourly dates for one year based on the year 2001
ShF=zeros(length(t),size(domainDTM,1),size(domainDTM,2),'uint8');
[HZ,Z]=HorizonAnglePolar(domainDTM,100,5,10000);
progressbar('Setting Sun variables')
for i=1:length(t)
    [Sun.altitude(i),Sun.declination(i),Sun.azimuth(i),Sun.sunrise(i),Sun.sunset(i),Sun.dayLength(i)]=SetSunVariables([year(t(i)),month(t(i)),day(t(i)),hour(t(i))],1,8.4,46.8,-1,2);
    Sun.date(i)=t(i);
    Sun.E0(i)=findE0(t(i));
    progressbar(i/length(t))
end
Sun.sunrise=Sun.sunrise+1; % GMT correction
Sun.sunset=Sun.sunset+1;
HZ=single(HZ);
save('Files\Shading.mat','HZ','Z','Sun','-v7.3');
%% Compute 30 years base temperature (at 0 masl) for differential lapse rate computation
Generate_base_Ta; % In this case Ta is a cell array of 30 years of hourly mean areal temperature at 0 masl
%% Estimating lapse rate
% Lapse rate is elevation depended calculated based on a comparison between
% ground stations at different elevation and the temperatue reference
% fields (at 0 masl) for every month and hour (to preserve the daily cycle)
progressbar('Calculating monthly lapse rate')
for m=1:12
    for h=0:23
        tmp=[];
        for I=1:30
            tmp=cat(1,tmp,Ta{I}(month(t)==m&hour(t)==h));
        end
        for II=1:5 % 5 stations for the Engelberger catchment
            tmp2=Ta_Stations{II}{3}(month(Ta_Stations{II}{3}(:,2))==m&hour(Ta_Stations{II}{3}(:,2))==h,1);
            estLR{II}(m).h(h+1).miu=(mean(tmp2)-mean(tmp));
            estLR{II}(m).h(h+1).sigma=std(tmp);
        end
    end
    progressbar(m/12)
end
progressbar('Estimating monthly lapse rate parameters change with elevation')
for m=1:12
    for h=0:23
        x=[]; y=[];
        for II=1:5
            x(II,1)=Ta_Stations{II}{2}/1000;
            y(II,1)=estLR{II}(m).h(h+1).miu;
        end
        tmp=fit(x,y,'poly1'); % Linear lapse rate correction with eleveation
        genLRpar(m).h(h+1).miuA=tmp.p1;
        genLRpar(m).h(h+1).miuB=tmp.p2;
        genLRpar(m).h(h+1).sigma=1; % Unit std
    end
    progressbar(m/12)
end
clear Ta estLR x y
save('Files\genLRpar.mat','genLRpar','-v7.3');
%% Calculate Ea and Esat (ambient and saturated vapor pressure)
for i=1:length(Rsw_Stations)
    Rsw_Stations(i).Esat=611.*exp(17.27.*Rsw_Stations(i).Ta./(237.3+Rsw_Stations(i).Ta));
    Rsw_Stations(i).Ea=Rsw_Stations(i).RH./100.*Rsw_Stations(i).Esat;
end
%% Computing deterministic and stochastic vapor coefficients
for i=1:length(Rsw_Stations)
    for m=1:12
        [ Vap(i).par(m).a , Vap(i).dDem(m) , Vap(i).rhodDe(m) , Vap(i).sigmadDe(m) ] = vap_pre_parameter( Rsw_Stations(i).Ea(month(Rsw_Stations(i).Time)==m) , Rsw_Stations(i).Esat(month(Rsw_Stations(i).Time)==m) , Rsw_Stations(i).Ta(month(Rsw_Stations(i).Time)==m) , Rsw_Stations(i).Rsw(month(Rsw_Stations(i).Time)==m) );
    end
end
% Estimating parameters changes with elevation
progressbar('Estimating monthly vapor parameters change with elevation')
for m=1:12
    for i=1:length(Rsw_Stations)
        dDem(i,1)=Vap(i).dDem(m);
        rhodDe(i,1)=Vap(i).rhodDe(m);
        sigmadDe(i,1)=Vap(i).sigmadDe(m);
        a0(i,1)=Vap(i).par(m).a(1);
        a1(i,1)=Vap(i).par(m).a(2);
        a2(i,1)=Vap(i).par(m).a(3);
        a3(i,1)=Vap(i).par(m).a(4);
        z(i,1)=Rsw_Stations(i).elevation; % [m]
    end
    Vapor.dDem{m}=mean(dDem);
    Vapor.rhodDe{m}=mean(rhodDe);
    Vapor.sigmadDe{m}=fit(z,sigmadDe,'poly1');
    Vapor.A0{m}=fit(z,a0,'poly1');
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1000 1];
    Vapor.A1{m}=fit(z,a1,'power1',opts);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.2 9e-5 1000 -0.0007];
    Vapor.A2{m}=fit(z,a2,'exp2',opts);
    Vapor.A3{m}=fit(z,a3,'poly1');
    progressbar(m/12)
end
clear Vap
save('Files\Vapor.mat','Vapor','-v7.3');
%% Estimating monthly AR(1) component for ground pressure
i=2; % Choose Engelberg station as reference for this example
mdlP=arima(1,0,0);
for m=1:12
    Pressure(m).ARMdl=estimate(mdlP,Ps_Stations(i).Ps(month(Ps_Stations(i).Time)==m,1));
end
clear mdlP
save('Files\Pressure.mat','Pressure','-v7.3');
%% Estimating sun variables
[Slo_top,Aspect]=Slope_Aspect_indexes(domainDTM,100,'grad');
[SvF,Ct]=Sky_View_Factor(domainDTM,atand(Slo_top),Aspect,HZ,Z);
save('Files\TerrainRadiation.mat','Slo_top','Aspect','SvF','Ct','-v7.3');
%% Fixed radiation parameters
Radiation.Ro=1360.8; % Solar constant, after Kopp 2011 [W m^-2]
Radiation.omega_A1=0.92;
Radiation.omega_A2=0.84;
Radiation.uo=0.35;
Radiation.un=0.0002;
Radiation.rho_g=0.15;
Radiation.AODcorrection=[1.8;1.3;1.3;1.45;1.75;2;1.5;1.35;1.4;1.3;1;1.25]; % Aerosol optical depth elevation correction, after Ingold (2001)
Radiation.AERONETz=866; % Elevation [m] of the aeronet station data
Radiation.lwp=[200;110;65;75;95;105;170;145;150;85;170;210]; % LWP that were calibrated using data derived from CM SAF as initial guess
%% Estimating alpha and beta (Angstorm turbidity parameters)
% Using AERONET aerosol optical depth data
lambda=[1.02;0.87;0.675;0.50;0.44;0.38;0.34]; % Wavelength (fit to AOD data)
progressbar('Fitting Angstorm turbidity parameters')
for m=1:12
    [ Alpha , Beta ] = fitAODpar( m , lambda , AOD );
    Radiation.Angstrom(m).alpha=nanmean(Alpha);
    Radiation.Angstrom(m).beta=nanmean(Beta);
    progressbar(m/12)
end
save('Files\Radiation.mat','Radiation','-v7.3');
%% Near ground wind
Hagl=2; % Height above ground level [m]
vk=0.4; % Von Karman constant
domainLU=rasterread('engelberger_land_use.txt'); % Land Use
domainLU=flipud(domainLU);
z0=nan(size(domainLU)); % Roughness map [m]
Z0_assign
cF=repmat(linspace(46.76,47,size(domainLU,1))',1,size(domainLU,1)); % Coriolis force
cF=sin(deg2rad(cF)).*1.458e-4;
Z0=unique(z0);
Z0(isnan(Z0))=[];
options=optimset('Display','off');
[ betaS , xiS ] = findSlope( domainDTM , 100 ); % Slope in [deg]
[ omegaC ] = findCurv( domainDTM , 100 , 1000 );
betaP=tand(betaS).*100; % Slope in [%]
betaP(betaP>100)=100;
gammaS=0.5; % Weghting parameter between Slope and Curves
gammaC=0.5;
save('Files\Wind.mat','vk','Hagl','domainLU','cF','z0','Z0','gammaC','gammaS','betaP','omegaC','betaS','xiS','options','-v7.3');
WindZ0Calibration % Calibrating the z0 (roughness coefficients) to match the 500 hPa field
%% Processing
I=1; % Number of realizations to simulate (the processing stage can be looped)
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
[ tsMonth , tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , 10/(24*60) );
for i=1:length(tsWetDry) % For this example the 5-min temporal resolution is interpolated from the 10-min resolution
    tsWetDry5Min(i*2-1)=tsWetDry(i);
    tsWetDry5Min(i*2)=tsWetDry(i);
    tsMonth5Min(i*2-1)=tsMonth(i);
    tsMonth5Min(i*2)=tsMonth(i);
end
tsWetDry=tsWetDry5Min;
tsMonth=tsMonth5Min;
clear tsMonth5Min tsWetDry5Min
[ tsMatrix ] = tsTable( tsWetDry , tsMonth );
%% Varfima - simulating WAR, CAR and IMF
simWAR=zeros(tsMatrix(end),1);
simIMF=zeros(tsMatrix(end),1);
simCAR=zeros(tsMatrix(end),1);
[ simDis ] = WetEventDis( simCAR , tsWetDry , 5 );
parfor i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if tsMatrix(i,1)==1
        [ tsWARFIMA{i} ] = MultiTriVARFIMA( MaternCov , Precipitation , Cloud , Copula , m , tsMatrix(i,4)-tsMatrix(i,3)+1 , 5 );
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
        simCAR(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.CARreal(:,E);
        simIMF(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.IMFreal(:,E);
    end
end
clear tsWARFIMA m
simCAR(simCAR<simWAR)=simWAR(simCAR<simWAR);
progressbar('Generating cloud cover') % Init single bar
for i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if tsMatrix(i,1)==0
        if i>1 && i<size(tsMatrix,1)
            try
                if simCAR(tsMatrix(i-1,4))>f_johnson_inv(.99999,Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type)
                    YS=f_johnson_cdf(f_johnson_inv(.99999,Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type),Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type);
                else
                    YS=f_johnson_cdf(simCAR(tsMatrix(i-1,4)),Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type);
                end
            catch
                a=simCAR(tsMatrix(i-1,4));
                YS=Cloud.Dry(m).fit.cdf(a);
            end
            YS=(YS-Cloud.ar(m).disMiu(0))/Cloud.ar(m).disStd;
            [ S ] = simulate(Cloud.ar(m).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1,'NumPaths',50,'Y0',repmat(YS,Cloud.ar(m).mdl.P,1));
            normClt=bsxfun(@plus,S.*Cloud.ar(m).disStd,Cloud.ar(m).disMiu(simDis(tsMatrix(i,3):tsMatrix(i,4))));
            normClt=normcdf(normClt,repmat(mean(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1),repmat(std(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1));
            try
                invS=Cloud.Dry(m).fit.icdf(normClt);
            catch
                jType=Cloud.Dry(m).fit.type;
                invS=f_johnson_inv(normClt(:),Cloud.Dry(m).fit.coef,jType);
                invS=reshape(invS,tsMatrix(i,4)-tsMatrix(i,3)+1,50);
            end
            invS(invS>1)=1;
            invS(invS<0)=0;
            SV=simCAR(tsMatrix(i-1,4)); EV=simCAR(tsMatrix(i+1,3));
            E=abs(SV-invS(1,:))+abs(EV-invS(end,:));
            E=find(min(E)==E);
            simCAR(tsMatrix(i,3):tsMatrix(i,4),1)=invS(:,E(1));
        end
        if i==1
            [ S ] = simulate(Cloud.ar(m).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1,'NumPaths',50);
            normClt=bsxfun(@plus,S.*Cloud.ar(m).disStd,Cloud.ar(m).disMiu(simDis(tsMatrix(i,3):tsMatrix(i,4))));
            normClt=normcdf(normClt,repmat(mean(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1),repmat(std(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1));
            try
                invS=Cloud.Dry(m).fit.icdf(normClt);
            catch
                jType=Cloud.Dry(m).fit.type;
                invS=f_johnson_inv(normClt(:),Cloud.Dry(m).fit.coef,jType);
                invS=reshape(invS,tsMatrix(i,4)-tsMatrix(i,3)+1,50);
            end
            invS(invS>1)=1;
            invS(invS<0)=0;
            EV=simCAR(tsMatrix(i+1,3));
            E=abs(EV-invS(end,:));
            E=find(min(E)==E);
            simCAR(tsMatrix(i,3):tsMatrix(i,4),1)=invS(:,E(1));
        end
        if i==size(tsMatrix,1)
            try
                if simCAR(tsMatrix(i-1,4))>f_johnson_inv(.99999,Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type)
                    YS=f_johnson_cdf(f_johnson_inv(.99999,Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type),Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type);
                else
                    YS=f_johnson_cdf(simCAR(tsMatrix(i-1,4)),Cloud.Dry(m).fit.coef,Cloud.Dry(m).fit.type);
                end
            catch
                a=simCAR(tsMatrix(i-1,4));
                YS=Cloud.Dry(m).fit.cdf(a);
            end
            YS=(YS-Cloud.ar(m).disMiu(0))/Cloud.ar(m).disStd;
            [ S ] = simulate(Cloud.ar(m).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1,'NumPaths',50,'Y0',repmat(YS,Cloud.ar(m).mdl.P,1));
            normClt=bsxfun(@plus,S.*Cloud.ar(m).disStd,Cloud.ar(m).disMiu(simDis(tsMatrix(i,3):tsMatrix(i,4))));
            normClt=normcdf(normClt,repmat(mean(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1),repmat(std(normClt),tsMatrix(i,4)-tsMatrix(i,3)+1,1));
            try
                invS=Cloud.Dry(m).fit.icdf(normClt);
            catch
                jType=Cloud.Dry(m).fit.type;
                invS=f_johnson_inv(normClt(:),Cloud.Dry(m).fit.coef,jType);
                invS=reshape(invS,tsMatrix(i,4)-tsMatrix(i,3)+1,50);
            end
            invS(invS>1)=1;
            invS(invS<0)=0;
            SV=simCAR(tsMatrix(i-1,4));
            E=abs(SV-invS(1,:));
            E=find(min(E)==E);
            simCAR(tsMatrix(i,3):tsMatrix(i,4),1)=invS(:,E(1));
        end
    end
    progressbar(i/size(tsMatrix,1)) % Update progress bar
end
%% Simulating advection
simU=zeros(tsMatrix(end),1);
simV=zeros(tsMatrix(end),1);
parfor i=1:size(tsMatrix,1)
    Sd{i}=vgxsim(Advection.Dry(tsMatrix(i,2)).mdl,20000);
    Sw{i}=vgxsim(Advection.Wet(tsMatrix(i,2)).mdl,20000);
end
progressbar('Generating advection components') % Init single bar
for i=1:size(tsMatrix,1)
    m=tsMatrix(i,2);
    if i==1
        switch tsMatrix(i,1)
            case 0
                Ss=vgxsim(Advection.Dry(tsMatrix(i,2)).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1);
                YS=Ss(end,:);
                Ss=normcdf(Ss);
                Ss=round(Ss.*1000)./1000;
                Ss(Ss==1)=0.999;
                Ss(Ss==0)=0.001;
                Ss(:,1)=norminv(Ss(:,1),Advection.Dry(m).unorm(1),Advection.Dry(m).unorm(2));
                Ss(:,2)=norminv(Ss(:,2),Advection.Dry(m).vnorm(1),Advection.Dry(m).vnorm(2));
            case 1
                Ss=vgxsim(Advection.Wet(tsMatrix(i,2)).mdl,tsMatrix(i,4)-tsMatrix(i,3)+1);
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
simRain=zeros(105120,13,13);
progressbar('Generating rain fields') % Init single bar
for m=1:12
    mmin(m,1)=min(tsMatrix(tsMatrix(:,2)==m,3));
    mmax(m,1)=max(tsMatrix(tsMatrix(:,2)==m,4));
    [ QField(mmin(m,1):mmax(m,1),:,:) ] = quantileFieldGen( [26 26] , ARMA , 5 , Precipitation.data(m).SpatialAlpha , [mmin(m,1) mmax(m,1)] , simU , simV , 2 , 5 );
    %% Non-homogeneous probability of precipitation occurrence
    for i=mmin(m,1):mmax(m,1)
        [ simRain(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simWAR(i) , obsOCC{m} );
    end
    %% Assigning rainfall intensity for the Gaussian field
    for i=1:size(tsMatrix,1)
        M=tsMatrix(i,2);
        if tsMatrix(i,1)==1 && m==M
            [ simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) ] = invLN2( simRain(tsMatrix(i,3):tsMatrix(i,4),:,:) , simWAR(tsMatrix(i,3):tsMatrix(i,4)) , simIMF(tsMatrix(i,3):tsMatrix(i,4)) , CV(m));
        end
    end
    progressbar(m/12) % Update progress bar
end
%% Non-homogeneous spatial rainfall accumulation
% Simulation is divided to daily time steps to match the correction records
t=(datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,23,55,0))'; % One year dates based on the year 2001
tmp=squeeze(sum(reshape(simRain,[288 365 13 13]))).*5/60;
for i=1:13*13
    for d=1:365
        if tmp(d,i)<thresholdI(i)
            k=tmp(d,i)/thresholdI(i);
            simRain(1+288*(d-1):288*d,i)=simRain(1+288*(d-1):288*d,i).*k; % Below the 1 mm threshold is consider as "drizzle"
            tmp(d,i)=0;
        end
    end
end
for d=1:365
    m=month(t(d));
    tmp2=squeeze(tmp(d,:,:));
    gpCDF=gpcdf(tmp2,gpFitS{1}{m},gpFitS{2}{m},thresholdI);
    gpCDF(gpCDF>0.999)=0.999;
    gpINV=gpinv(gpCDF,gpFitO{1}{m},gpFitO{2}{m},1);
    k=gpINV./tmp2; % Daily rainfall correction factor
    k(isinf(k))=1;
    simRain(1+288*(d-1):288*d,:,:)=simRain(1+288*(d-1):288*d,:,:).*repmat(reshape(k,[1 13 13]),288,1,1);
end
simRain(isnan(simRain))=0;
simRain=single(simRain);
%% Non-homogeneous probability of cloud cccurrence
simCloud=zeros(105120,13,13);
progressbar('Generating cloud fields') % Init single bar
for m=1:12
    for i=mmin(m,1):mmax(m,1)
        [ simCloud(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simCAR(i) , obsOCC{m} ); % Using same occurence filter as was used for the WAR, but for CAR variable
    end
    progressbar(m/12) % Update progress bar
end
simCloud(simCloud>0)=1;
%% Saving output
save('Files\Output_Pr_Clt_500UV.mat','simCloud','simRain','simIMF','simCAR','simWAR','simU','simV','-v7.3'); % Matlab
clear simRain
%% Converting simulated cloud cover from 5-min resolution to 1-h resolution
Nsim=mean(reshape(simCAR,12,[]))';
%% Generating near surface air temperature (2-m) fields
HZr=HZ*pi/180;
t=(datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,23,0,0))'; % Hourly dates for one year based on the year 2001
Tmax=zeros(length(Nsim),1); Tmin=zeros(length(Nsim),1);
Ta=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2));
mdl=arima('Constant',0,'Variance',1,'AR',{mean(ARcoeffT)}); % Hourly lapse rate time series
v = simulate(mdl,length(Ta));
v = normcdf(v,mean(v),std(v));
progressbar('Generating near surface air temperature fields') % Init single bar
for i=1:length(Nsim)
    m=month(t(i));
    h=hour(t(i));
    GenLR=norminv(v(i),genLRpar(m).h(h+1).miuA.*(domainDTM./1000)+genLRpar(m).h(h+1).miuB,genLRpar(m).h(h+1).sigma); % Lapse rate correction
    ShF=ShadowEffect(domainDTM,Sun.altitude(i),Sun.azimuth(i),HZr,Z);
    if i==1
        % Can also compute using ComputeAirTemperatureFast without the
        % shading effect (see example file:
        % Processing_without_shading_effect.m
        [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},zeros(size(domainDTM)),dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),0,ones(size(domainDTM))*Ti(month(t(i))),ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
        Ta(i,:,:)=Tsim+GenLR;
        TI=ones(size(domainDTM))*Ti(month(t(i)));
    else
        if hour(t(i))==0
            TI=T_tilde;
            I2=zeros(size(domainDTM));
            I3=zeros(size(domainDTM));
            I4=zeros(size(domainDTM));
        end
        [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},qt,dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),dT,TI,I2,I3,I4,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
        Ta(i,:,:)=Tsim+GenLR;
    end
    Tmax(i,1)=max(max(max(Tsim)));
    Tmin(i,1)=min(min(min(Tsim)));
    progressbar(i/length(Nsim)) % Update progress bar
end
%% Calculating the daily temperature amplitude
dTmax=datevec(t);
dTmax=datenum(dTmax(:,1),dTmax(:,2),dTmax(:,3),0,0,0);
dTmax2=unique(dTmax);
for i=1:length(dTmax2)
    dailyAmp(i,1)=max(Tmax(dTmax==dTmax2(i)));
    dailyAmp(i,2)=min(Tmax(dTmax==dTmax2(i)));
    dailyAmp(i,3)=dailyAmp(i,1)-dailyAmp(i,2);
end
%% Thermal inversion - correcting for positive lapse rate
x=zeros(length(Nsim),1); invLR=zeros(length(Nsim),1); u=zeros(length(Nsim),1);
for i=1:length(dTmax2)
    x(dTmax2(i)==dTmax,1)=dailyAmp(i,3);
end
v = simulate(mdl,length(x)); % Simulating hourly lapse rate time series copula with daily temperature amplitude
v = norminv(normcdf(v,mean(v),std(v)),0,1);
for i=1:length(x)
    m=month(t(i));
    h=hour(t(i));
    u(i,1) = norminv(normcdf(x(i),invLRcopula(m).h(h+1).miu1,invLRcopula(m).h(h+1).sigma1),0,1);
    invLR(i,1) = invLRcopula(m).h(h+1).sigma2.*(invLRcopula(m).h(h+1).corr.*u(i,1)+sqrt(1-invLRcopula(m).h(h+1).corr.^2).*v(i,1))+invLRcopula(m).h(h+1).miu2;
end
for i=1:length(Nsim)
    if invLR(i)>0 % In case of thermal inversion
        tmp=squeeze(Ta(i,:,:));
        tmp(domainDTM<=Z1)=tmp(domainDTM<=Z1)-invLR(i).*((Z1(domainDTM<=Z1)-domainDTM(domainDTM<=Z1))./1000);
        Ta(i,:,:)=tmp;
    end
end
clear x mdl v u m h dailyAmp dTmax dTmax2 Tmax Tmin
Ta=single(Ta);
%% Saving output
save('Files\Output_Ta.mat','Ta','-v7.3'); % Matlab
%% Computing pressure at ground level
Ps=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Atmoshperic pressure field [Pa] at ground level
for m=1:12 % Generating pressure at reference point
    p_sim(month(t)==m)=simulate(Pressure(m).ARMdl,sum(month(t)==m));
end
progressbar('Generating pressure at ground level fields') % Init single bar
for i=1:length(Nsim)
    [ P ] = computePressureField( p_sim(i) , Ps_Stations(2).elevation , squeeze(Ta(i,:,:)) , domainDTM ); % Reference point for this example, Engelberg station [Ps_Stations(2)]
    Ps(i,:,:)=P;
    progressbar(i/length(Nsim)) % Update progress bar
end
Ps=single(Ps);
%% Saving output
save('Files\Output_Ps.mat','Ps','-v7.3'); % Matlab format
clear Ps
%% Generating vapor pressure and radiation (dependency among them)
dDe=zeros(length(Nsim),1);
Rsw_tm1=zeros(size(domainDTM));
Rsw_tm2=zeros(size(domainDTM));
HZr=HZ*pi/180;
M=0;
ea=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Vapor pressure [Pa] field at ground level
ea=single(ea);
ssr=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Surface shortwave incoming radaiation [W m^-2] field at ground level
ssr=single(ssr);
Rsw=zeros(size(domainDTM));
progressbar('Generating vapor pressure and incoming global radiation fields') % Init single bar
for i=1:length(Nsim)
    Ta_i=squeeze(Ta(i,:,:));
    m=month(t(i));
    %% Vapor pressure Generator
    esat_s=611.*exp(17.27.*Ta_i./(237.3+Ta_i));
    if m~=M
        %% Generating monthly vapor pressure data
        a0=reshape(Vapor.A0{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
        a1=reshape(Vapor.A1{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
        a2=reshape(Vapor.A2{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
        a3=reshape(Vapor.A3{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
        dDem=Vapor.dDem{m};
        rhodDe=Vapor.rhodDe{m};
        sigmadDe=reshape(Vapor.sigmadDe{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
        M=m;
    end
    if i<3
        [ ea(i,:,:) , dDe(i) ] = ComputeVapPressure2( esat_s , Ta_i , 0 , 0 , 0 , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
    else
        [ ea(i,:,:) , dDe(i) ] = ComputeVapPressure2( esat_s , Ta_i , Rsw_tm1 , Rsw_tm2 , dDe(i-1) ,a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
    end
    U=squeeze(ea(i,:,:))./esat_s;
    G_g=17.27.*Ta_i./(237.7+Ta_i)+log(U);
    Tdew=237.7.*G_g./(17.27-G_g);
    Tdew(isinf(Tdew))=NaN;
    %% Radiation Generator
    if Sun.altitude(i)>0 || i==1
        [ cos_ST ] = SolarIlluminationAngle( Slo_top , Sun.altitude(i) , Sun.azimuth(i) , Aspect );
        [ S ] = Shadow_Effect2( domainDTM , pi/2-Sun.altitude(i) , Sun.azimuth(i) , HZr , Z );
        [SB,SD,SAB1,SAB2,SAD1,SAD2,PARB,PARD]=ComputeRadiationForcings( Sun.altitude(i) , Sun.E0(i) , Radiation.Ro , Radiation.lwp(m) , Nsim(i) , domainDTM , Tdew , Radiation.Angstrom(m).beta , Radiation.Angstrom(m).alpha , Radiation.omega_A1 ,Radiation.omega_A2 , Radiation.uo , Radiation.un , Radiation.rho_g , cos_ST , S , SvF , Ct , Radiation.AODcorrection(m) , Radiation.AERONETz );
        SB(SB<0)=0;
        SD(SD<0)=0;
        Rsw=SB+SD;
    else
        Rsw=zeros(size(domainDTM));
    end
    Rsw_tm2=Rsw_tm1;
    Rsw_tm1=Rsw;
    ssr(i,:,:)=Rsw;
    progressbar(i/length(Nsim)) % Update progress bar
end
%% Saving output
save('Files\Output_ssr.mat','ssr','-v7.3'); % Matlab format
clear Rsw
for i=1:length(Nsim)
    Rsw(i,1)=nanmean(nanmean(ssr(i,:,:)));
end
clear ssr
%% Relative humidity
hur=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Relative humidity [%] field at ground level
hur=single(hur);
progressbar('Generating relative humidity fields') % Init single bar
for i=1:length(Nsim)
    esat_s=611.*exp(17.27.*squeeze(Ta(i,:,:))./(237.3+squeeze(Ta(i,:,:))));
    hur(i,:,:)=squeeze(ea(i,:,:))./esat_s;
    progressbar(i/length(Nsim)) % Update progress bar
end
%% Saving output
save('Files\Output_hur.mat','hur','-v7.3'); % Matlab format
save('Files\Output_Ea.mat','ea','-v7.3'); % Matlab format
clear ea
%% Dew point temperature
Tdew=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Dew point temperature [°C] field at ground level
Tdew=single(Tdew);
progressbar('Generating dew point temperature fields') % Init single bar
for i=1:length(Nsim)
    hur_lim=squeeze(hur(i,:,:));
    hur_lim(hur_lim<0.01)=0.01; % Minimum relative humidity for the log
    G_g=17.27.*squeeze(Ta(i,:,:))./(237.7+squeeze(Ta(i,:,:)))+log(hur_lim);
    Tdew(i,:,:)=237.7.*G_g./(17.27-G_g);
    progressbar(i/length(Nsim)) % Update progress bar
end
Tdew(isinf(Tdew))=NaN;
save('Files\Output_Tdew.mat','Tdew','-v7.3'); % Matlab format
clear hur Tdew Ta
%% Near surface wind speed
uRot=single(zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)));
vRot=single(zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)));
mcF=mean2(cF);
G=sqrt(simU.^2+simV.^2); % Geostrophic wind speed [m s^-1]
G=mean(reshape(G,[],8760))'; % Geostrophic wind speed for 1-h intervals [m s^-1]
[ classP ] = findPasquillClass( Nsim , Rsw ); % Find potential Pasquill classes for each hour
for i=1:length(Z0)
    z_cond{i}=(z0==Z0(i));
end
for m=1:12
    Z0_corr{m}=w2z{m}(z2w{m}(Z0));
end
%%  Monin-Obukov lenght L
for m=1:12
    z0_corr{m}=reshape(w2z{m}(z2w{m}(z0)),size(domainDTM,1),size(domainDTM,2));
    MOL1{m}=-11.4.*z0_corr{m}.^0.1;
    MOL2{m}=-26.*z0_corr{m}.^0.17;
    MOL3{m}=-123.*z0_corr{m}.^0.3;
    MOL5{m}=123.*z0_corr{m}.^0.3;
    MOL6{m}=26.*z0_corr{m}.^0.17;
end
%% Iterations
Ws=zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)); % Wind speed [m s^-1] field at near ground level (2-m elevation)
tF=0.01:0.01:5;
for i=1:length(Nsim)
    m=month(t(i));
    cP=classP{i}';
    FVi=[];
    if length(cP)>1
        for cPi=1:length(cP)
            switch cP(cPi)
                case 4 % Neutral atmosphere
                    FV=nan(size(domainDTM));
                    bU=4;
                    b_U=-4.5;
                    aU=10;
                    a_U=-5.5;
                    AU=4.5;
                    BU=0.5;
                    for j=1:length(Z0)
                        [ fun ] = nF4(tF,vk,G(i),AU,BU,mcF,Z0_corr{m}(j));
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    MOL=0; % Monin-Obukov lenght L
                    H=0.3.*(FV./abs(cF)); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                case 5 % Stable atmosphere
                    FV=nan(size(domainDTM));
                    MOL=123.*Z0_corr{m}.^0.3;
                    for j=1:length(Z0)
                        [ fun ] = nF5( tF,vk,G(i),mcF,Z0_corr{m}(j),MOL(j) );
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    coeffU=(vk.*FV)./(abs(cF).*MOL5{m});
                    bU=4+10.2.*sqrt(coeffU);
                    b_U=-4.5-7.65.*sqrt(coeffU);
                    aU=10;
                    a_U=-5.5+1.765.*sqrt(coeffU);
                    H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                case 6 % Very stable atmosphere
                    FV=nan(size(domainDTM));
                    MOL=26.*Z0_corr{m}.^0.17;
                    for j=1:length(Z0)
                        fun = nF6(tF,vk,G(i),mcF,Z0_corr{m}(j),MOL(j));
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    coeffU=(vk.*FV)./(abs(cF).*MOL6{m});
                    bU=4+10.2.*sqrt(coeffU);
                    b_U=-4.5-7.65.*sqrt(coeffU);
                    aU=10;
                    a_U=-5.5+1.765.*sqrt(coeffU);
                    H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                case 1 % Extreme unstable atmosphere
                    FV=nan(size(domainDTM));
                    MOL=-11.4.*Z0_corr{m}.^0.1;
                    for j=1:length(Z0)
                        fun = nF1(tF,vk,G(i),mcF,Z0_corr{m}(j),MOL(j));
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    coeffU=(vk.*FV)./(abs(cF).*MOL1{m});
                    bU=-34+38./(1+0.027.*sqrt(-coeffU));
                    b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                    aU=10./(1+1.581.*sqrt(-coeffU));
                    a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                    H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                case 2 % Very unstable atmosphere
                    FV=nan(size(domainDTM));
                    MOL=-26.*Z0_corr{m}.^0.17;
                    for j=1:length(Z0)
                        fun = nF2(tF,vk,G(i),mcF,Z0_corr{m}(j),MOL(j));
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    coeffU=(vk.*FV)./(abs(cF).*MOL2{m});
                    bU=-34+38./(1+0.027.*sqrt(-coeffU));
                    b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                    aU=10./(1+1.581.*sqrt(-coeffU));
                    a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                    H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                case 3 % Unstable atmosphere
                    FV=nan(size(domainDTM));
                    MOL=-123.*Z0_corr{m}.^0.3;
                    for j=1:length(Z0)
                        fun = nF3(tF,vk,G(i),mcF,Z0_corr{m}(j),MOL(j));
                        tmp=tF(min(abs(fun))==abs(fun));
                        FV(z_cond{j})=tmp(1);
                    end
                    FVi{cP(cPi)}=FV;
                    coeffU=(vk.*FV)./(abs(cF).*MOL3{m});
                    bU=-34+38./(1+0.027.*sqrt(-coeffU));
                    b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                    aU=10./(1+1.581.*sqrt(-coeffU));
                    a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                    H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
                    cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
            end
        end
    end
    if length(cP)>1
        if i==1 || i==8760
            if Nsim(i)<0.5
                cP(cP(:,2)<=2,3)=6;
                cP(cP(:,2)>2&cP(:,2)<=3,3)=6;
                cP(cP(:,2)>3&cP(:,2)<=5,3)=5;
                cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                cP(cP(:,2)>=6,3)=4;
            else
                cP(cP(:,2)<=2,3)=6;
                cP(cP(:,2)>2&cP(:,2)<=3,3)=5;
                cP(cP(:,2)>3&cP(:,2)<=5,3)=4;
                cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                cP(cP(:,2)>=6,3)=4;
            end
        else
            if Rsw(i)>0 % Day
                if Rsw(i)<50
                    cP(cP(:,2)<=2,3)=3;
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=3;
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=3;
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                    cP(cP(:,2)>=6,3)=4;
                elseif Rsw(i)>=50 && Rsw(i)<300
                    cP(cP(:,2)<=2,3)=2;
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=3;
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=3;
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                    cP(cP(:,2)>=6,3)=4;
                elseif Rsw(i)>=300 && Rsw(i)<600
                    cP(cP(:,2)<=2,3)=round(unifrnd(1,2,1,1));
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=2;
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=round(unifrnd(2,3,1,1));
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=round(unifrnd(3,4,1,1));
                    cP(cP(:,2)>=6,3)=4;
                else
                    cP(cP(:,2)<=2,3)=1;
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=round(unifrnd(1,2,1,1));
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=2;
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=3;
                    cP(cP(:,2)>=6,3)=3;
                end
            elseif Rsw(i)==0 && Rsw(i-1)==0 && Rsw(i+1)==0 % Night
                if Nsim(i)<0.5
                    cP(cP(:,2)<=2,3)=6;
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=6;
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=5;
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                    cP(cP(:,2)>=6,3)=4;
                else
                    cP(cP(:,2)<=2,3)=6;
                    cP(cP(:,2)>2&cP(:,2)<=3,3)=5;
                    cP(cP(:,2)>3&cP(:,2)<=5,3)=4;
                    cP(cP(:,2)>5&cP(:,2)<=6,3)=4;
                    cP(cP(:,2)>=6,3)=4;
                end
            end
        end
        cP(cP==0)=nan;
        classP{i}=mode(cP(:,3));
    end
    %% Calculate near ground wind (2-m) uv components
    switch classP{i}
        case 4 % Neutral atmosphere
            bU=4;
            b_U=-4.5;
            aU=10;
            a_U=-5.5;
            AU=4.5;
            BU=0.5;
            try
                FV=FVi{classP{i}};
            catch
                FV=nan(size(domainDTM));
                for j=1:length(Z0)
                    fun = nF4(tF,vk,G(i),AU,BU,mcF,Z0_corr{m}(j));
                    tmp=tF(min(abs(fun))==abs(fun));
                    FV(z_cond{j})=tmp(1);
                end
            end
            MOL=0; % Monin-Obukov lenght L
            H=0.3.*(FV./abs(cF)); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
        case 5 % Stable atmosphere
            FV=FVi{classP{i}};
            coeffU=(vk.*FV)./(abs(cF).*MOL5{m});
            bU=4+10.2.*sqrt(coeffU);
            b_U=-4.5-7.65.*sqrt(coeffU);
            aU=10;
            a_U=-5.5+1.765.*sqrt(coeffU);
            H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
            AU=4.5+1.765.*sqrt(coeffU);
            BU=0.5-2.55.*sqrt(coeffU);
        case 6 % Very stable atmosphere
            FV=FVi{classP{i}};
            coeffU=(vk.*FV)./(abs(cF).*MOL6{m});
            bU=4+10.2.*sqrt(coeffU);
            b_U=-4.5-7.65.*sqrt(coeffU);
            aU=10;
            a_U=-5.5+1.765.*sqrt(coeffU);
            H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
            AU=4.5+1.765.*sqrt(coeffU);
            BU=0.5-2.55.*sqrt(coeffU);
        case 1 % Extreme unstable atmosphere
            FV=FVi{classP{i}};
            coeffU=(vk.*FV)./(abs(cF).*MOL1{m});
            bU=-34+38./(1+0.027.*sqrt(-coeffU));
            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
            aU=10./(1+1.581.*sqrt(-coeffU));
            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
            AU=4.5./(1+1.581.*sqrt(-coeffU));
            BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
        case 2 % Very unstable atmosphere
            FV=FVi{classP{i}};
            coeffU=(vk.*FV)./(abs(cF).*MOL2{m});
            bU=-34+38./(1+0.027.*sqrt(-coeffU));
            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
            aU=10./(1+1.581.*sqrt(-coeffU));
            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
            AU=4.5./(1+1.581.*sqrt(-coeffU));
            BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
        case 3 % Unstable atmosphere
            FV=FVi{classP{i}};
            coeffU=(vk.*FV)./(abs(cF).*MOL3{m});
            bU=-34+38./(1+0.027.*sqrt(-coeffU));
            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
            aU=10./(1+1.581.*sqrt(-coeffU));
            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
            u=(FV./vk).*(log(Hagl./z0_corr{m})+bU.*((Hagl-z0_corr{m})./H)+b_U.*((Hagl-z0_corr{m})./H).^2);
            v=(-FV./vk).*(aU.*((Hagl-z0_corr{m})./H)+a_U.*((Hagl-z0_corr{m})./H).^2);
            AU=4.5./(1+1.581.*sqrt(-coeffU));
            BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
    end
    %% Transfrom to Cartesian coordinate system
    alphaWind=-atand(AU./(log(FV./(abs(cF).*z0_corr{m}))-BU)); % Angle between geostrophic wind vector and the axis of the friction velocity [deg]
    if simV(i)>0 && simU(i)>0 % Angle between the geostrophic wind and X-axis (eastward)
        betaWind=atand(simV(i)./simU(i));
    elseif simV(i)<0 && simU(i)>0
        betaWind=atand(simV(i)./simU(i));
    elseif simV(i)<0 && simU(i)<0
        betaWind=-(atand(simV(i)./simU(i))+90);
    elseif simV(i)>0 && simU(i)<0
        betaWind=atand(simV(i)./simU(i))-90;
    elseif simV(i)==0 && simU(i)>0
        betaWind=atand(simV(i)./simU(i));
    elseif simV(i)==0 && simU(i)<0
        betaWind=atand(simV(i)./simU(i))-180;
    elseif simV(i)>0 && simU(i)==0
        betaWind=atand(simV(i)./simU(i));
    elseif simV(i)<0 && simU(i)==0
        betaWind=atand(simV(i)./simU(i));
    end
    rotAngle=deg2rad(alphaWind+betaWind); % Positive - counterclock
    uRot(:,:,i)=u.*cos(rotAngle)-v.*sin(rotAngle);
    vRot(:,:,i)=u.*sin(rotAngle)+v.*cos(rotAngle);
    progressbar(i/length(Nsim)) % Update progress bar
end
%% Correction of slope and curve
parfor i=1:length(Nsim)
    [ W , tet ] = z2p( squeeze(uRot(:,:,i)) , squeeze(vRot(:,:,i)) );
    [ omegaS ] = findOmegaS( betaS , tet , xiS );
    Ww=1+gammaS.*omegaS+gammaC.*omegaC;
    Ws(:,:,i)=Ww.*W;
end
clear uRot vRot
%% Saving output
save('Files\Output_Ws.mat','Ws','-v7.3'); % Matlab format
clear Ws
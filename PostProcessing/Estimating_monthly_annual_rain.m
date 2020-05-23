%Estimating_monthly_annual_rain This script run the rainfall simulator for
%30 years and compute the monthly and annual rainfall. The simulated fields
%are compared to the observed and plotted
addpath(genpath(cd))
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
load('thresholdI.mat')
load('gpFitSO.mat') % Generalized Pareto parameters for the observed and simulated daily data
%% Processing
%% Initilazing
rng('shuffle');
t=(datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,23,55,0))'; % Dummy - One year dates based on the year 2001
for I=1:30
    for m=1:12
        mR{I}{m}=[]; % Rain monthly
        aR{I}=[]; % Rain annualy
    end
end
I=30; % Number of years to simulate
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
progressbar('Generating rain field ensemble') % Init single bar
for I=1:30
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
    %% Advection - constant 5 m s^-1 toward the East
    simU=ones(tsMatrix(end),1).*5;
    simV=zeros(tsMatrix(end),1);
    simU=single(simU);
    simV=single(simV);
    %% Varfima - simulating WAR, CAR and IMF
    simWAR=zeros(tsMatrix(end),1);
    simCAR=zeros(tsMatrix(end),1);
    simIMF=zeros(tsMatrix(end),1);
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
            simIMF(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA{i}.IMFreal(:,E);
        end
    end
    clear tsWARFIMA m
    %% Generating rain fields
    rField=zeros(105120,13,13); mmin=zeros(12,1); mmax=zeros(12,1); QField=zeros(105120,13,13);
    for m=1:12
        mmin(m,1)=min(tsMatrix(tsMatrix(:,2)==m,3));
        mmax(m,1)=max(tsMatrix(tsMatrix(:,2)==m,4));
        [ QField(mmin(m,1):mmax(m,1),:,:) ] = quantileFieldGen( [26 26] , ARMA , 5 , Precipitation.data(m).SpatialAlpha , [mmin(m,1) mmax(m,1)] , simU , simV , 2 , 5 );
        %% Non-homogeneous probability of precipitation cccurrence
        for i=mmin(m,1):mmax(m,1)
            [ rField(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simWAR(i) , obsOCC{m} );
        end
        %% Assigning rainfall intensity for the Gaussian field
        for i=1:size(tsMatrix,1)
            M=tsMatrix(i,2);
            if tsMatrix(i,1)==1 && m==M
                [ rField(tsMatrix(i,3):tsMatrix(i,4),:,:) ] = invLN2( rField(tsMatrix(i,3):tsMatrix(i,4),:,:) , simWAR(tsMatrix(i,3):tsMatrix(i,4)) , simIMF(tsMatrix(i,3):tsMatrix(i,4)) , CV(m));
            end
        end
    end
    %% Saving the daily data per pixel
    tmp=squeeze(sum(reshape(rField,[288 365 13 13]))).*5/60;
    for i=1:13*13
        for d=1:365
            if tmp(d,i)<thresholdI(i)
                k=tmp(d,i)/thresholdI(i);
                rField(1+288*(d-1):288*d,i)=rField(1+288*(d-1):288*d,i).*k; % Below the 1 mm threshold is consider as "drizzle"
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
        rField(1+288*(d-1):288*d,:,:)=rField(1+288*(d-1):288*d,:,:).*repmat(reshape(k,[1 13 13]),288,1,1);
    end
    rField(isnan(rField))=0;
    tmp=squeeze(sum(reshape(rField,[288 365 13 13]))); % For daily setup: 288 5-min time steps, 365 days and 13 x 13 domain size
    aR{I}=squeeze(sum(tmp).*5/60);
    for m=1:12
        mR{I}{m}=squeeze(sum(tmp(month(t)==m,:,:)).*5/60);
    end
    progressbar(I/30) % Update progress bar
end
%% Save results
save('Files\Est_mon_ann_rain.mat','mR','aR','-v7.3');
%% Load observations
load('ObsAnnualRain.mat')
ri=zeros(13,13); % Mean annual rain [mm]
ri_std=zeros(13,13); % Standard deviation annual rain [mm]
ri_min=ones(13,13).*1e6; % Min annual rain [mm]
ri_max=zeros(13,13); % Max annual rain [mm]
for i=1:32
    ri=ri+ObsAnnualRain{i};
    ri_min(ri_min>ObsAnnualRain{i})=ObsAnnualRain{i}(ri_min>ObsAnnualRain{i});
    ri_max(ri_max<ObsAnnualRain{i})=ObsAnnualRain{i}(ri_max<ObsAnnualRain{i});
end
ri=ri./32;
for i=1:13*13
    clear tmp
    for j=1:32
        tmp(j)=ObsAnnualRain{j}(i);
    end
    ri_std(i)=std(tmp);
end
%% Postprocess simulated data
riS=zeros(13,13); % Mean annual rain [mm]
riS_std=zeros(13,13); % Standard deviation annual rain [mm]
riS_min=ones(13,13).*1e6; % Min annual rain [mm]
riS_max=zeros(13,13); % Max annual rain [mm]
for i=1:30
    riS=riS+aR{i};
    riS_min(riS_min>aR{i})=aR{i}(riS_min>aR{i});
    riS_max(riS_max<aR{i})=aR{i}(riS_max<aR{i});
end
riS=riS./30;
for i=1:13*13
    clear tmp
    for j=1:30
        tmp(j)=aR{j}(i);
    end
    riS_std(i)=std(tmp);
end
% monthly
for m=1:12
    riM{m}=zeros(13,13); % Mean monthly rain [mm]
    riM_std{m}=zeros(13,13); % Standard deviation monthly rain [mm]
end
for m=1:12
    for i=1:30
        riM{m}=riM{m}+mR{i}{m};
    end
end
for m=1:12
    riM{m}=riM{m}./30;
end
for m=1:12
    for i=1:13*13
        clear tmp
        for j=1:30
            tmp(j)=mR{j}{m}(i);
        end
        riM_std{m}(i)=std(tmp);
    end
end
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Plot annual maps
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(2,2,1)
imagesc(ri)
colorbar
axis('square')
title('\mu_{obs} [mm]')
caxis([1400 2200])
subplot(2,2,2)
imagesc(ri_std)
colorbar
axis('square')
title('\sigma_{obs} [mm]')
caxis([150 350])
subplot(2,2,3)
imagesc(riS)
colorbar
axis('square')
title('\mu_{sim} [mm]')
caxis([1400 2200])
subplot(2,2,4)
imagesc(riS_std)
colorbar
axis('square')
title('\sigma_{sim} [mm]')
caxis([150 350])
export_fig Figs\tif\Estimate_ann_rain.tif -tif -transparent -r300
savefig('Figs\Estimate_ann_rain.fig'); close all
%% Plot monthly rainfall at a point
% In this example data from a point location (e.g. Engelberg station) is
% compared to an areal rainfall data (2 x 2 km grid at a specific location).
% This comparison is simply done for a rough estimation of the model.
% obs:
ENGs=[89;89;108;113;155;178;196;190;130;101;108;103]; % Monthly data from Meteo Swiss climate normals (1981-2010)
% sim:
clear tmp
for m=1:12
    tmp(m,1)=riM{m}(6,4); % Closest location to the Engelberg station, valley
    tmp(m,2)=riM_std{m}(6,4);
end
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,2)
hold on
plot(1:12,ENGs(:,1),'.')
errorbar(1:12,tmp(:,1),tmp(:,2),'red.')
grid
box('on')
axis([0.5 12.5 0 300])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Rainfall [mm]')
legend('Obs.','Sim.')
title('Engelberg station [1,035 m a.s.l]')
export_fig Figs\tif\Monthly_rain_point.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Monthly_rain_at_a_point.fig'); close(gcf);
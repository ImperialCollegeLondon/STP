%Estimating_pressure This script run the pressure simulator for
%30 years and compute some statistics for the validation process
addpath(genpath(cd))
%% Load data
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
load('Zlevels.mat')
load('Shading.mat')
load('MeanArealTempParameters.mat')
load('ARcoeffT.mat')
load('genLRpar.mat')
load('LR.mat')
load('Pressure.mat')
load('Ps_Stations.mat')
%% Processing
%% Initilazing
domainDTM=rasterread('engelberger_DTM.txt');
domainDTM=flipud(domainDTM);
rng('shuffle');
td=(datenum(2001,1,1,0,0,0):datenum(2001,12,31,23,0,0))'; % Dates for one year based on the year 2001
I=1; % Number of realizations to simulate (the processing stage can be looped)
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
for I=1:1
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
        end
    end
    clear tsWARFIMA m
    simCAR(simCAR<simWAR)=simWAR(simCAR<simWAR);
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
                normClt(normClt>1)=.99999;
                normClt(normClt<0)=0.0001;
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
                normClt(normClt>1)=.99999;
                normClt(normClt<0)=0.0001;
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
                normClt(normClt>1)=.99999;
                normClt(normClt<0)=0.0001;
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
    end
    %% Advection - fixed for 5 m s^-1, westrlly winds
    simU=zeros(tsMatrix(end),1)+5;
    simV=zeros(tsMatrix(end),1);
    %% Converting simulated cloud cover from 5-min resolution to 1-h resolution
    Nsim=mean(reshape(simCAR,12,[]))';
    %% Generating the air temperature (2-m) field
    HZr=HZ*pi/180;
    t=(datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,23,0,0))'; % Hourly dates for one year based on the year 2001
    Tmax=zeros(length(Nsim),1); Tmin=zeros(length(Nsim),1);
    Ta=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2));
    mdl=arima('Constant',0,'Variance',1,'AR',{mean(ARcoeffT)}); % Hourly lapse rate time series
    v = simulate(mdl,length(Ta));
    v = normcdf(v,mean(v),std(v));
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
end
progressbar('Generating pressure ensemble') % Init single bar
for I=1:30
    %% Computing pressure at ground level
    Ps=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Pressure field at ground level
    for m=1:12 % Generating pressure at reference point
        p_sim(month(t)==m)=simulate(Pressure(m).ARMdl,sum(month(t)==m));
    end
    parfor i=1:length(Nsim)
        [ P ] = computePressureField( p_sim(i) , Ps_Stations(2).elevation , squeeze(Ta(i,:,:)) , domainDTM ); % Reference point for this example, Engelberg station [Ps_Stations(2)]
        Ps(i,:,:)=P;
    end
    %% Save data
    annualP{I}=squeeze(mean(Ps));
    ENG_vg{I}=Ps(:,71,112); % Engelberger virtual gauge
    LUZ_vg{I}=Ps(:,236,198); % "Luzern" virtual gauge; not the same location as Luzern is outside of the domain, but same elevation
    PIL_vg{I}=Ps(:,166,151); % "Pilatus" virtual gauge; not the same location as Pilatus is outside of the domain, but same elevation
    progressbar(I/30) % Update progress bar
end
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Plot annual maps
annualMap_mean=zeros(250,250);
annualMap_std=zeros(250,250);
for j=1:250*250
    tmp=zeros(30,1);
    for I=1:30
        tmp(I,1)=annualP{I}(j);
    end
    annualMap_mean(j)=mean(tmp);
    annualMap_std(j)=std(tmp);
end
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(annualMap_mean)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Mean pressure (30-years simulation)')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
subplot(1,2,2)
imagesc(annualMap_std)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Pressure standard deviation (30-years simulation)')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
export_fig Figs\tif\Annual_Pressure_One_Reaslization.tif -tif -transparent -r300
savefig(gcf,'Figs\Annual_Pressure_One_Reaslization.fig'); close(gcf);
%% Plot PDF
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I});
end
oGeng=Ps_Stations(2).Ps(:,1);
pts = (round(min([min(oGeng) min(vGeng)])):2:round(max([max(oGeng) max(vGeng)])));
[fveng,xieng] = ksdensity(vGeng,pts,'function','pdf');
[foeng,xieng] = ksdensity(oGeng,pts,'function','pdf');
vGluz=[];
for I=1:30
    vGluz=cat(1,vGluz,LUZ_vg{I});
end
oGluz=Ps_Stations(3).Ps(:,1);
pts = (round(min([min(oGluz) min(vGluz)])):2:round(max([max(oGluz) max(vGluz)])));
[fvluz,xiluz] = ksdensity(vGluz,pts,'function','pdf');
[foluz,xiluz] = ksdensity(oGluz,pts,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I});
end
oGpil=Ps_Stations(4).Ps(:,1);
pts = (round(min([min(oGpil) min(vGpil)])):2:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil] = ksdensity(vGpil,pts,'function','pdf');
[fopil,xipil] = ksdensity(oGpil,pts,'function','pdf');

figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,2)
hold on
bar(xieng,foeng,'red')
plot(xieng,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Engelberg (\mu_{obs}=',num2str(nanmean(oGeng),4),'; \mu_{sim}=',num2str(nanmean(vGeng),4),'; \sigma_{obs}=',num2str(nanstd(oGeng),2),'; \sigma_{sim}=',num2str(nanstd(vGeng),2),')'])
xlabel('Pressure [hPa]')
ylabel('Frequency')
legend('Obs.','Sim.')
box('on')
axis('square')
subplot(1,3,3)
hold on
bar(xipil,fopil,'red')
plot(xipil,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Pilatus^* (\mu_{obs}=',num2str(nanmean(oGpil),4),'; \mu_{sim}=',num2str(nanmean(vGpil),4),'; \sigma_{obs}=',num2str(nanstd(oGpil),2),'; \sigma_{sim}=',num2str(nanstd(vGpil),2),')'])
xlabel('Pressure [hPa]')
ylabel('Frequency')
box('on')
axis('square')
subplot(1,3,1)
hold on
bar(xiluz,foluz,'red')
plot(xiluz,fvluz,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Luzern^* (\mu_{obs}=',num2str(nanmean(oGluz),4),'; \mu_{sim}=',num2str(nanmean(vGluz),4),'; \sigma_{obs}=',num2str(nanstd(oGluz),2),'; \sigma_{sim}=',num2str(nanstd(vGluz),2),')'])
xlabel('Pressure [hPa]')
ylabel('Frequency')
box('on')
axis('square')
export_fig Figs\tif\Pressure_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\Pressure_PDF.fig'); close(gcf);
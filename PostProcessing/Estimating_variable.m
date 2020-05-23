%Estimating_variable This script run the model for
%30 years and compute statistics of the following variables: temperaure,
%vapor pressure, relative humidity, dew temperature and shortwave incoming
%radiation.
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
load('Vapor')
load('Radiation.mat')
load('TerrainRadiation.mat')
load('Wind.mat')
load('Advection.mat')
Advection.Dry=Dry;
Advection.Wet=Wet;
%% Processing
%% Initilazing
domainDTM=rasterread('engelberger_DTM.txt');
domainDTM=flipud(domainDTM);
rng('shuffle');
td=(datenum(2001,1,1,0,0,0):datenum(2001,12,31,23,0,0))'; % Dates for one year based on the year 2001
I=30; % Number of realizations to simulate (the processing stage can be looped)
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
progressbar('Generating 30 years ensemble') % Init single bar
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
    %% Simulating advection
    clear S Ss Sd Sw
    simU=zeros(tsMatrix(end),1);
    simV=zeros(tsMatrix(end),1);
    parfor i=1:size(tsMatrix,1)
        Sd{i}=vgxsim(Advection.Dry(tsMatrix(i,2)).mdl,20000);
        Sw{i}=vgxsim(Advection.Wet(tsMatrix(i,2)).mdl,20000);
    end
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
    end
    simU=single(simU);
    simV=single(simV);
    clear S Ss Sd Sw
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
    for i=1:length(Nsim)
        m=month(t(i));
        %% Vapor pressure Generator
        esat_s=611.*exp(17.27.*squeeze(Ta(i,:,:))./(237.3+squeeze(Ta(i,:,:))));
        if m~=M
            %% Generating monthly vapor pressure data
            a0=reshape(Vapor.A0{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
            a1=reshape(Vapor.A1{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
            a2=Vapor.A2{m};
            a3=Vapor.A3{m};
            dDem=reshape(Vapor.dDem{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
            rhodDe=reshape(Vapor.rhodDe{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
            rhodDe(rhodDe>0.99)=0.99;
            sigmadDe=reshape(Vapor.sigmadDe{m}(domainDTM),size(domainDTM,1),size(domainDTM,2));
            M=m;
        end
        if i<3
            [ ea(i,:,:) , dDe(i) ] = ComputeVapPressure2( esat_s , squeeze(Ta(i,:,:)) , 0 , 0 , 0 , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
        else
            [ ea(i,:,:) , dDe(i) ] = ComputeVapPressure2( esat_s , squeeze(Ta(i,:,:)) , Rsw_tm1 , Rsw_tm2 , dDe(i-1) ,a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
        end
        U=squeeze(ea(i,:,:))./esat_s;
        G_g=17.27.*squeeze(Ta(i,:,:))./(237.7+squeeze(Ta(i,:,:)))+log(U);
        Tdew=237.7.*G_g./(17.27-G_g);
        Tdew(isinf(Tdew))=NaN;
        %% Radiation Generator
        if Sun.altitude(i)>0 || i==1
            [ cos_ST ] = SolarIlluminationAngle( Slo_top , Sun.altitude(i) , Sun.azimuth(i) , Aspect );
            [ S ] = Shadow_Effect2( domainDTM , pi/2-Sun.altitude(i) , Sun.azimuth(i) , HZr , Z );
            [SB,SD,SAB1,SAB2,SAD1,SAD2,PARB,PARD]=ComputeRadiationForcings( Sun.altitude(i) , Sun.E0(i) , Radiation.Ro , Radiation.lwp(m) , Nsim(i) , domainDTM , Tdew , Radiation.Angstrom(m).beta , Radiation.Angstrom(m).alpha , Radiation.omega_A1 ,Radiation.omega_A2 , Radiation.uo , Radiation.un , Radiation.rho_g , cos_ST , S , SvF , Ct , Radiation.AODcorrection(m) , Radiation.AERONETz );
            SB(SB<0)=0;
            SD(SD<0)=0;
            Ws=SB+SD;
        else
            Rsw=zeros(size(domainDTM));
        end
        Rsw_tm2=Rsw_tm1;
        Rsw_tm1=Rsw;
        ssr(i,:,:)=Rsw;
    end
    %% Save data
    annualRsw{I}=squeeze(mean(ssr));
    ENG_vg{I,1}=ssr(:,71,112); % Engelberger virtual gauge
    TIT_vg{I,1}=ssr(:,17,128); % Titlis virtual gauge
    PIL_vg{I,1}=ssr(:,166,151); % "Pilatus" virtual gauge; not the same location as Pilatus is outside of the domain, but same elevation
    LUZ_vg{I,1}=ssr(:,236,198); % "Luzern" virtual gauge; not the same location as Luzern is outside of the domain, but same elevation
    clear Rsw
    for i=1:length(Nsim)
        Rsw(i,1)=nanmean(nanmean(ssr(i,:,:)));
    end
    clear ssr
    ENG_vg{I,2}=ea(:,71,112);
    TIT_vg{I,2}=ea(:,17,128);
    PIL_vg{I,2}=ea(:,166,151);
    LUZ_vg{I,2}=ea(:,236,198);
    %% Relative humidity
    hur=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Relative humidity [%] field at ground level
    hur=single(hur);
    for i=1:length(Nsim)
        esat_s=611.*exp(17.27.*squeeze(Ta(i,:,:))./(237.3+squeeze(Ta(i,:,:))));
        hur(i,:,:)=squeeze(ea(i,:,:))./esat_s;
    end
    ENG_vg{I,3}=hur(:,71,112);
    TIT_vg{I,3}=hur(:,17,128);
    PIL_vg{I,3}=hur(:,166,151);
    LUZ_vg{I,3}=hur(:,236,198);
    clear ea
    %% Dew point temperature
    Tdew=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Dew point temperature [°C] field at ground level
    Tdew=single(Tdew);
    for i=1:length(Nsim)
        G_g=17.27.*squeeze(Ta(i,:,:))./(237.7+squeeze(Ta(i,:,:)))+log(squeeze(hur(i,:,:)));
        Tdew(i,:,:)=237.7.*G_g./(17.27-G_g);
    end
    Tdew(isinf(Tdew))=NaN;
    ENG_vg{I,4}=Tdew(:,71,112);
    TIT_vg{I,4}=Tdew(:,17,128);
    PIL_vg{I,4}=Tdew(:,166,151);
    LUZ_vg{I,4}=Tdew(:,236,198);
    clear Tdew hur
    %% Save data
    annualT{I}=squeeze(mean(Ta));
    for m=1:12
        monthT{I,m}=squeeze(mean(Ta(month(t)==m,:,:)));
    end
    ENG_vg{I,5}=Ta(:,71,112);
    TIT_vg{I,5}=Ta(:,17,128);
    PIL_vg{I,5}=Ta(:,166,151);
    clear Ta
    %% Near surface wind speed
    G=sqrt(simU.^2+simV.^2); % Geostrophic wind speed [m s^-1]
    G=mean(reshape(G,[],8760))'; % Geostrophic wind speed for 1-h intervals [m s^-1]
    [ classP ] = findPasquillClass( Nsim , Rsw ); % Find potential Pasquill classes for each hour
    %% Iterations
    Ws=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Wind speed [m s^-1] field at near ground level (2-m elevation)
    Ws_no_correction=zeros(length(Nsim),size(domainDTM,1),size(domainDTM,2)); % Wind speed [m s^-1] field at near ground level (2-m elevation)
    tF=0.01:0.01:5;
    for i=1:length(Nsim)
        cP=classP{i}';
        clear FVi
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
                        nF = @(fv,vk,G,AU,BU,cF,z0) (vk.*G)./sqrt((log(0.3.*(fv./abs(cF))./z0)-BU).^2+AU.^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),AU,BU,mean2(cF),Z0(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=0; % Monin-Obukov lenght L
                        H=0.3.*(FV./abs(cF)); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
                    case 5 % Stable atmosphere
                        FV=nan(size(domainDTM));
                        MOL=123.*Z0.^0.3;
                        nF = @(fv,vk,G,cF,z0,MOL) (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1./(1+0.882.*(sqrt((vk.*fv)./(abs(cF).*MOL))))))./z0)-0.5+2.55.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2+(4.5+1.765.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),mean2(cF),Z0(j),MOL(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=123.*z0.^0.3; % Monin-Obukov lenght L
                        coeffU=(vk.*FV)./(abs(cF).*MOL);
                        bU=4+10.2.*sqrt(coeffU);
                        b_U=-4.5-7.65.*sqrt(coeffU);
                        aU=10;
                        a_U=-5.5+1.765.*sqrt(coeffU);
                        H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
                    case 6 % Very stable atmosphere
                        FV=nan(size(domainDTM));
                        MOL=26.*Z0.^0.17;
                        nF = @(fv,vk,G,cF,z0,MOL) (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1./(1+0.882.*(sqrt((vk.*fv)./(abs(cF).*MOL))))))./z0)-0.5+2.55.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2+(4.5+1.765.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),mean2(cF),Z0(j),MOL(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=26.*z0.^0.17; % Monin-Obukov lenght L
                        coeffU=(vk.*FV)./(abs(cF).*MOL);
                        bU=4+10.2.*sqrt(coeffU);
                        b_U=-4.5-7.65.*sqrt(coeffU);
                        aU=10;
                        a_U=-5.5+1.765.*sqrt(coeffU);
                        H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
                    case 1 % Extreme unstable atmosphere
                        FV=nan(size(domainDTM));
                        MOL=-11.4.*Z0.^0.1;
                        nF = @(fv,vk,G,cF,z0,MOL) (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))./z0)-10+(9.5./(1+0.027.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))).^2+(4.5./(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL)))).^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),mean2(cF),Z0(j),MOL(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=-11.4.*z0.^0.1; % Monin-Obukov lenght L
                        coeffU=(vk.*FV)./(abs(cF).*MOL);
                        bU=-34+38./(1+0.027.*sqrt(-coeffU));
                        b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                        aU=10./(1+1.581.*sqrt(-coeffU));
                        a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                        H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
                    case 2 % Very unstable atmosphere
                        FV=nan(size(domainDTM));
                        MOL=-26.*Z0.^0.17;
                        nF = @(fv,vk,G,cF,z0,MOL) (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))./z0)-10+(9.5./(1+0.027.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))).^2+(4.5./(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL)))).^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),mean2(cF),Z0(j),MOL(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=-26.*z0.^0.17; % Monin-Obukov lenght L
                        coeffU=(vk.*FV)./(abs(cF).*MOL);
                        bU=-34+38./(1+0.027.*sqrt(-coeffU));
                        b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                        aU=10./(1+1.581.*sqrt(-coeffU));
                        a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                        H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
                    case 3 % Unstable atmosphere
                        FV=nan(size(domainDTM));
                        MOL=-123.*Z0.^0.3;
                        nF = @(fv,vk,G,cF,z0,MOL) (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))./z0)-10+(9.5./(1+0.027.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))).^2+(4.5./(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL)))).^2)-fv;
                        for j=1:length(Z0)
                            fun = @(fv) nF(fv,vk,G(i),mean2(cF),Z0(j),MOL(j));
                            tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                            FV(z0==Z0(j))=tmp;
                        end
                        FVi{cP(cPi)}=FV;
                        MOL=-123.*z0.^0.3; % Monin-Obukov lenght L
                        coeffU=(vk.*FV)./(abs(cF).*MOL);
                        bU=-34+38./(1+0.027.*sqrt(-coeffU));
                        b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                        aU=10./(1+1.581.*sqrt(-coeffU));
                        a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                        H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                        u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                        v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                        cP(cPi,2)=nanmean(nanmean(sqrt(u.^2+v.^2)));
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
                    nF = @(fv,vk,G,AU,BU,cF,z0) (vk.*G)./sqrt((log(0.3.*(fv./abs(cF))./z0)-BU).^2+AU.^2)-fv;
                    for j=1:length(Z0)
                        fun = @(fv) nF(fv,vk,G(i),AU,BU,mean2(cF),Z0(j));
                        tmp=tF(min(abs(fun(tF)))==abs(fun(tF)));
                        FV(z0==Z0(j))=tmp;
                    end
                end
                MOL=0; % Monin-Obukov lenght L
                H=0.3.*(FV./abs(cF)); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
            case 5 % Stable atmosphere
                FV=FVi{classP{i}};
                MOL=123.*z0.^0.3; % Monin-Obukov lenght L
                coeffU=(vk.*FV)./(abs(cF).*MOL);
                bU=4+10.2.*sqrt(coeffU);
                b_U=-4.5-7.65.*sqrt(coeffU);
                aU=10;
                a_U=-5.5+1.765.*sqrt(coeffU);
                H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                AU=4.5+1.765.*sqrt(coeffU);
                BU=0.5-2.55.*sqrt(coeffU);
            case 6 % Very stable atmosphere
                FV=FVi{classP{i}};
                MOL=26.*z0.^0.17; % Monin-Obukov lenght L
                coeffU=(vk.*FV)./(abs(cF).*MOL);
                bU=4+10.2.*sqrt(coeffU);
                b_U=-4.5-7.65.*sqrt(coeffU);
                aU=10;
                a_U=-5.5+1.765.*sqrt(coeffU);
                H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                AU=4.5+1.765.*sqrt(coeffU);
                BU=0.5-2.55.*sqrt(coeffU);
            case 1 % Extreme unstable atmosphere
                FV=FVi{classP{i}};
                MOL=-11.4.*z0.^0.1; % Monin-Obukov lenght L
                coeffU=(vk.*FV)./(abs(cF).*MOL);
                bU=-34+38./(1+0.027.*sqrt(-coeffU));
                b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                aU=10./(1+1.581.*sqrt(-coeffU));
                a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                AU=4.5./(1+1.581.*sqrt(-coeffU));
                BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
            case 2 % Very unstable atmosphere
                FV=FVi{classP{i}};
                MOL=-26.*z0.^0.17; % Monin-Obukov lenght L
                coeffU=(vk.*FV)./(abs(cF).*MOL);
                bU=-34+38./(1+0.027.*sqrt(-coeffU));
                b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                aU=10./(1+1.581.*sqrt(-coeffU));
                a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                AU=4.5./(1+1.581.*sqrt(-coeffU));
                BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
            case 3 % Unstable atmosphere
                FV=FVi{classP{i}};
                MOL=-123.*z0.^0.3; % Monin-Obukov lenght L
                coeffU=(vk.*FV)./(abs(cF).*MOL);
                bU=-34+38./(1+0.027.*sqrt(-coeffU));
                b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                aU=10./(1+1.581.*sqrt(-coeffU));
                a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                AU=4.5./(1+1.581.*sqrt(-coeffU));
                BU=10-(9.5./(1+0.027.*sqrt(-coeffU)));
        end
        %% Transfrom to Cartesian coordinate system
        alphaWind=-atand(AU./(log(FV./(abs(cF).*z0))-BU)); % Angle between geostrophic wind vector and the axis of the friction velocity [deg]
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
        uRot=u.*cos(rotAngle)-v.*sin(rotAngle);
        vRot=u.*sin(rotAngle)+v.*cos(rotAngle);
        %% Correction of slope and curve
        [ W , tet ] = z2p( uRot , vRot );
        [ omegaS ] = findOmegaS( betaS , tet , xiS );
        Ww=1+gammaS.*omegaS+gammaC.*omegaC;
        Wt=Ww.*W;
        Ws(i,:,:)=Wt;
        %     tetD=-0.225.*betaP.*sind(2.*(xiS-tet)); % Required if wished to save the U and V wind components
        %     tetT=tet+tetD;
        %     [ U , V ] = p2z( tetT , Wt );
        clear uRot vRot
    end
    parfor i=1:length(Nsim) % Smooth the field to correct for the spatial roughness variability
        Ws(i,:,:)=smoothn(squeeze(Ws(i,:,:)),5);
    end
    %% Save data
    annualWs{I}=squeeze(mean(Ws));
    for m=1:12
        monthWs{I,m}=squeeze(mean(Ws(month(t)==m,:,:)));
    end
    tmp1=zeros(length(Nsim),1);
    tmp2=zeros(length(Nsim),1);
    tmp3=zeros(length(Nsim),1);
    tmp4=zeros(length(Nsim),1);
    for i=1:length(Nsim)
        tmp1(i,1)=mean2(squeeze(Ws(i,70:72,111:113)));
        tmp2(i,1)=mean2(squeeze(Ws(i,16:18,127:129)));
        tmp3(i,1)=mean2(squeeze(Ws(i,165:167,150:152)));
        tmp4(i,1)=mean2(squeeze(Ws(i,235:237,197:199)));
    end
    ENG_vg{I,6}=tmp1;
    TIT_vg{I,6}=tmp2;
    PIL_vg{I,6}=tmp3;
    LUZ_vg{I,6}=tmp4;
    clear Ws
    progressbar(I/30) % Update progress bar    
end
%% Load Engelberger stations data
load('Ta_Stations.mat') % Temperature data from ground stations for ENG, TIT and PIL locations
load('Rsw_Stations.mat') % Shortwave radiation and relative humidity data for stations
load('Ws_Stations.mat') % Near ground level (2-m) wind speed data for stations
load('SIS.mat') % Gridded (2-km) annual shortwave radiation mean
%% Plot temperature annual maps
annualMap_mean=zeros(250,250);
annualMap_std=zeros(250,250);
for j=1:250*250
    tmp=zeros(30,1);
    for I=1:30
        tmp(I,1)=annualT{I}(j);
    end
    annualMap_mean(j)=mean(tmp);
    annualMap_std(j)=std(tmp);
end
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(2,1,1)
imagesc(annualMap_mean)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Mean temperature (30-years simulation)' )
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25]);
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
subplot(2,1,2)
imagesc(annualMap_std)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Temperature standard deviation (30-years simulation)')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
export_fig Figs\tif\Annual_Temperature_One_Reaslization.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Annual_Temperature_One_Reaslization.fig'); close(gcf);
%% Plot temperature monthly maps
for m=1:12
    monthlyMap_mean{m}=zeros(250,250);
    monthlyMap_std{m}=zeros(250,250);
end
for m=1:12
    for j=1:250*250
        tmp=zeros(30,1);
        for I=1:30
            tmp(I,1)=monthT{I,m}(j);
        end
        monthlyMap_mean{m}(j)=mean(tmp);
        monthlyMap_std{m}(j)=std(tmp);
    end
end
figure('units','normalized','outerposition',[0 0 1 1])
for m=1:12
    subplot(3,4,m)
    imagesc(monthlyMap_mean{m})
    colorbar
    axis('square')
    set(gca,'YDir','normal');
    title(datestr(datenum(1990,m,1),'mmm'))
    set(gca,'XTick',[] );
    set(gca,'XTickLabel',[] );
    set(gca,'YTick',[] );
    set(gca,'YTickLabel',[] );
    caxis([-10 18])
end
export_fig Figs\tif\Month_Temperature_One_Reaslization.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Temperature_One_Reaslization.fig'); close(gcf);
figure('units','normalized','outerposition',[0 0 1 1])
for m=1:12
    subplot(3,4,m)
    imagesc(monthlyMap_std{m})
    colorbar
    axis('square')
    set(gca,'YDir','normal');
    title(datestr(datenum(1990,m,1),'mmm'))
    set(gca,'XTick',[] );
    set(gca,'XTickLabel',[] );
    set(gca,'YTick',[] );
    set(gca,'YTickLabel',[] );
end
export_fig Figs\tif\Month_Temperature_One_Reaslization_std.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Temperature_One_Reaslization_std.fig'); close(gcf);
%% Monthly average temperature at stations
for m=1:12
    tmp=[];
    for y=1993:2014
        tmp=cat(1,tmp,nanmean(Ta_Stations{1}{3}(month(Ta_Stations{1}{3}(:,2))==m&year(Ta_Stations{1}{3}(:,2))==y,1)));
    end
    TIT(m,1)=nanmean(tmp); % Obs
    TIT(m,2)=nanstd(tmp);
    tmp=[];
    for y=1993:2014
        tmp=cat(1,tmp,nanmean(Ta_Stations{2}{3}(month(Ta_Stations{2}{3}(:,2))==m&year(Ta_Stations{2}{3}(:,2))==y,1)));
    end
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=1993:2014
        tmp=cat(1,tmp,nanmean(Ta_Stations{4}{3}(month(Ta_Stations{4}{3}(:,2))==m&year(Ta_Stations{4}{3}(:,2))==y,1)));
    end
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    tit=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,nanmean(ENG_vg{I,5}(month(t)==m)));
        tit=cat(1,tit,nanmean(TIT_vg{I,5}(month(t)==m)));
        pil=cat(1,pil,nanmean(PIL_vg{I,5}(month(t)==m)));
    end
    TITs(m,1)=nanmean(tit);
    TITs(m,2)=nanstd(tit);
    PILs(m,1)=nanmean(pil);
    PILs(m,2)=nanstd(pil);
    ENGs(m,1)=nanmean(eng);
    ENGs(m,2)=nanstd(eng);
end
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,1)
box('on')
hold on
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
axis([0.5 12.5 -10 25])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Temperature [^oC]')
legend('Obs.','Sim.')
title('Engelberg station [1,035 m a.s.l]')
subplot(3,1,2)
box('on')
hold on
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
axis([0.5 12.5 -15 15])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Temperature [^oC]')
title('Pilatus^* station [2,106 m a.s.l]')
subplot(3,1,3)
box('on')
hold on
errorbar(1:12,TIT(:,1),TIT(:,2),'blue.')
errorbar(1:12,TITs(:,1),TITs(:,2),'r.')
grid
axis([0.5 12.5 -20 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Temperature [^oC]')
title('Titlis station [3,040 m a.s.l]')
export_fig Figs\tif\Month_Temperature_One_Reaslization_stations.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Temperature_One_Reaslization_stations.fig'); close(gcf);
%% Plot Temp. distribution
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I,5});
end
oGeng=Ta_Stations{4}{3}(:,1);
pts = (round(min([min(oGeng) min(vGeng)])):2:round(max([max(oGeng) max(vGeng)])));
[fveng,xieng] = ksdensity(vGeng,pts,'function','pdf');
[foeng,xieng] = ksdensity(oGeng,pts,'function','pdf');
vGtit=[];
for I=1:30
    vGtit=cat(1,vGtit,TIT_vg{I,5});
end
oGtit=Ta_Stations{1}{3}(:,1);
pts = (round(min([min(oGtit) min(vGtit)])):2:round(max([max(oGtit) max(vGtit)])));
[fvtit,xitit] = ksdensity(vGtit,pts,'function','pdf');
[fotit,xitit] = ksdensity(oGtit,pts,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I,5});
end
oGpil=Ta_Stations{2}{3}(:,1);
pts = (round(min([min(oGpil) min(vGpil)])):2:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil] = ksdensity(vGpil,pts,'function','pdf');
[fopil,xipil] = ksdensity(oGpil,pts,'function','pdf');
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
box('on')
hold on
bar(xieng,foeng,'red')
plot(xieng,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Engelberg (\mu_{obs}=',num2str(nanmean(oGeng),2),'; \mu_{sim}=',num2str(nanmean(vGeng),2),'; \sigma_{obs}=',num2str(nanstd(oGeng),2),'; \sigma_{sim}=',num2str(nanstd(vGeng),2),')'])
xlabel('Temperature [^oC]')
ylabel('Frequency')
xlim([-30 30])
axis('square')
subplot(1,3,2)
hold on
box('on')
bar(xipil,fopil,'red')
plot(xipil,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Pilatus^* (\mu_{obs}=',num2str(nanmean(oGpil),2),'; \mu_{sim}=',num2str(nanmean(vGpil),2),'; \sigma_{obs}=',num2str(nanstd(oGpil),2),'; \sigma_{sim}=',num2str(nanstd(vGpil),2),')'])
xlabel('Temperature [^oC]')
ylabel('Frequency')
legend('Observed','Simulated')
xlim([-30 30])
axis('square')
subplot(1,3,3)
box('on')
hold on
bar(xitit,fotit,'red')
plot(xitit,fvtit,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Titlis (\mu_{obs}=',num2str(nanmean(oGtit),2),'; \mu_{sim}=',num2str(nanmean(vGtit),2),'; \sigma_{obs}=',num2str(nanstd(oGtit),2),'; \sigma_{sim}=',num2str(nanstd(vGtit),2),')'])
xlabel('Temperature [^oC]')
ylabel('Frequency')
xlim([-30 30])
axis('square')
export_fig Figs\tif\Temp_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Temp_PDF.fig'); close(gcf);
%% Plot average temperature daily cycle
for h=0:23
    dailyCycle(h+1,1)=nanmean(Ta_Stations{4}{3}(hour(Ta_Stations{4}{3}(:,2))==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,ENG_vg{I,5}(hour(t)==h));
    end
    dailyCycle(h+1,2)=nanmean(tmp);
    dailyCycle(h+1,3)=nanmean(Ta_Stations{2}{3}(hour(Ta_Stations{2}{3}(:,2))==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,PIL_vg{I,5}(hour(t)==h));
    end
    dailyCycle(h+1,4)=nanmean(tmp);
    dailyCycle(h+1,5)=nanmean(Ta_Stations{1}{3}(hour(Ta_Stations{1}{3}(:,2))==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,TIT_vg{I,5}(hour(t)==h));
    end
    dailyCycle(h+1,6)=nanmean(tmp);
end
dailyCycle(25,:)=dailyCycle(1,:);
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
plot(0:24,dailyCycle(:,1))
hold on
plot(0:24,dailyCycle(:,2),'r')
title('Engelberg [1,035 m a.s.l]')
xlim([0 24])
axis('square')
xlabel('Hour')
ylabel('Temperature [^oC]')
set(gca,'XTick',[0 6 12 18 24] );
grid
subplot(1,3,2)
plot(0:24,dailyCycle(:,3))
hold on
plot(0:24,dailyCycle(:,4),'r')
title('Pilatus^* [2,106 m a.s.l]')
xlim([0 24])
axis('square')
xlabel('Hour')
ylabel('Temperature [^oC]')
legend('Observed','Simulated')
set(gca,'XTick',[0 6 12 18 24] );
grid
subplot(1,3,3)
plot(0:24,dailyCycle(:,5))
hold on
plot(0:24,dailyCycle(:,6),'r')
title('Titlis [3,040 m a.s.l]')
xlim([0 24])
axis('square')
xlabel('Hour')
ylabel('Temperature [^oC]')
set(gca,'XTick',[0 6 12 18 24] );
grid
export_fig Figs\tif\TemperatureDailyCycle.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\TemperatureDailyCycle.fig'); close(gcf);
%% Plot annual observed and simulated incoming shortwave radiation
year_mean=zeros(size(domainDTM,1),size(domainDTM,2));
for I=1:30
    year_mean=year_mean+annualRsw{I};
end
year_mean=year_mean./I;
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(2,1,1)
imagesc(SIS)
caxis([60 180])
colorbar
axis('square')
set(gca,'YDir','normal');
title('Observed global radiation [W m^{-2}], 2-km resolution')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
subplot(2,1,2)
imagesc(year_mean)
caxis([60 180])
colorbar
axis('square')
set(gca,'YDir','normal');
title('Simulated global radiation [W m^{-2}], 0.1-km resolution')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
export_fig Figs\tif\Rsw_Annual_Fields.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Rsw_Annual_Fields.fig'); close(gcf);
%% Monthly average incoming shortwave radiation at stations
for m=1:12
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(3).Rsw(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y)));
    end
    LUZ(m,1)=nanmean(tmp); % Obs
    LUZ(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(4).Rsw(month(Rsw_Stations(4).Time)==m&year(Rsw_Stations(4).Time)==y)));
    end
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(2).Rsw(month(Rsw_Stations(2).Time)==m&year(Rsw_Stations(2).Time)==y)));
    end
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    luz=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,mean(ENG_vg{I,1}(month(t)==m)));
        luz=cat(1,luz,mean(LUZ_vg{I,1}(month(t)==m)));
        pil=cat(1,pil,mean(PIL_vg{I,1}(month(t)==m)));
    end
    LUZs(m,1)=mean(luz);
    LUZs(m,2)=std(luz);
    PILs(m,1)=mean(pil);
    PILs(m,2)=std(pil);
    ENGs(m,1)=mean(eng);
    ENGs(m,2)=std(eng);
end
figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(3,1,2)
hold on
box('on')
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
axis([0.5 12.5 0 300])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Global radiation [W m^{-2}]')
title('Engelberg station [1,035 m a.s.l]')
subplot(3,1,3)
hold on
box('on')
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
axis([0.5 12.5 0 300])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Global radiation [W m^{-2}]')
title('Pilatus^* station [2,106 m a.s.l]')
subplot(3,1,1)
hold on
box('on')
errorbar(1:12,LUZ(:,1),LUZ(:,2),'blue.')
errorbar(1:12,LUZs(:,1),LUZs(:,2),'r.')
grid
axis([0.5 12.5 0 300])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Global radiation [W m^{-2}]')
title('Luzern^* station [454 m a.s.l]')
legend('Obs.','Sim.')
export_fig Figs\tif\Month_Rsw_One_Reaslization_stations.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Rsw_One_Reaslization_stations.fig'); close(gcf);
%% Plot Rsw distribution
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I,1});
end
oGeng=Rsw_Stations(2).Rsw;
[fveng,xieng1] = ksdensity(vGeng,'function','pdf');
[foeng,xieng2] = ksdensity(oGeng,'function','pdf');
vGluz=[];
for I=1:30
    vGluz=cat(1,vGluz,LUZ_vg{I,1});
end
oGluz=Rsw_Stations(3).Rsw;
pts = (round(min([min(oGluz) min(vGluz)])):100:round(max([max(oGluz) max(vGluz)])));
[fvluz,xiluz1] = ksdensity(vGluz,'function','pdf');
[foluz,xiluz2] = ksdensity(oGluz,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I,1});
end
oGpil=Rsw_Stations(4).Rsw;
pts = (round(min([min(oGpil) min(vGpil)])):100:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil1] = ksdensity(vGpil,'function','pdf');
[fopil,xipil2] = ksdensity(oGpil,'function','pdf');
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,2)
hold on
bar(xieng2,foeng,'red')
plot(xieng1,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Engelberg'];[];['\mu_{o}=',num2str(nanmean(oGeng),3),' \mu_{s}=',num2str(nanmean(vGeng),3),' \sigma_{o}=',num2str(nanstd(oGeng),3),' \sigma_{s}=',num2str(nanstd(vGeng),3)]})
xlabel('Global radiation [W m^{-2}]')
ylabel('Frequency')
xlim([0 1000])
axis('square')
box('on')
legend('Observed','Simulated')
subplot(1,3,3)
hold on
bar(xipil2,fopil,'red')
plot(xipil1,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Pilatus^*'],;[];['\mu_{o}=',num2str(nanmean(oGpil),3),' \mu_{s}=',num2str(nanmean(vGpil),3),' \sigma_{o}=',num2str(nanstd(oGpil),3),' \sigma_{s}=',num2str(nanstd(vGpil),3)]})
xlabel('Global radiation [W m^{-2}]')
ylabel('Frequency')
xlim([0 1000])
axis('square')
box('on')
subplot(1,3,1)
hold on
bar(xiluz2,foluz,'red')
plot(xiluz1,fvluz,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Luzern^*'];[];['\mu_{o}=',num2str(nanmean(oGluz),3),' \mu_{s}=',num2str(nanmean(vGluz),3),' \sigma_{o}=',num2str(nanstd(oGluz),3),' \sigma_{s}=',num2str(nanstd(vGluz),3)]})
xlabel('Global radiation [W m^{-2}]')
ylabel('Frequency')
xlim([0 1000])
axis('square')
box('on')
export_fig Figs\tif\Rsw_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Rsw_PDF.fig'); close(gcf);
%% Plot average Rsw daily cycle
clear dailyCycle
for h=0:23
    dailyCycle(h+1,1)=nanmean(Rsw_Stations(2).Rsw(hour(Rsw_Stations(2).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,ENG_vg{I,1}(hour(t)==h));
    end
    dailyCycle(h+1,2)=nanmean(tmp);
    dailyCycle(h+1,3)=nanmean(Rsw_Stations(4).Rsw(hour(Rsw_Stations(4).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,PIL_vg{I,1}(hour(t)==h));
    end
    dailyCycle(h+1,4)=nanmean(tmp);
    dailyCycle(h+1,5)=nanmean(Rsw_Stations(3).Rsw(hour(Rsw_Stations(3).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,LUZ_vg{I,1}(hour(t)==h));
    end
    dailyCycle(h+1,6)=nanmean(tmp);
end
dailyCycle(:,1)=circshift(dailyCycle(:,1),1);
dailyCycle(:,3)=circshift(dailyCycle(:,3),1);
dailyCycle(:,5)=circshift(dailyCycle(:,5),1);
dailyCycle(25,:)=dailyCycle(1,:);
figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,3,2)
plot(0:24,dailyCycle(:,1))
hold on
plot(0:24,dailyCycle(:,2),'r')
title('Engelberg [1,035 m a.s.l]')
axis([0 24 0 500])
axis('square')
xlabel('Hour')
ylabel('Global radiation [W m^{-2}]')
set(gca,'XTick',[0 6 12 18 24] );
legend('Obs.','Sim.')
box('on')
grid
subplot(1,3,3)
plot(0:24,dailyCycle(:,3))
hold on
plot(0:24,dailyCycle(:,4),'r')
title('Pilatus^* [2,106 m a.s.l]')
axis([0 24 0 500])
axis('square')
xlabel('Hour')
ylabel('Global radiation [W m^{-2}]')
set(gca,'XTick',[0 6 12 18 24] );
box('on')
grid
subplot(1,3,1)
plot(0:24,dailyCycle(:,5))
hold on
plot(0:24,dailyCycle(:,6),'r')
title('Luzern^* [454 m a.s.l]')
axis([0 24 0 500])
axis('square')
xlabel('Hour')
ylabel('Global radiation [W m^{-2}]')
set(gca,'XTick',[0 6 12 18 24] );
box('on')
grid
export_fig Figs\tif\RswDailyCycle.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\RswDailyCycle.fig'); close(gcf);
%% Generating vapor pressure for the observed data
for i=1:length(Rsw_Stations)
    esat_s=611.*exp(17.27.*Rsw_Stations(i).Ta./(237.3+Rsw_Stations(i).Ta));
    ea(:,i)=esat_s.*(Rsw_Stations(i).RH./100);
end
%% Dew point temperature
for i=1:length(Rsw_Stations)
    g_g(:,1)=17.27.*Rsw_Stations(i).Ta./(237.7+Rsw_Stations(i).Ta)+log(Rsw_Stations(i).RH./100); % Observed
    Tdew(i,:)=237.7.*g_g./(17.27-g_g); % Observed
    Tdew(isinf(Tdew))=NaN;
end
Tdew=Tdew';
%% Ploting monthly mean vapor pressure
for m=1:12
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(ea(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,3)));
    end
    LUZ(m,1)=nanmean(tmp); % Obs
    LUZ(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(ea(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,4)));
    end
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(ea(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,2)));
    end
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    luz=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,mean(ENG_vg{I,2}(month(t)==m)));
        luz=cat(1,luz,mean(LUZ_vg{I,2}(month(t)==m)));
        pil=cat(1,pil,mean(PIL_vg{I,2}(month(t)==m)));
    end
    LUZs(m,1)=mean(luz);
    LUZs(m,2)=std(luz);
    PILs(m,1)=mean(pil);
    PILs(m,2)=std(pil);
    ENGs(m,1)=mean(eng);
    ENGs(m,2)=std(eng);
end
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,2)
hold on
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 2000])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Vapor pressure [Pa]')
legend('Obs.','Sim.')
title('Engelberg station [1,035 m a.s.l]')
subplot(3,1,3)
hold on
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 2000])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Vapor pressure [Pa]')
title('Pilatus^* station [2,106 m a.s.l]')
subplot(3,1,1)
hold on
errorbar(1:12,LUZ(:,1),LUZ(:,2),'blue.')
errorbar(1:12,LUZs(:,1),LUZs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 2000])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Vapor pressure [Pa]')
title('Luzern^* station [454 m a.s.l]')
export_fig Figs\tif\Monthly_avg_vp.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Monthly_avg_vp.fig'); close(gcf);
%% Relative humidity daily cycle
clear dailyCycle
for h=0:23
    dailyCycle(h+1,1)=nanmean(Rsw_Stations(2).RH(hour(Rsw_Stations(2).Time)==h,1)./100);
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,ENG_vg{I,3}(hour(t)==h));
    end
    dailyCycle(h+1,2)=nanmean(tmp);
    dailyCycle(h+1,3)=nanmean(Rsw_Stations(4).RH(hour(Rsw_Stations(4).Time)==h,1)./100);
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,PIL_vg{I,3}(hour(t)==h));
    end
    dailyCycle(h+1,4)=nanmean(tmp);
    dailyCycle(h+1,5)=nanmean(Rsw_Stations(3).RH(hour(Rsw_Stations(3).Time)==h,1)./100);
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,LUZ_vg{I,3}(hour(t)==h));
    end
    dailyCycle(h+1,6)=nanmean(tmp);
end
dailyCycle(:,1)=circshift(dailyCycle(:,1),1);
dailyCycle(:,3)=circshift(dailyCycle(:,3),1);
dailyCycle(:,5)=circshift(dailyCycle(:,5),1);
dailyCycle(25,:)=dailyCycle(1,:);
figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,3,2)
plot(0:24,dailyCycle(:,1))
hold on
plot(0:24,dailyCycle(:,2),'r')
title('Engelberg [1,035 m a.s.l]')
axis([0 24 0.5 1])
axis('square')
box('on')
xlabel('Hour')
ylabel('Relative humidity [-]')
set(gca,'XTick',[0 6 12 18 24] );
grid
legend('Obs.','Sim.')
subplot(1,3,3)
plot(0:24,dailyCycle(:,3))
hold on
plot(0:24,dailyCycle(:,4),'r')
title('Pilatus^* [2,106 m a.s.l]')
axis([0 24 0.5 1])
axis('square')
box('on')
xlabel('Hour')
ylabel('Relative humidity [-]')
set(gca,'XTick',[0 6 12 18 24] );
grid
subplot(1,3,1)
plot(0:24,dailyCycle(:,5))
hold on
plot(0:24,dailyCycle(:,6),'r')
title('Luzern^* [454 m a.s.l]')
axis([0 24 0.5 1])
axis('square')
xlabel('Hour')
box('on')
ylabel('Relative humidity [-]')
set(gca,'XTick',[0 6 12 18 24] );
grid
export_fig Figs\tif\RH_Daily_Cycle.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\RH_Daily_Cycle.fig');  close(gcf);
%% Plot VP distribution
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I,2});
end
oGeng=ea(:,2);
pts = (round(min([min(oGeng) min(vGeng)])):200:round(max([max(oGeng) max(vGeng)])));
[fveng,xieng] = ksdensity(vGeng,pts,'function','pdf');
[foeng,xieng] = ksdensity(oGeng,pts,'function','pdf');
vGluz=[];
for I=1:30
    vGluz=cat(1,vGluz,LUZ_vg{I,2});
end
oGluz=ea(:,3);
pts = (round(min([min(oGluz) min(vGluz)])):200:round(max([max(oGluz) max(vGluz)])));
[fvluz,xiluz] = ksdensity(vGluz,pts,'function','pdf');
[foluz,xiluz] = ksdensity(oGluz,pts,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I,2});
end
oGpil=ea(:,4);
pts = (round(min([min(oGpil) min(vGpil)])):200:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil] = ksdensity(vGpil,pts,'function','pdf');
[fopil,xipil] = ksdensity(oGpil,pts,'function','pdf');
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,2)
hold on
bar(xieng,foeng,'red')
plot(xieng,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Engelberg (\mu_{obs}=',num2str(nanmean(oGeng),3),'; \mu_{sim}=',num2str(nanmean(vGeng),3),'; \sigma_{obs}=',num2str(nanstd(oGeng),3),'; \sigma_{sim}=',num2str(nanstd(vGeng),3),')'];[]})
xlabel('Vapor pressure [Pa]')
ylabel('Frequency')
xlim([0 3000])
box('on')
axis('square')
legend('Obs.','Sim.')
subplot(1,3,3)
hold on
bar(xipil,fopil,'red')
plot(xipil,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Pilatus^* (\mu_{obs}=',num2str(nanmean(oGpil),3),'; \mu_{sim}=',num2str(nanmean(vGpil),3),'; \sigma_{obs}=',num2str(nanstd(oGpil),3),'; \sigma_{sim}=',num2str(nanstd(vGpil),3),')'];[]})
xlabel('Vapor pressure [Pa]')
ylabel('Frequency')
xlim([0 3000])
box('on')
axis('square')
subplot(1,3,1)
hold on
bar(xiluz,foluz,'red')
plot(xiluz,fvluz,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title({['Luzern^* (\mu_{obs}=',num2str(nanmean(oGluz),3),'; \mu_{sim}=',num2str(nanmean(vGluz),3),'; \sigma_{obs}=',num2str(nanstd(oGluz),3),'; \sigma_{sim}=',num2str(nanstd(vGluz),3),')'];[]})
xlabel('Vapor pressure [Pa]')
ylabel('Frequency')
xlim([0 3000])
box('on')
axis('square')
export_fig Figs\tif\VP_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\VP_PDF.fig'); close(gcf);
%% Plot DP distribution
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I,4});
end
oGeng=Tdew(:,2);
pts = (round(min([min(oGeng) min(vGeng)])):2:round(max([max(oGeng) max(vGeng)])));
[fveng,xieng] = ksdensity(vGeng,pts,'function','pdf');
[foeng,xieng] = ksdensity(oGeng,pts,'function','pdf');
vGluz=[];
for I=1:30
    vGluz=cat(1,vGluz,LUZ_vg{I,4});
end
oGluz=Tdew(:,3);
pts = (round(min([min(oGluz) min(vGluz)])):2:round(max([max(oGluz) max(vGluz)])));
[fvluz,xiluz] = ksdensity(vGluz,pts,'function','pdf');
[foluz,xiluz] = ksdensity(oGluz,pts,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I,4});
end
oGpil=Tdew(:,4);
pts = (round(min([min(oGpil) min(vGpil)])):2:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil] = ksdensity(vGpil,pts,'function','pdf');
[fopil,xipil] = ksdensity(oGpil,pts,'function','pdf');
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,2)
hold on
bar(xieng,foeng,'red')
plot(xieng,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Engelberg (\mu_{obs}=',num2str(nanmean(oGeng),3),'; \mu_{sim}=',num2str(nanmean(vGeng),3),'; \sigma_{obs}=',num2str(nanstd(oGeng),3),'; \sigma_{sim}=',num2str(nanstd(vGeng),3),')'])
xlabel('Dew point temperature [^oC]')
ylabel('Frequency')
xlim([-30 30])
axis('square')
box('on')
legend('Obs.','Sim.')
subplot(1,3,3)
hold on
bar(xipil,fopil,'red')
plot(xipil,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Pilatus^* (\mu_{obs}=',num2str(nanmean(oGpil),3),'; \mu_{sim}=',num2str(nanmean(vGpil),3),'; \sigma_{obs}=',num2str(nanstd(oGpil),3),'; \sigma_{sim}=',num2str(nanstd(vGpil),3),')'])
xlabel('Dew point temperature [^oC]')
ylabel('Frequency')
xlim([-30 30])
axis('square')
box('on')
subplot(1,3,1)
hold on
bar(xiluz,foluz,'red')
plot(xiluz,fvluz,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Luzern^* (\mu_{obs}=',num2str(nanmean(oGluz),3),'; \mu_{sim}=',num2str(nanmean(vGluz),3),'; \sigma_{obs}=',num2str(nanstd(oGluz),3),'; \sigma_{sim}=',num2str(nanstd(vGluz),3),')'])
xlabel('Dew point temperature [^oC]')
ylabel('Frequency')
xlim([-30 30])
axis('square')
box('on')
export_fig Figs\tif\DP_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\DP_PDF.fig'); close(gcf);
%% Plot RH distribution
vGeng=[];
for I=1:30
    vGeng=cat(1,vGeng,ENG_vg{I,3});
end
oGeng=Rsw_Stations(2).RH./100;
pts = (round(min([min(oGeng) min(vGeng)])):0.05:round(max([max(oGeng) max(vGeng)])));
[fveng,xieng] = ksdensity(vGeng,pts,'function','pdf');
[foeng,xieng] = ksdensity(oGeng,pts,'function','pdf');
vGluz=[];
for I=1:30
    vGluz=cat(1,vGluz,LUZ_vg{I,3});
end
oGluz=Rsw_Stations(3).RH./100;
pts = (round(min([min(oGluz) min(vGluz)])):0.05:round(max([max(oGluz) max(vGluz)])));
[fvluz,xiluz] = ksdensity(vGluz,pts,'function','pdf');
[foluz,xiluz] = ksdensity(oGluz,pts,'function','pdf');
vGpil=[];
for I=1:30
    vGpil=cat(1,vGpil,PIL_vg{I,3});
end
oGpil=Rsw_Stations(4).RH./100;
pts = (round(min([min(oGpil) min(vGpil)])):0.05:round(max([max(oGpil) max(vGpil)])));
[fvpil,xipil] = ksdensity(vGpil,pts,'function','pdf');
[fopil,xipil] = ksdensity(oGpil,pts,'function','pdf');
figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,3,2)
hold on
bar(xieng,foeng,'red')
plot(xieng,fveng,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Engelberg (\mu_{obs}=',num2str(nanmean(oGeng),3),'; \mu_{sim}=',num2str(nanmean(vGeng),3),'; \sigma_{obs}=',num2str(nanstd(oGeng),3),'; \sigma_{sim}=',num2str(nanstd(vGeng),3),')'])
xlabel('Relative humidity [-]')
ylabel('Frequency')
xlim([-0.05 1.05])
axis('square')
box('on')
legend('Obs.','Sim.')
subplot(1,3,3)
hold on
bar(xipil,fopil,'red')
plot(xipil,fvpil,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Pilatus^* (\mu_{obs}=',num2str(nanmean(oGpil),3),'; \mu_{sim}=',num2str(nanmean(vGpil),3),'; \sigma_{obs}=',num2str(nanstd(oGpil),3),'; \sigma_{sim}=',num2str(nanstd(vGpil),3),')'])
xlabel('Relative humidity [-]')
ylabel('Frequency')
xlim([-0.05 1.05])
axis('square')
box('on')
subplot(1,3,1)
hold on
bar(xiluz,foluz,'red')
plot(xiluz,fvluz,'-blue','LineWidth',2,'MarkerSize',20,'Marker','.')
title(['Luzern^* (\mu_{obs}=',num2str(nanmean(oGluz),3),'; \mu_{sim}=',num2str(nanmean(vGluz),3),'; \sigma_{obs}=',num2str(nanstd(oGluz),3),'; \sigma_{sim}=',num2str(nanstd(vGluz),3),')'])
xlabel('Relative humidity [-]')
ylabel('Frequency')
xlim([-0.05 1.05])
axis('square')
box('on')
export_fig Figs\tif\RH_PDF.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\RH_PDF.fig'); close(gcf);
%% Ploting monthly mean RH
for m=1:12
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(3).RH(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,1)));
    end
    tmp=tmp./100;
    LUZ(m,1)=nanmean(tmp); % Obs
    LUZ(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(4).RH(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,1)));
    end
    tmp=tmp./100;
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Rsw_Stations(2).RH(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,1)));
    end
    tmp=tmp./100;
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    luz=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,mean(ENG_vg{I,3}(month(t)==m)));
        luz=cat(1,luz,mean(LUZ_vg{I,3}(month(t)==m)));
        pil=cat(1,pil,mean(PIL_vg{I,3}(month(t)==m)));
    end
    LUZs(m,1)=mean(luz);
    LUZs(m,2)=std(luz);
    PILs(m,1)=mean(pil);
    PILs(m,2)=std(pil);
    ENGs(m,1)=mean(eng);
    ENGs(m,2)=std(eng);
end
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,2)
hold on
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0.5 1])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Relative humidity [-]')
title('Engelberg station [1,035 m a.s.l]')
subplot(3,1,3)
hold on
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0.5 1])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Relative humidity [-]')
title('Pilatus^* station [2,106 m a.s.l]')
legend('Obs.','Sim.')
subplot(3,1,1)
hold on
errorbar(1:12,LUZ(:,1),LUZ(:,2),'blue.')
errorbar(1:12,LUZs(:,1),LUZs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0.5 1])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Relative humidity [-]')
title('Luzern^* station [454 m a.s.l]')
export_fig Figs\tif\Monthly_avg_rh.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Monthly_avg_rh.fig'); close(gcf);
%% Ploting monthly mean DP
for m=1:12
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Tdew(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,3)));
    end
    LUZ(m,1)=nanmean(tmp); % Obs
    LUZ(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Tdew(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,4)));
    end
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=2000:2014
        tmp=cat(1,tmp,nanmean(Tdew(month(Rsw_Stations(3).Time)==m&year(Rsw_Stations(3).Time)==y,2)));
    end
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    luz=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,mean(ENG_vg{I,4}(month(t)==m)));
        luz=cat(1,luz,mean(LUZ_vg{I,4}(month(t)==m)));
        pil=cat(1,pil,mean(PIL_vg{I,4}(month(t)==m)));
    end
    LUZs(m,1)=mean(luz);
    LUZs(m,2)=std(luz);
    PILs(m,1)=mean(pil);
    PILs(m,2)=std(pil);
    ENGs(m,1)=mean(eng);
    ENGs(m,2)=std(eng);
end
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,2)
hold on
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 -20 20])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Dew point [^oC]')
title('Engelberg station [1,035 m a.s.l]')
legend('Obs.','Sim.')
subplot(3,1,3)
hold on
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 -20 20])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Dew point [^oC]')
title('Pilatus^* station [2,106 m a.s.l]')
subplot(3,1,1)
hold on
errorbar(1:12,LUZ(:,1),LUZ(:,2),'blue.')
errorbar(1:12,LUZs(:,1),LUZs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 -20 20])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Dew point [^oC]')
title('Luzern^* station [454 m a.s.l]')
export_fig Figs\tif\Monthly_avg_Tdew.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Monthly_avg_Tdew.fig'); close(gcf);
%% Plot near surface wind annual maps
annualMap_mean=zeros(250,250);
annualMap_std=zeros(250,250);
for j=1:250*250
    tmp=zeros(30,1);
    for I=1:30
        tmp(I,1)=annualWs{I}(j);
    end
    annualMap_mean(j)=mean(tmp);
    annualMap_std(j)=std(tmp);
end
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(2,1,1)
imagesc(annualMap_mean)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Mean sim. near ground wind speed' )
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25]);
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
subplot(2,1,2)
imagesc(annualMap_std)
colorbar
axis('square')
set(gca,'YDir','normal');
title('Sim. near ground wind speed std.')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
export_fig Figs\tif\Annual_Ws_One_Reaslization.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Annual_Ws_One_Reaslization.fig'); close(gcf);
%% Plot wind speed (near ground) monthly maps
for m=1:12
    monthlyMap_mean{m}=zeros(250,250);
    monthlyMap_std{m}=zeros(250,250);
end
for m=1:12
    for j=1:250*250
        tmp=zeros(30,1);
        for I=1:30
            tmp(I,1)=monthWs{I,m}(j);
        end
        monthlyMap_mean{m}(j)=mean(tmp);
        monthlyMap_std{m}(j)=std(tmp);
    end
end
figure('units','normalized','outerposition',[0 0 1 1])
for m=1:12
    subplot(3,4,m)
    imagesc(monthlyMap_mean{m})
    colorbar
    axis('square')
    set(gca,'YDir','normal');
    title(datestr(datenum(1990,m,1),'mmm'))
    set(gca,'XTick',[] );
    set(gca,'XTickLabel',[] );
    set(gca,'YTick',[] );
    set(gca,'YTickLabel',[] );
    caxis([0 10])
end
export_fig Figs\tif\Month_Ws_One_Reaslization.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Ws_One_Reaslization.fig'); close(gcf);
figure('units','normalized','outerposition',[0 0 1 1])
for m=1:12
    subplot(3,4,m)
    imagesc(monthlyMap_std{m})
    colorbar
    axis('square')
    set(gca,'YDir','normal');
    title(datestr(datenum(1990,m,1),'mmm'))
    set(gca,'XTick',[] );
    set(gca,'XTickLabel',[] );
    set(gca,'YTick',[] );
    set(gca,'YTickLabel',[] );
end
export_fig Figs\tif\Month_Ws_One_Reaslization_std.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Month_Ws_One_Reaslization_std.fig'); close(gcf);
%% Wind speed daily cycle
clear dailyCycle
for h=0:23
    dailyCycle(h+1,1)=nanmean(Ws_Stations(2).Ws(hour(Ws_Stations(2).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,ENG_vg{I,6}(hour(t)==h));
    end
    dailyCycle(h+1,2)=nanmean(tmp);
    dailyCycle(h+1,3)=nanmean(Ws_Stations(4).Ws(hour(Ws_Stations(4).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,PIL_vg{I,6}(hour(t)==h));
    end
    dailyCycle(h+1,4)=nanmean(tmp);
    dailyCycle(h+1,5)=nanmean(Ws_Stations(1).Ws(hour(Ws_Stations(1).Time)==h,1));
    tmp=[];
    for I=1:30
        tmp=cat(1,tmp,TIT_vg{I,6}(hour(t)==h));
    end
    dailyCycle(h+1,6)=nanmean(tmp);
end
dailyCycle(:,1)=circshift(dailyCycle(:,1),1);
dailyCycle(:,3)=circshift(dailyCycle(:,3),1);
dailyCycle(:,5)=circshift(dailyCycle(:,5),1);
dailyCycle(25,:)=dailyCycle(1,:);
figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,3,1)
plot(0:24,dailyCycle(:,1))
hold on
plot(0:24,dailyCycle(:,2),'r')
title('Engelberg [1,035 m a.s.l]')
axis([0 24 0 10])
axis('square')
box('on')
xlabel('Hour')
ylabel('Wind speed [m s^{-1}]')
set(gca,'XTick',[0 6 12 18 24] );
grid
subplot(1,3,2)
plot(0:24,dailyCycle(:,3))
hold on
plot(0:24,dailyCycle(:,4),'r')
title('Pilatus^* [2,106 m a.s.l]')
axis([0 24 0 10])
axis('square')
box('on')
xlabel('Hour')
ylabel('Wind speed [m s^{-1}]')
set(gca,'XTick',[0 6 12 18 24] );
grid
legend('Obs.','Sim.')
subplot(1,3,3)
plot(0:24,dailyCycle(:,5))
hold on
plot(0:24,dailyCycle(:,6),'r')
title('Titlis [3,040 m a.s.l]')
axis([0 24 0 10])
axis('square')
xlabel('Hour')
box('on')
ylabel('Wind speed [m s^{-1}]')
set(gca,'XTick',[0 6 12 18 24] );
grid
export_fig Figs\tif\Ws_Daily_Cycle.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Ws_Daily_Cycle.fig');  close(gcf);
%% Ploting monthly mean wind speed
for m=1:12
    tmp=[];
    for y=1999:2012
        tmp=cat(1,tmp,nanmean(Ws_Stations(1).Ws(month(Ws_Stations(1).Time)==m&year(Ws_Stations(1).Time)==y,1)));
    end
    tmp=tmp;
    TIT(m,1)=nanmean(tmp); % Obs
    TIT(m,2)=nanstd(tmp);
    tmp=[];
    for y=1999:2012
        tmp=cat(1,tmp,nanmean(Ws_Stations(4).Ws(month(Ws_Stations(3).Time)==m&year(Ws_Stations(3).Time)==y,1)));
    end
    tmp=tmp;
    PIL(m,1)=nanmean(tmp);
    PIL(m,2)=nanstd(tmp);
    tmp=[];
    for y=1999:2012
        tmp=cat(1,tmp,nanmean(Ws_Stations(2).Ws(month(Ws_Stations(3).Time)==m&year(Ws_Stations(3).Time)==y,1)));
    end
    tmp=tmp;
    ENG(m,1)=nanmean(tmp);
    ENG(m,2)=nanstd(tmp);
    tit=[]; pil=[]; eng=[]; % Sim
    for I=1:30
        eng=cat(1,eng,mean(ENG_vg{I,6}(month(t)==m)));
        tit=cat(1,tit,mean(TIT_vg{I,6}(month(t)==m)));
        pil=cat(1,pil,mean(PIL_vg{I,6}(month(t)==m)));
    end
    TITs(m,1)=mean(tit);
    TITs(m,2)=std(tit);
    PILs(m,1)=mean(pil);
    PILs(m,2)=std(pil);
    ENGs(m,1)=mean(eng);
    ENGs(m,2)=std(eng);
end
figure('units','normalized','outerposition',[0 0 .5 1])
subplot(3,1,1)
hold on
errorbar(1:12,ENG(:,1),ENG(:,2),'blue.')
errorbar(1:12,ENGs(:,1),ENGs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Wind speed [m s^{-1}]')
title('Engelberg station [1,035 m a.s.l]')
subplot(3,1,2)
hold on
errorbar(1:12,PIL(:,1),PIL(:,2),'blue.')
errorbar(1:12,PILs(:,1),PILs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Wind speed [m s^{-1}]')
title('Pilatus^* station [2,106 m a.s.l]')
legend('Obs.','Sim.')
subplot(3,1,3)
hold on
errorbar(1:12,TIT(:,1),TIT(:,2),'blue.')
errorbar(1:12,TITs(:,1),TITs(:,2),'r.')
grid
box('on')
axis([0.5 12.5 0 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('Jan');('Feb');('Mar');('Apr');('May');('Jun');('Jul');('Aug');('Sep');('Oct');('Nov');('Dec')] );
ylabel('Wind speed [m s^{-1}]')
title('Titlis station [3,040 m a.s.l]')
export_fig Figs\tif\Monthly_avg_Ws.tif -tif -transparent -r300
savefig(gcf,'Figs\matlab\Monthly_avg_Ws.fig'); close(gcf);
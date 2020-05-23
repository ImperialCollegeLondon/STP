%Calibrate z0 for the wind component
addpath(genpath(cd))
%% Load data
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
load('Precipitation.mat');
Precipitation.data=data;
Precipitation.WAR=WAR;
Precipitation.IMF=IMF;
load('Ws_Stations_monthly.mat')
%% Initilazing
domainDTM=rasterread('engelberger_DTM.txt');
domainDTM=flipud(domainDTM);
rng('shuffle');
td=(datenum(2001,1,1,0,0,0):datenum(2001,12,31,23,0,0))'; % Dates for one year based on the year 2001
%% Processing
for R=1
    I=30; % Number of years to simulate (the processing stage can be looped)
    %% Storm arrival process
    [ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
    progressbar('Generating 30 years')
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
        % Save
        tsMatrix(:,5)=tsMatrix(:,4)-tsMatrix(:,3)+1;
        for m=1:12
            eventSum{R}{I,m}=sum(tsMatrix(tsMatrix(:,2)==m,1));
            wetSum{R}{I,m}=sum(tsMatrix(tsMatrix(:,1)==1&tsMatrix(:,2)==m,5))*5;
        end
        %% Varfima - simulating WAR, CAR and IMF
        simWAR=zeros(105120,1);
        simIMF=zeros(105120,1);
        simCAR=zeros(105120,1);
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
        end
        % Save
        sCAR{R}{I}=simCAR;
        sWAR{R}{I}=simWAR;
        sIMF{R}{I}=simIMF;
        %% Simulating advection
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
                        E=find(min((Sd{i}(1:5000,1)-YS(1)).^2+(Sd{i}(1:5000,2)-YS(2)).^2)==(Sd{i}(1:5000,1)-YS(1)).^2+(Sd{i}(1:5000,2)-YS(2)).^2);
                        E=E(1);
                        S=Sd{i}(E:E+tsMatrix(i,4)-tsMatrix(i,3),:);
                        YS=S(end,:);
                        S=normcdf(S);
                        S(:,1)=norminv(S(:,1),Advection.Dry(m).unorm(1),Advection.Dry(m).unorm(2));
                        S(:,2)=norminv(S(:,2),Advection.Dry(m).vnorm(1),Advection.Dry(m).vnorm(2));
                    case 1
                        E=find(min((Sw{i}(1:5000,1)-YS(1)).^2+(Sw{i}(1:5000,2)-YS(2)).^2)==(Sw{i}(1:5000,1)-YS(1)).^2+(Sw{i}(1:5000,2)-YS(2)).^2);
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
        clear S Ss
        % Save
        sU{R}{I}=simU;
        sV{R}{I}=simV;
        %% Converting simulated cloud cover from 5-min resolution to 1-h resolution
        Nsim=mean(reshape(simCAR,12,[]))';
        %% Generating near surface air temperature (2-m) fields
        HZr=HZ*pi/180;
        t=(datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,23,0,0))'; % Hourly dates for one year based on the year 2001
        Tmax=zeros(length(Nsim),1); Tmin=zeros(length(Nsim),1);
        Ta=zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim));
        mdl=arima('Constant',0,'Variance',1,'AR',{mean(ARcoeffT)}); % Hourly lapse rate time series
        v = simulate(mdl,length(Ta));
        v = normcdf(v,mean(v),std(v));
        for i=1:length(Nsim)
            m=month(t(i));
            h=hour(t(i));
            GenLR=norminv(v(i),genLRpar(m).h(h+1).miuA.*(domainDTM./1000)+genLRpar(m).h(h+1).miuB,genLRpar(m).h(h+1).sigma); % Lapse rate correction
            ShF=ShadowEffect(domainDTM,Sun.altitude(i),Sun.azimuth(i),HZr,Z);
            if i==1
                [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},zeros(size(domainDTM)),dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),0,ones(size(domainDTM))*Ti(month(t(i))),ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
                Ta(:,:,i)=Tsim+GenLR;
                TI=ones(size(domainDTM))*Ti(month(t(i)));
            else
                if hour(t(i))==0
                    TI=T_tilde;
                    I2=zeros(size(domainDTM));
                    I3=zeros(size(domainDTM));
                    I4=zeros(size(domainDTM));
                end
                [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},qt,dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),dT,TI,I2,I3,I4,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
                Ta(:,:,i)=Tsim+GenLR;
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
                tmp=squeeze(Ta(:,:,i));
                tmp(domainDTM<=Z1)=tmp(domainDTM<=Z1)-invLR(i).*((Z1(domainDTM<=Z1)-domainDTM(domainDTM<=Z1))./1000);
                Ta(:,:,i)=tmp;
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
        ea=zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)); % Vapor pressure [Pa] field at ground level
        ea=single(ea);
        ssr=zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)); % Surface shortwave incoming radaiation [W m^-2] field at ground level
        ssr=single(ssr);
        Rsw=zeros(size(domainDTM));
        for i=1:length(Nsim)
            Ta_i=squeeze(Ta(:,:,i));
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
                [ ea(:,:,i) , dDe(i) ] = ComputeVapPressure2( esat_s , Ta_i , 0 , 0 , 0 , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
            else
                [ ea(:,:,i) , dDe(i) ] = ComputeVapPressure2( esat_s , Ta_i , Rsw_tm1 , Rsw_tm2 , dDe(i-1) ,a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
            end
            U=squeeze(ea(:,:,i))./esat_s;
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
            ssr(:,:,i)=Rsw;
        end
        clear Rsw
        for i=1:length(Nsim)
            Rsw(i,1)=nanmean(nanmean(ssr(:,:,i)));
        end
        clear ssr
        clear ea
        clear hur Tdew Ta 
        %% Near surface wind speed
        uRot=single(zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)));
        vRot=single(zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)));
        mcF=mean2(cF);
        G=sqrt(simU.^2+simV.^2); % Geostrophic wind speed [m s^-1]
        G=mean(reshape(G,[],8760))'; % Geostrophic wind speed for 1-h intervals [m s^-1]
        [ classP ] = findPasquillClass( Nsim , Rsw ); % Find potential Pasquill classes for each hour
        % For the calibration of the rougnhess:
        Z0=linspace(0.00001,2,50)';
        z0=reshape(repmat(Z0,1250,1),size(domainDTM,1),size(domainDTM,2));
        for i=1:length(Z0)
            z_cond{i}=(z0==Z0(i));
        end
        %%  Monin-Obukov lenght L
        MOL1=-11.4.*z0.^0.1;
        MOL2=-26.*z0.^0.17;
        MOL3=-123.*z0.^0.3;
        MOL5=123.*z0.^0.3;
        MOL6=26.*z0.^0.17;
        %% Iterations
        Ws=zeros(size(domainDTM,1),size(domainDTM,2),length(Nsim)); % Wind speed [m s^-1] field at near ground level (2-m elevation)
        tF=0.01:0.01:5;
        for i=1:length(Nsim)
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
                                [ fun ] = nF4(tF,vk,G(i),AU,BU,mcF,Z0(j));
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            MOL=0; % Monin-Obukov lenght L
                            H=0.3.*(FV./abs(cF)); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                            cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                        case 5 % Stable atmosphere
                            FV=nan(size(domainDTM));
                            MOL=123.*Z0.^0.3;
                            for j=1:length(Z0)
                                [ fun ] = nF5( tF,vk,G(i),mcF,z0(j),MOL(j) );
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            coeffU=(vk.*FV)./(abs(cF).*MOL5);
                            bU=4+10.2.*sqrt(coeffU);
                            b_U=-4.5-7.65.*sqrt(coeffU);
                            aU=10;
                            a_U=-5.5+1.765.*sqrt(coeffU);
                            H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                            cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                        case 6 % Very stable atmosphere
                            FV=nan(size(domainDTM));
                            MOL=26.*Z0.^0.17;
                            for j=1:length(Z0)
                                fun = nF6(tF,vk,G(i),mcF,Z0(j),MOL(j));
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            coeffU=(vk.*FV)./(abs(cF).*MOL6);
                            bU=4+10.2.*sqrt(coeffU);
                            b_U=-4.5-7.65.*sqrt(coeffU);
                            aU=10;
                            a_U=-5.5+1.765.*sqrt(coeffU);
                            H=0.3.*(FV./abs(cF)).*(1./(1+0.882.*sqrt(coeffU))); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                            cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                        case 1 % Extreme unstable atmosphere
                            FV=nan(size(domainDTM));
                            MOL=-11.4.*Z0.^0.1;
                            for j=1:length(Z0)
                                fun = nF1(tF,vk,G(i),mcF,Z0(j),MOL(j));
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            coeffU=(vk.*FV)./(abs(cF).*MOL1);
                            bU=-34+38./(1+0.027.*sqrt(-coeffU));
                            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                            aU=10./(1+1.581.*sqrt(-coeffU));
                            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                            cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                        case 2 % Very unstable atmosphere
                            FV=nan(size(domainDTM));
                            MOL=-26.*Z0.^0.17;
                            for j=1:length(Z0)
                                fun = nF2(tF,vk,G(i),mcF,Z0(j),MOL(j));
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            coeffU=(vk.*FV)./(abs(cF).*MOL2);
                            bU=-34+38./(1+0.027.*sqrt(-coeffU));
                            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                            aU=10./(1+1.581.*sqrt(-coeffU));
                            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                            cP(cPi,2)=nanmean(nanmean(sqrt(u.*u+v.*v)));
                        case 3 % Unstable atmosphere
                            FV=nan(size(domainDTM));
                            MOL=-123.*Z0.^0.3;
                            for j=1:length(Z0)
                                fun = nF3(tF,vk,G(i),mcF,Z0(j),MOL(j));
                                tmp=tF(min(abs(fun))==abs(fun));
                                FV(z_cond{j})=tmp(1);
                            end
                            FVi{cP(cPi)}=FV;
                            coeffU=(vk.*FV)./(abs(cF).*MOL3);
                            bU=-34+38./(1+0.027.*sqrt(-coeffU));
                            b_U=24-28.5./(1+0.027.*sqrt(-coeffU));
                            aU=10./(1+1.581.*sqrt(-coeffU));
                            a_U=-5.5./(1+1.581.*sqrt(-coeffU));
                            H=0.3.*(FV./abs(cF)).*(1+1.581.*sqrt(-coeffU)); % PBL height h
                            u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                            v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
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
                            fun = nF4(tF,vk,G(i),AU,BU,mcF,Z0(j));
                            tmp=tF(min(abs(fun))==abs(fun));
                            FV(z_cond{j})=tmp(1);
                        end
                    end
                    MOL=0; % Monin-Obukov lenght L
                    H=0.3.*(FV./abs(cF)); % PBL height h
                    u=(FV./vk).*(log(Hagl./z0)+bU.*((Hagl-z0)./H)+b_U.*((Hagl-z0)./H).^2);
                    v=(-FV./vk).*(aU.*((Hagl-z0)./H)+a_U.*((Hagl-z0)./H).^2);
                case 5 % Stable atmosphere
                    FV=FVi{classP{i}};
                    coeffU=(vk.*FV)./(abs(cF).*MOL5);
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
                    coeffU=(vk.*FV)./(abs(cF).*MOL6);
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
                    coeffU=(vk.*FV)./(abs(cF).*MOL1);
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
                    coeffU=(vk.*FV)./(abs(cF).*MOL2);
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
                    coeffU=(vk.*FV)./(abs(cF).*MOL3);
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
            uRot(:,:,i)=u.*cos(rotAngle)-v.*sin(rotAngle);
            vRot(:,:,i)=u.*sin(rotAngle)+v.*cos(rotAngle);
        end
        for i=1:50
            for m=1:12
                Y{i}(I,m)=mean(sqrt(squeeze(uRot(i,1,month(t)==m)).^2+squeeze(vRot(i,1,month(t)==m)).^2));
            end
        end
        progressbar(I/30)
    end
end
%% Seasonality fitting roughness (z0 - y) to mean monthly wind speed (x) derived from geostrophics wind ("simulated")
clear x y
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
for m=1:12
    for i=1:50
        x(i,1)=Z0(i);
        y(i,1)=mean(Y{i}(:,m));
    end
    x2=y(1:49);
    y2=x(1:49);
    w2z{m}=fit(x2,y2,'exp2',opts); % fit between the mean monthly wind at 2 m and the expected rougness (for the 500 hPa layer)
end
%% Seasonality fitting roughness (z0 - x) to mean monthly wind speed (y) derived from stations ("observed")
clear x y
for m=1:12
    for i=1:6
        x(i,1)=Ws_Stations_monthly(i).z0;
        y(i,1)=Ws_Stations_monthly(i).Wind(m);
    end
    z2w{m}=fit(x,y,'power2'); 
end
%% Save
save('Files\z0_Corr.mat','z2w','w2z','-v7.3');
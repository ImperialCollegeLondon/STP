%Generate_base_Ta This script run the temperature simulator for
%30 years and and save a generated base temperature (elevation==0 m) for
%lapse rate computation
%% If this script is called outside of the main body
% addpath(genpath(cd))
% % Load data
% load('Precipitation.mat');
% Precipitation.data=data;
% Precipitation.WAR=WAR;
% Precipitation.IMF=IMF;
% load('triNormTrans.mat')
% load('Cloud.mat')
% Cloud=load('Cloud.mat','ar');
% Cloud.cloudDis=cloudDis;
% Cloud.Dry=Dry;
% Cloud.Wet=Wet;
% load('Shading.mat')
% load('MeanArealTempParameters.mat')
% load('MeanArealTemp.mat')
%% Processing
%% Initilazing
domainDTM=rasterread('engelberger_DTM.txt');
domainDTM=flipud(domainDTM);
rng('shuffle');
td=(datenum(2001,1,1,0,0,0):datenum(2001,12,31,23,0,0))'; % Dates for one year based on the year 2001
I=30; % Number of realizations to simulate (the processing stage can be looped)
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
progressbar('Generating base level temperature (30 years ensemble)') % Init single bar
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
    %% Advection - fixed for 5 m s^-1, westrlly winds
    simU=zeros(tsMatrix(end),1)+5;
    simV=zeros(tsMatrix(end),1);
    %% Converting simulated cloud cover from 5-min resolution to 1-h resolution
    Nsim=mean(reshape(simCAR,12,[]))';
    %% Generating the air temperature (2-m) field
    HZr=HZ*pi/180;
    t=(datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,23,0,0))'; % Hourly dates for one year based on the year 2001
    Ta{I}=zeros(length(Nsim),1);
    for i=1:length(Nsim)
        ShF=ShadowEffect(domainDTM,Sun.altitude(i),Sun.azimuth(i),HZr,Z);
        if i==1
            [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},zeros(size(domainDTM)),dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),0,ones(size(domainDTM))*Ti(month(t(i))),ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*0,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
            Ta{I}(i,1)=mean2(Tsim);
            TI=ones(size(domainDTM))*Ti(month(t(i)));
        else
            if hour(t(i))==0
                TI=T_tilde;
                I2=zeros(size(domainDTM));
                I3=zeros(size(domainDTM));
                I4=zeros(size(domainDTM));
            end
            [Tsim,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(t(i),1,8.4,46.8,ones(size(domainDTM))*Nsim(i),Tb{month(t(i))},qt,dTbar(month(t(i)),1+hour(t(i))),dTrho(month(t(i))),dTsigma(month(t(i)),1+hour(t(i))),dT,TI,I2,I3,I4,ones(size(domainDTM))*Ti(month(t(i))),ShF,0.75);
            Ta{I}(i,1)=mean2(Tsim);
        end
    end
    progressbar(I/30) % Update progress bar
end
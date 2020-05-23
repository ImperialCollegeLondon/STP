function [ rmse ] = AODcalibration( Tas , beta , AODcorr , lwp , mmm , domainDTM , t , N , Sun , Aspect , Slo_top , HZr , Z , SvF , Ct , Radiation , Vapor , monthlyGlobalRadZ )
domainSize=size(domainDTM,1);
Rsw_tm1=zeros(domainSize,domainSize);
Rsw_tm2=zeros(domainSize,domainSize);
monthly_run=zeros(domainSize,domainSize);
M=0; C=1;
for i=1:8760 % Hourly time step
    m=month(t(i));
    if m==mmm
        Ta=squeeze(Tas(i,:,:));
        esat_s=611.*exp(17.27.*Ta./(237.3+Ta));
        if Sun.altitude(i)>0 || i==1
            [ cos_ST ] = SolarIlluminationAngle( Slo_top , Sun.altitude(i) , Sun.azimuth(i) , Aspect );
            cos_ST(cos_ST<0)=0;
            [ S ] = Shadow_Effect2( domainDTM , pi/2-Sun.altitude(i) , Sun.azimuth(i) , HZr , Z );
        end
        %% Vapor pressure Generator
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
        if C<3
            [ ea , dDe(i) ] = ComputeVapPressure2( esat_s , Ta , 0 , 0 , 0 , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
            C=C+1;
        else
            [ ea , dDe(i) ] = ComputeVapPressure2( esat_s , Ta , Rsw_tm1 , Rsw_tm2 , dDe(i-1) , a0 , a1 , a2 , a3 , dDem , rhodDe , sigmadDe );
        end
        U=ea./esat_s;
        G_g=17.27.*Ta./(237.7+Ta)+log(U);
        Tdew=237.7.*G_g./(17.27-G_g);
        Tdew(isinf(Tdew))=NaN;
        %% Radiation Generator
        if Sun.altitude(i)<=0 && i>1
            Rsw=zeros(domainSize,domainSize);
        else
            [SB,SD]=ComputeRadiationForcings( Sun.altitude(i) , Sun.E0(i) , Radiation.Ro , lwp , N(i) , domainDTM , Tdew , beta , Radiation.Angstrom(m).alpha , Radiation.omega_A1 ,Radiation.omega_A2 , Radiation.uo , Radiation.un , Radiation.rho_g , cos_ST , S , SvF , Ct , AODcorr , Radiation.AERONETz );
            SB(SB<0)=0;
            SD(SD<0)=0;
            Rsw=SB+SD;
        end
        Rsw_tm2=Rsw_tm1;
        Rsw_tm1=Rsw;
        %% Save data
        monthly_run=monthly_run+Rsw;
    end
end
tmp=monthly_run./(eomday(2001,mmm)*24);
c=1;
for d=500:100:3100
    x(c,1)=d+50;
    y(c,1)=max(tmp(domainDTM>=d&domainDTM<d+100));
    c=c+1;
end
a=0;
for i=1:length(y)
    if ~isnan(y(i))
        a=a+(Radiation.monthlyAltitude{mmm}(x(i))-y(i)).^2;
    end
end
rmse=double(sqrt(a./sum(~isnan(y))));
%% Plot
% figure,
% plot(monthlyGlobalRadZ,Radiation.monthlyAltitude{mmm}(monthlyGlobalRadZ),'.');
% hold on
% plot(x,y,'r.')
% title(['m=',num2str(mmm)])
end
function [ b , dTbar , dTrho , dTsigma , Ti , R2 ] = temperature_parameter( Ta , N , d , DeltaGMT , Lon , Lat )
% temperature_parameter A function that evaluate the deterministics and stochastic parameters
% for the temperature simulations on monthly basis
%   Inputs:
%   Ta - temperature [°C]
%	N - cloud fraction [-]
%   d - date and time, matlab format
%   DeltaGMT [h]
%	Lon [°]
%	Lat [°]
%	Outputs:
%	b - regression coefficients of deterministic component [b0...b4]
%	dTbar - mean random temperature deviate [°C]
%	dTrho - lag-1 autocorrelation random temperature deviate
%	dTsigma - standard deviation random temperature deviate
%	Ti - first step temperature [°C]
%   R2 - the r-square of the linear regression
%% Initilazing
Ta(isnan(Ta))=nanmean(Ta);
N(isnan(N))=nanmean(N);
n=length(Ta);
Ta=reshape(Ta,1,n);
N=reshape(N,1,n);
%% Coefficients
Y=diff(Ta); % Hourly temperature change [°C]
t=hour(d);
K=1-0.75*N.^(3.4);
Kn=(1+0.17*N.^2);
qt=(5.6704e-8.*Kn.*(Ta + 273.15).^4)/1000; % Estimate incoming long wave radiation [W/m^2]
jDay=d-datenum(year(d),0,0); % Julian day
delta_S=23.45*pi/180*cos(2*pi/365*(172-jDay)); % Declination of the sun [rad]
T_sunrise=180/(15*pi)*(2*pi-acos(-tan(delta_S)*tan(deg2rad(Lat))))-12;
T_sunset=180/(15*pi)*acos(-tan(delta_S)*tan(deg2rad(Lat)))+12;
if Lon<0 % Difference between standard and local meridian
    Delta_TSL = -1/15*(15*abs(DeltaGMT) - abs(Lon));
else
    Delta_TSL = 1/15*(15*abs(DeltaGMT) - abs(Lon));
end
Is=zeros(1,n-1); Ir=zeros(1,n-1); X1=zeros(1,n-1); X2=zeros(1,n-1); X3=zeros(1,n-1); X4=zeros(1,n-1);
Tas=zeros(1,length(t)); T_tilde=zeros(1,length(t));
%% Computation coefficient
for i=2:n
    [Is(i-1),Ir(i-1)]= Computation_int_rs2(t(i),delta_S(i),deg2rad(Lat),T_sunrise(i),T_sunset(i),Delta_TSL);
    X1(i-1)=Ta(i-1);
    X2(i-1)=((K(i)+K(i-1))/2)*Is(i-1);
    X3(i-1)=((K(i)+K(i-1))/2)*Ir(i-1);
    X4(i-1)= qt(i-1);
end
XX=[ones(1,n-1); X1 ; X2; X3; X4];
[a,~,~,~,STATS]=regress(Y',XX');
R2=STATS(1);
b(2)=-log(1+a(2));
b(1)=(b(2)/(-a(2)))*a(1);
b(3:5)=(b(2)/(-a(2)))*a(3:5);
Tave=mean(Ta);
b=num2cell(b);
for i=1:length(t)
    if i==1
        Ti=Ta(1);
        [Tas(i),T_tilde(i),~,~,I2,I3,I4]=ComputeAirTemperature(d(i),DeltaGMT,Lon,Lat,N(i),b,0,0,0,0,0,Ti,0,0,0,Tave,0,1);
    else
        if t(i)==0
            Ti=T_tilde(i-1);
            I2=0; I3=0; I4=0;
        end
        [Tas(i),T_tilde(i),~,~,I2,I3,I4]=ComputeAirTemperature(d(i),DeltaGMT,Lon,Lat,N(i),b,qt(i-1)*1000,0,0,0,0,Ti,I2,I3,I4,Tave,0,1);
    end
end
clear I2 I3 I4 qtS T_tilde dT
dT = Ta- Tas;  dTbar=zeros(1,24); dTsigma=zeros(1,24);
for i=0:23
    dTh=dT(t==i);
    dTh=dTh(not(isnan(dTh)));
    dTbar(i+1)=mean(dTh); %mean random temperature deviate
    dTsigma(i+1) =std(dTh);%% standard deviation random temperature deviate
    clear  dTh
end
dT=dT(not(isnan(dT)));
R=xcov(dT,dT,10,'coeff');
dTrho=R(12); %lag-1 autocorrelation random temperature deviate
clear R
dTrho(dTrho>0.96)=0.96;
Ti = nanmean(Ta);
%% Nested function 
    function[Is,Ir]=Computation_int_rs2(t1,delta_S,phi,T_sunrise,T_sunset,DeltaTSL)
        T0 = 0 - DeltaTSL;
        T23 = 23.0 - DeltaTSL;
        TL  = t1 - DeltaTSL;
        aa = 0.0005;
        if (DeltaTSL >= 0.0)
            Rho = floor(T_sunrise+1.0) - DeltaTSL;
            if (Rho < T_sunrise)
                Rho=Rho+1;
            end
            Sigma = floor(T_sunset+1.0) - DeltaTSL;
            if (Sigma < T_sunset)
                Sigma=Sigma+1;
            end
            T12 = 13 - DeltaTSL;
        elseif (DeltaTSL < 0.0)
            Rho = floor(T_sunrise) -  DeltaTSL;
            if (Rho < T_sunrise)
                Rho=Rho+1;
            end
            Sigma = floor(T_sunset) -  DeltaTSL;
            if (Sigma < T_sunset)
                Sigma=Sigma+1;
            end
            T12 = 12.0 - DeltaTSL;
        end
        if (T0 <= TL && TL < T_sunrise)
            x2 = 0.0;
            x3 = 0.0;
        elseif (Rho-aa <= TL && Rho+aa >= TL)
            A = pi*T_sunrise/12.0;
            B = pi*Rho/12.0;
            x2 = (Rho-T_sunrise)*sin(phi)*sin(delta_S);
            x2 = x2- ((12.0/pi)*cos(delta_S)*cos(phi)*(sin(B)-sin(A)));
            x3 = cos(delta_S)*cos(phi)*(cos(A)-cos(B));
        elseif (Rho+aa <= TL && TL <= 12)
            A = pi*TL/12.0;
            B = pi*(TL-1)/12.0;
            x2 = sin(phi)*sin(delta_S);
            x2 = x2-((12.0/pi)*cos(delta_S)*cos(phi)*(sin(A)-sin(B)));
            x3 = cos(delta_S)*cos(phi)*(cos(B)-cos(A));
        elseif (T12-aa <= TL && T12+aa >= TL)
            A = pi*TL/12.0;
            B = pi*(TL-1)/12.0;
            C = pi*(T12-1)/12.0;
            x2 = sin(phi)*sin(delta_S);
            x2 = x2- ((12.0/pi)*cos(delta_S)*cos(phi)*(sin(A)-sin(B)));
            x3 = cos(delta_S)*cos(phi)*(cos(C)+1);
        elseif (T12+aa <= TL && TL < T_sunset)
            A = pi*TL/12.0;
            B = pi*(TL-1)/12.0;
            x2 = sin(phi)*sin(delta_S);
            x2 = x2- ((12.0/pi)*cos(delta_S)*cos(phi)*(sin(A)-sin(B)));
            x3 = 0.0;
        elseif (Sigma-aa <= TL && Sigma+aa >= TL)
            A = pi*T_sunset/12.0;
            B = pi*(Sigma-1)/12.0;
            x2 = (T_sunset-Sigma+1.0)*sin(phi)*sin(delta_S);
            x2 = x2+ ((12.0/pi)*cos(delta_S)*cos(phi)*(sin(B)-sin(A)));
            x3 = 0.0;
        elseif (Sigma+aa <= TL && TL <= T23)
            x2 = 0.0;
            x3 = 0.0;
        end
        Is=x2;
        Ir=x3;
    end
end
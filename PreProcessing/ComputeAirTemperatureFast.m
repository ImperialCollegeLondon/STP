function [ Ta , T_tilde , dT , qt , I2 , I3 , I4 ] = ComputeAirTemperatureFast( t , DeltaGMT , Lon , Lat , N , b , qtm1 , dTbar , dTrho ,dTsigma , dTtm1 , T_tildetm1 , I2tm1 , I3tm1 , I4tm1 , Tave )
% ComputeAirTemperatueFast Compute deterministic and stochasic components without shading effect, see details in AWE-GEN techniqal reference
%% Initilazing
jDay=t-datenum(year(t),0,0); % Julian day
delta_S=23.45*pi/180*cos(2*pi/365*(172-jDay)); % Declination of the sun [rad]
T_sunrise=180/(15*pi)*(2*pi-acos(-tan(delta_S)*tan(deg2rad(Lat))))-12;
T_sunset=180/(15*pi)*acos(-tan(delta_S)*tan(deg2rad(Lat)))+12;
if Lon<0 % Difference between standard and local meridian
    Delta_TSL = -1/15*(15*abs(DeltaGMT) - abs(Lon));
else
    Delta_TSL = 1/15*(15*abs(DeltaGMT) - abs(Lon));
end
Kn=1+0.17*N^2;
%% Compute deterministic component
[T_tilde,I2,I3,I4] = ComputeDeterministicT(T_tildetm1,delta_S,deg2rad(Lat),hour(t),T_sunrise,T_sunset,N,b,I2tm1,I3tm1,I4tm1,qtm1,Delta_TSL);
if b{2}<0.001
    T_tilde=Tave;
    I2=0;
    I3=0;
    I4=0;
end
%% Compute random deviate of air temperature
Ta=Inf(size(Tave)); kjx=0;
while all(all((Ta-Tave)>25)) && (kjx<100)
    dT = ComputeStochasticT(dTbar,dTrho,dTtm1,dTsigma);
    Ta = T_tilde + dT;
    kjx=kjx+1;
end
qt=5.6704e-8*Kn*(Ta+273.15)^4;
%% Nested functions
%% Compute stochastic component
    function dT = ComputeStochasticT( dTbar , rhodT , dTtm1 , sigmadT )
        epsT = normrnd(0,1);
        dT = dTbar + rhodT.*(dTtm1 - dTbar) + epsT.*sigmadT.*sqrt(1 - rhodT.^2);
    end
%% Compute deterministic component
    function [ T_tilde , I2 , I3 , I4 ] = ComputeDeterministicT( Ti , delta_S , phi , t , T_sunrise , T_sunset , N , b , I2tm1 , I3tm1 , I4tm1 , qtm1 , DeltaTSL )
        K =  1-0.75*N^3.4;
        TL  = t-DeltaTSL; % Current time step
        TP  = -DeltaTSL-1; % Initial time step
        T12 = 12-DeltaTSL;
        p  = pi/12;
        K1 = b{1}./b{2};
        K2 = (b{3}./b{2}).*sin(delta_S).*sin(phi);
        K3 = ((b{2}.*b{3})./(b{2}.^2 + p.^2)).*cos(delta_S).*cos(phi);
        K4 = ((p.*b{3})./(b{2}.^2 + p.^2)).*cos(delta_S).*cos(phi);
        K5 = ((p.^2.*b{4})./(b{2}.^2 + p.^2)).*cos(delta_S).*cos(phi);
        K6 = ((p.*b{2}.*b{4})./(b{2}.^2 + p.^2)).*cos(delta_S).*cos(phi);
        I1 = K1.*(exp(b{2}.*TL) - exp(b{2}.*TP)); %%% Integral first expression
        I4 = (b{5}./b{2}).*(qtm1./1000).*(1 - exp(-b{2})).*exp(b{2}.*TL) + I4tm1; % Integral 4th expression
        % Integral 2nd expression
        if (TL >= T_sunrise) && (TL <= T_sunset)
            [I2] = integral2(TL,TL-1,b{2},K,K2,K3,K4,I2tm1);
            if TL-1 <= T_sunrise
                [I2] = integral2(TL,T_sunrise,b{2},K,K2,K3,K4,I2tm1);
            end
            if TL+1 >= T_sunset
                [I2] = integral2(T_sunset,TL-1,b{2},K,K2,K3,K4,I2tm1);
            end
            if (TL-1 <= T_sunrise) && (TL+1 >= T_sunset)
                [I2] = integral2(T_sunset,T_sunrise,b{2},K,K2,K3,K4,I2tm1);
            end
        else
            I2=I2tm1;
        end
        % Integral 3rd expression
        if (TL >= T_sunrise) && (TL <= T12)
            [I3] = integral3(TL,TL-1,b{2},K,K5,K6,I3tm1);
            if TL-1 <= T_sunrise
                [I3] = integral3(TL,T_sunrise,b{2},K,K5,K6,I3tm1);
            end
            if TL+1 >= T12
                [I3] = integral3(T12,TL-1,b{2},K,K5,K6,I3tm1);
            end
            if (TL-1 <= T_sunrise) && (TL+1 >= T12)
                [I3] = integral3(T12,T_sunrise,b{2},K,K5,K6,I3tm1);
            end
        else
            I3=I3tm1;
        end
        G = I1 + I2 + I3 + I4;  % Auxiliary variable for the solution
        T_tilde = Ti.*exp(-b{2}.*(TL - TP)) + G.*exp(-b{2}.*TL); % Deterministic component of temperature [°C]
    end
%% Term I2
    function I2 = integral2( t , t2 , b1 , K , K2 , K3 , K4 , I2tm1 )
        I2 = K.*(K2.*(exp(b1.*t)-exp(b1.*(t2)))-K3.*exp(b1.*t).*cos(pi.*t./12)-K4.*exp(b1.*t).*sin(pi.*t./12)+K3.*exp(b1.*(t2)).*cos(pi.*(t2)./12)+K4.*exp(b1.*(t2)).*sin(pi.*(t2)./12))+I2tm1;
    end
%% Term I3
    function I3 = integral3( t , t2 , b1 , K , K5 , K6 , I3tm1 )
        I3 = K.*(K6.*exp(b1.*t).*sin(pi.*t./12)-K5.*exp(b1.*t).*cos(pi.*t./12)-K6.*exp(b1.*(t2)).*sin(pi.*(t2)./12)+K5.*exp(b1.*(t2)).*cos(pi.*(t2)./12)) + I3tm1;
    end
end
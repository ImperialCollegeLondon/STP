function [ Ta , T_tilde , dT , qt , I2 , I3 , I4 ] = ComputeAirTemperature( t , DeltaGMT , Lon , Lat , N , b , qtm1 , dTbar , dTrho ,dTsigma , dTtm1 , T_tildetm1 , I2tm1 , I3tm1 , I4tm1 , Tave , ShF , ShFPar )
% ComputeAirTemperatue Compute deterministic and stochasic components, see details in AWE-GEN techniqal reference
%   Inputs:
%   t - matlab time (year,month,day,hour)
%   DeltaGMT [h]
%   Lon [°]
%   Lat [°]
%   N - cloudiness / shading effect [-]
%   b - regression coefficients of deterministic component [b0 , ... , b4]
%   qtm1 - estimate incoming long wave radiation [W m^-2] from previous time step
%	dTbar - mean random temperature deviate
%	dTrho - lag-1 autocorrelation random temperature deviate
%	dTsigma -  standard deviation random temperature deviate
%	dTtm1 - mean random temperature deviate at pervious time step
%	T_tildetm1 - deterministic component of temperature [°C] at pervious time step
%	I2tm1 , I3tm1 and I4tm1 - integral components at pervious time step
%	Tave - deterministic component of temperature [°C] for the month average
%   ShF - Shading effect, 0 - pixel is shaded and 1 - pixel is not shaded
%   ShFPar - shading effect parameter, [%] of long wave radiation blocakge
%	Outputs:
%	Ta - air temperature [°C]
%	T_tilde - deterministic component of temperature [°C]
%	dT - mean random temperature deviate [°C]
%	qt - estimate incoming long wave radiation [W m^-2]
%	I2,I3,I4 - integral components
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
Kn=1+0.17.*N.^2;
%% Compute deterministic component
% In the following function N is given by its first value only to speed up
% the computation (assuming equal cloud distribution)
[T_tilde,I2,I3,I4] = ComputeDeterministicT(T_tildetm1,delta_S,deg2rad(Lat),hour(t),T_sunrise,T_sunset,N(1),b,I2tm1,I3tm1,I4tm1,qtm1,Delta_TSL,ShF,ShFPar);
T_tilde(b{2}<0.001)=Tave(b{2}<0.001);
I2(b{2}<0.001)=0;
I3(b{2}<0.001)=0;
I4(b{2}<0.001)=0;
%% Compute random deviate of air temperature
Ta=Inf(size(Tave)); kjx=0;
while all(all((Ta-Tave)>25)) && (kjx<100)
    dT = ComputeStochasticT(dTbar,dTrho,dTtm1,dTsigma);
    Ta = T_tilde + dT;
    kjx=kjx+1;
end
% qt=5.6704e-8.*Kn.*(Ta+273.15).^4;
% To speed up the calculation the eqaution as changed, considering some
% minot calculations differences (10^-6 magnitude)
qt=5.6704e-8.*Kn.*(Ta+273.15).*(Ta+273.15).*(Ta+273.15).*(Ta+273.15);
%% Nested functions
%% Compute stochastic component
    function dT = ComputeStochasticT( dTbar , rhodT , dTtm1 , sigmadT )
        [X,Y]=size(dTbar);
        epsT = normrnd(0,1,X,Y);
        dT = dTbar + rhodT.*(dTtm1 - dTbar) + epsT.*sigmadT.*sqrt(1 - rhodT.^2);
    end
%% Compute deterministic component
    function [ T_tilde , I2 , I3 , I4 ] = ComputeDeterministicT( Ti , delta_S , phi , t , T_sunrise , T_sunset , N , b , I2tm1 , I3tm1 , I4tm1 , qtm1 , DeltaTSL , shf , ShFPar )
%         K =  1-0.75.*N.^(3.4); % If cloud are not equally distributed
        K =  1-0.75*N(1)^3.4*ones(size(shf));
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
        if ShFPar<1
            I2(shf==0)=I2tm1(shf==0)+ShFPar.*(I2(shf==0)-I2tm1(shf==0)); % Taking into account the shading effect for the 2nd and 3rd component (direct sun)
            I3(shf==0)=I3tm1(shf==0)+ShFPar.*(I3(shf==0)-I3tm1(shf==0));
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
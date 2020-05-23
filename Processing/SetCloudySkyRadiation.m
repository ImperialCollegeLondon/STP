%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Subfunction SetCloudySkyRadiation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EB1,EB2,ED1,ED2] = SetCloudySkyRadiation(LWP,h_S,Eb1,Eb2,Edp1,Edp2,N,rho_s1,rho_s2,rho_g)
%%%% INPUT
%%h_S,[rad]  solar altitude
%LWP % Liquid Water Path whit cloudiness  [g/m^2]
%%%Eb1 [W/m^2]  clear sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%Eb2  [W/m^2]  clear sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%Edp1 [W/m^2]  clear sky incidente irradiance  diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%Edp2   [W/m^2] clear sky incidente irradiance diffuse flux at the ground NIR band [0.70 um -4.0um ]
%%% N [0-1] cloudiness
%%% rho_s1 rho_s2   clear sky albedo 
%%% rho_g ground albedo 
%%% OUTPUT
%%EB1 [W/m^2]  cloud sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%EB2  [W/m^2]  cloud sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%ED1 [W/m^2]   cloud sky total diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%ED2   [W/m^2]  cloud sky  total diffuse flux at the ground NIR band [0.70 um - 4.0um ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if h_S>0 % if there is sunshine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Partition in Two principals band 
    K = [0.4651 0.5195 0.5195 0.5195]; 
    %% 4 band approach of Slingo 1989
    %% Slingo (1989) considered four spectral  bands, one in UV/VIS and three in NIR wavelength
    %% intervals: [0.25um ÷ 0.69um], [0.69um ÷ 1.19um], [1.19um ÷ 2.38um], [2.38?m ÷ 4.0um]
    %% Parameterization of Slingo (1989)
    k = [0.460 0.326 0.181 0.033]; %% respective fraction of solar irradiance at the top of the atmosphere
    a = [2.817 2.682 2.264 1.281]*1e-2;  %%[m^2/g]
    b = [1.305 1.346 1.454 1.641]; %% [um m^2 /g]
    c = [-5.62e-8 -6.94e-6 4.64e-4 2.01e-1]; %%[]
    d = [1.63e-7 2.35e-5 1.24e-3 7.56e-3]; %%[1/um]
    e = [0.829 0.794 0.754 0.826];  %[]
    f = [2.482 4.226 6.560 4.353]*1e-3; %%[1/um]
    %%% Estimate effective radius of drop-size distribution
    tn1=(10^(0.2633 + 1.7095*log(log10(LWP)))); %% optical thickness band [0.29 um - 0.70 um ]
    tn2=(10^(0.3492 + 1.6518*log(log10(LWP)))); %% optical thickness band [0.70 um - 4.0um ]
    re1 = 1.5*LWP/tn1;  %% [um] effective radius of cloud-droplet size distribution  band [0.29 um - 0.70 um ]
    re2 = 1.5*LWP/tn2;  %% [um] effective radius of cloud-droplet size distribution [0.70 um - 4.0um ]
    re  = [re2 re1 re1 re1]; %%[um] effective radius of cloud-droplet size distribution 4 bands
    re(re < 4.2) = 4.2;  %% range re adjustament
    re(re > 16.6) = 16.6;  %% range re adjustament
    %%%% Compute cloudy sky direct beam irradiance and cloudy sky diffuse beam irradiance
    tau         = LWP.*(a + b./re);  %%  cloud optical depth 4 bands
    omega_tilde = 1 - (c + d.*re);  %% single  scatter  albedo
    g           = e + f.*re; %%  asymmetry parameter
    %%%%%%
    beta0   = (3/7).*(1 - g); %% fraction of the scattered diffuse radiation
    betah   = 0.5 - (3.*sin(h_S).*g)./(4.*(1 + g));  %% fraction of the scattered direct radiation
    f2      = g.^2;
    U1      = 7/4; %% reciprocals of the effective cosines for the diffuse upward flux
    U2      = (7/4).*(1 - (1 - omega_tilde)./(7.*omega_tilde.*beta0));  % reciprocals of the effective cosines for the diffuse  downward flux
    alpha1  = U1.*(1 - omega_tilde.*(1 - beta0));
    alpha2  = U2.*omega_tilde.*beta0;
    alpha3  = (1 - f2).*omega_tilde.*betah;
    alpha4  = (1 - f2).*omega_tilde.*(1 - betah);
    eps     = sqrt(alpha1.^2 - alpha2.^2);
    M       = alpha2./(alpha1 + eps);
    E       = exp(-eps.*tau);
    gamma1  = ((1 - omega_tilde.*f2).*alpha3 - sin(h_S).*(alpha1.*alpha3 + alpha2.*alpha4))./((1 - omega_tilde.*f2).^2 - eps.^2.*sin(h_S).^2);
    gamma2  = (-(1 - omega_tilde.*f2).*alpha4 - sin(h_S).*(alpha1.*alpha4 + alpha2.*alpha3))./((1 - omega_tilde.*f2).^2 - eps.^2.*sin(h_S).^2);
    %%%%%%%%%
    TDB         = exp(-(1 - omega_tilde.*f2).*tau./sin(h_S)); %% cloud transmissivity for the direct beam flux in 4 bands j
    TDB(TDB < 0.0) = 0.0;
    RDIF    = M.*(1 - E.^2)./(1 - E.^2.*M.^2);  %%% diffuse reflectivity for diffuse incident radiation
    RDIF(RDIF < 0.0) = 0.0;
    TDIF    = E.*(1 - M.^2)./(1 - E.^2.*M.^2);  %%%  diffuse transmissivity for diffuse incident radiation
    TDIF(TDIF < 0.0) = 0.0;
    RDIR    = (-gamma2.*RDIF - gamma1.*TDB.*TDIF + gamma1); %%% diffuse reflectivity for direct incident radiation
    TDIR    = (-gamma2.*TDIF - gamma1.*TDB.*RDIF + gamma2.*TDB);  %%% diffuse transmissivity for direct incident radiation
    %%%%%%%%%%%%
    for j=1:4  %%% check TDB + TDIR <= 1
        if  (TDB(j) + TDIR(j)) > 1
            SDB_DIR = TDB(j) + TDIR(j);
            TDIR(j) = TDIR(j)/SDB_DIR;
            TDB(j)  = TDB(j)/SDB_DIR;
        end
    end
    %%%%%%% cloud sky beam irradiance
    EB1     = Eb1.*((1 - N) + TDB(1).*N).*k(1)./K(1);
    EB2_multi=(TDB(2:4).*k(2:4))./K(2:4);
    EB2_all = Eb2.*EB2_multi(1)+Eb2.*EB2_multi(2)+Eb2.*EB2_multi(3);
    EB2     = (1 - N).*Eb2 + N.*EB2_all;
    %%%%%% cloud sky  total diffuse flux
    EDp1     = (1-N).*Edp1 + N.*(TDIR(1).*Eb1 + TDIF(1).*Edp1).*k(1)./K(1);  %%
    EDp2_multi=k(2:4)./K(2:4);
    EDp2_a1 =TDIR(2).*Eb2+TDIR(3).*Eb2+TDIR(4).*Eb2;
    EDp2_a2 =TDIF(2).*Edp2+TDIF(3).*Edp2+TDIF(4).*Edp2;
    EDp2_all = (EDp2_a1+EDp2_a2).*EDp2_multi(1)+(EDp2_a1+EDp2_a2).*EDp2_multi(2)+(EDp2_a1+EDp2_a2).*EDp2_multi(3);
    EDp2     = (1 - N).*Edp2 + N.*EDp2_all;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sky albedo for cloudy sky 
    rho_csb1 = (1-N)*rho_s1 + N*RDIR(1)*k(1)/K(1); 
    rho_csd1 = (1-N)*rho_s1 + N*RDIF(1)*k(1)/K(1); 
    rho_csb2 = (1-N)*rho_s2 + N*sum(RDIR(2:4).*k(2:4)./K(2:4)); 
    rho_csd2 = (1-N)*rho_s2 + N*sum(RDIF(2:4).*k(2:4)./K(2:4)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Backscattered diffuse component  after sky attenuation 
    EDd1= rho_g.*rho_csb1.*(EB1.*sin(h_S))./(1-rho_g.*rho_csb1) +  rho_g.*rho_csd1.*(EDp1)./(1-rho_g.*rho_csd1);
    EDd2= rho_g.*rho_csb2.*(EB2.*sin(h_S))./(1-rho_g.*rho_csb2) +  rho_g.*rho_csd2.*(EDp2)./(1-rho_g.*rho_csd2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Diffuse radiation 
    ED1= EDp1 +EDd1; 
    ED2= EDp2 +EDd2; 
else % not sunshine
    EB1 = zeros(size(Eb1));
    EB2 = zeros(size(Eb1));
    ED1 = zeros(size(Eb1));
    ED2 = zeros(size(Eb1));
end
end
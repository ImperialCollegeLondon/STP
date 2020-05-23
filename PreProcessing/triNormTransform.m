function [ Norm , Copula , MaternCov ] = triNormTransform( IsEvent , rIMF , rWAR , mCAR , rTime , prdata , clddata )
%trinormTransform A function that: (1) normalize the WAR, CAR and IMF; (2) Copula
%between WAR, CAR IMF and event duration; (3) fits the VARFIMA parameters.
%   Inputs:
%   IsEvent - ID for each wet event (0 is dry interval)
%   rIMF - mean areal rain intensity [mm h^-1] per each time step
%   rWAR - fraction of the wet area ratio [-]
%   mCAR - fraction of the cloud cover [-]
%   rTime - matlab time series
%   prdata - IMF and WAR parameters
%   clddata - CAR parameters
%   Outputs:
%   Norm - structure that saves the normalizes coefficient and also WARg
%   and IMFg
%   Copula - coefficients for the T-copula
%   MartenCov - store the coefficients for the Marten covariance function
%   for the VARFIMA
%% Initilazing
c(1:12)=1;
for i=1:max(IsEvent)
    if sum(IsEvent==i)>=6
        IMF=rIMF(IsEvent==i);
        imf=smooth(IMF);
        WAR=rWAR(IsEvent==i);
        war=smooth(WAR);
        CAR=mCAR(IsEvent==i);
        car=smooth(CAR);
        T=rTime(IsEvent==i);
        m=month(T(1));
        %% Normalization transform
        rank_war = Rankings(war,'competition2');
        rank_car = Rankings(car,'competition2');
        rank_imf = Rankings(imf,'competition2');
        Norm(m).war(c(m)).data = norminv(rank_war/(length(rank_war)+1),0,1);
        Norm(m).car(c(m)).data = norminv(rank_car/(length(rank_car)+1),0,1);
        Norm(m).imf(c(m)).data = norminv(rank_imf/(length(rank_imf)+1),0,1);
        Norm(m).duration(c(m)).data = sum(IsEvent==i);
        c(m)=c(m)+1;
    end
end
%% Copula
for m=1:12
    for i=1:size(Norm(m).duration,2)
        Norm(m).war(i).miu=mean(Norm(m).war(i).data);
        Norm(m).war(i).std=std(Norm(m).war(i).data);
        Norm(m).car(i).miu=mean(Norm(m).car(i).data);
        Norm(m).car(i).std=std(Norm(m).car(i).data);
        Norm(m).imf(i).miu=mean(Norm(m).imf(i).data);
        Norm(m).imf(i).std=std(Norm(m).imf(i).data);
    end
    u_q = ksdensity([Norm(m).duration.data],[Norm(m).duration.data],'function','cdf');
    v_q = ksdensity([Norm(m).war(:).miu],[Norm(m).war(:).miu],'function','cdf');
    w_q = ksdensity([Norm(m).war(:).std],[Norm(m).war(:).std],'function','cdf');
    z_q = ksdensity([Norm(m).imf(:).miu],[Norm(m).imf(:).miu],'function','cdf');
    q_q = ksdensity([Norm(m).imf(:).std],[Norm(m).imf(:).std],'function','cdf');
    vv_q = ksdensity([Norm(m).car(:).miu],[Norm(m).car(:).miu],'function','cdf');
    ww_q = ksdensity([Norm(m).car(:).std],[Norm(m).car(:).std],'function','cdf');
    [Rho,nu] = copulafit('t',[u_q' v_q' w_q' z_q' q_q' vv_q' ww_q'],'Method','ApproximateML');
    [a_swar, ~] = gamfit([Norm(m).war(:).std]);
    [a_scar, ~] = gamfit([Norm(m).car(:).std]);
    [a_simf, ~] = gamfit([Norm(m).imf(:).std]);
    Copula(m).MWAR1 = mean([Norm(m).war(:).miu]);
    Copula(m).MWAR2 = std([Norm(m).war(:).miu]);
    Copula(m).SWAR1 = a_swar(1);
    Copula(m).SWAR2 = a_swar(2);
    Copula(m).MCAR1 = mean([Norm(m).car(:).miu]);
    Copula(m).MCAR2 = std([Norm(m).car(:).miu]);
    Copula(m).SCAR1 = a_scar(1);
    Copula(m).SCAR2 = a_scar(2);
    Copula(m).MIMF1 = mean([Norm(m).imf(:).miu]);
    Copula(m).MIMF2 = std([Norm(m).imf(:).miu]);
    Copula(m).SIMF1 = a_simf(1);
    Copula(m).SIMF2 = a_simf(2);
    Copula(m).cop_rho = Rho;
    Copula(m).cop_nu = nu;
end
%% ARFIMA(q) calculation
for m=1:12
    [ N1 , N2 , N3 , A , N12 , N13 , N23 , RHO12 , RHO13 , RHO23 ] = TrivariateParsimoniousMaternFit( prdata.WAR(m).WAR , prdata.IMF(m).IMF , clddata.Wet(m).CAR , 100 );
    MaternCov(m).data = [ N1 , N2 , N3 , A , N12 , N13 , N23 , RHO12 , RHO13 , RHO23 ];
end
end
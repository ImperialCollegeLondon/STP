function [ tsWARFIMA ] = MultiTriVARFIMA( MaternCov , prFit , cldFit , Copula , m , eventDuration , K )
%% Generating mean WAR and IMF by copula
r = copularnd('t',Copula(m).cop_rho,Copula(m).cop_nu,1);
mwar = icdf('normal',r(:,2),Copula(m).MWAR1,Copula(m).MWAR2);
swar = icdf('gamma',r(:,3),Copula(m).SWAR1,Copula(m).SWAR2);
mcar = icdf('normal',r(:,6),Copula(m).MCAR1,Copula(m).MCAR2);
scar = icdf('gamma',r(:,7),Copula(m).SCAR1,Copula(m).SCAR2);
mimf = icdf('normal',r(:,4),Copula(m).MIMF1,Copula(m).MIMF2);
simf = icdf('gamma',r(:,5),Copula(m).SIMF1,Copula(m).SIMF2);
%% WAR - IMF dissagregation
coef = MaternCov(m).data;
N1 = coef(1);
N2 = coef(2);
N3 = coef(3);
A = coef(4);
N12 = coef(5);
N13 = coef(6);
N23 = coef(7);
RHO12 = coef(8);
RHO13 = coef(9);
RHO23 = coef(10);
[ WAR_norm , IMF_norm , CAR_norm ] = MultiTrivariateMaternSim( N1 , N2 , N3 , A , N12 , N13 , N23 , RHO12 , RHO13 , RHO23 , max(eventDuration,1024) , K );
IMF_norm=IMF_norm(1:eventDuration,:);
WAR_norm=WAR_norm(1:eventDuration,:);
CAR_norm=CAR_norm(1:eventDuration,:);
IMF_norm = bsxfun(@rdivide,bsxfun(@minus,IMF_norm,nanmean(IMF_norm)),nanstd(IMF_norm))*simf+mimf;
WAR_norm = bsxfun(@rdivide,bsxfun(@minus,WAR_norm,nanmean(WAR_norm)),nanstd(WAR_norm))*swar+mwar;
CAR_norm = bsxfun(@rdivide,bsxfun(@minus,CAR_norm,nanmean(CAR_norm)),nanstd(CAR_norm))*scar+mcar;
temp_imf_n=normcdf(IMF_norm);
temp_war_n=normcdf(WAR_norm);
temp_car_n=normcdf(CAR_norm);
try
    temp_war_n=prFit.WAR(m).fit.icdf(temp_war_n);
catch
    jType=prFit.WAR(m).fit.type; % Johnson dist.
    temp_war_n=multi_f_johnson_inv(temp_war_n,prFit.WAR(m).fit.coef,jType);
end
try
    temp_car_n=cldFit.Wet(m).fit.icdf(temp_car_n);
catch
    jType=cldFit.Wet(m).fit.type; % Johnson dist.
    temp_car_n=multi_f_johnson_inv(temp_car_n,cldFit.Wet(m).fit.coef,jType);
end
switch prFit.IMF(m).case
    case {1,2,3}
        temp_imf_n=prFit.IMF(m).fit.icdf(temp_imf_n);
    case 4
        icdf_expmixture = @(x,alphaX,a,mu1,mu2) a*expcdf(x,mu1) + (1-a)*expcdf(x,mu2)-alphaX;
        for III=1:K*eventDuration
            fun = @(x) icdf_expmixture(x,temp_imf_n(III),prFit.IMF(m).fit(1),prFit.IMF(m).fit(2),prFit.IMF(m).fit(3));
            temp_imf_n(III)=fzero(fun,0.9);
        end
    case 5
        icdfMGGP = @(x,alphaX,W,A,B,K,S) (1-W).*gamcdf(x,A,B) + W.*gpcdf(x,K,S)-alphaX;
        temp_imf_n2=arrayfun(@(i) fzero(@(x) icdfMGGP(x,temp_imf_n(i),prFit.IMF(m).fit(1),prFit.IMF(m).fit(2),prFit.IMF(m).fit(3),prFit.IMF(m).fit(4),prFit.IMF(m).fit(5)),0.9),1:numel(temp_imf_n));
        temp_imf_n=reshape(temp_imf_n2,size(temp_imf_n,1),size(temp_imf_n,2));
end
temp_war_n(temp_war_n<0)=0;
temp_war_n(temp_war_n>1)=1;
temp_car_n(temp_car_n>1)=1;
temp_car_n(isnan(temp_car_n))=1;
tsWARFIMA.IMFnorm = IMF_norm;
tsWARFIMA.IMFreal = temp_imf_n;
tsWARFIMA.WARnorm = WAR_norm;
tsWARFIMA.WARreal = temp_war_n;
tsWARFIMA.CARnorm = CAR_norm;
tsWARFIMA.CARreal = temp_car_n;
end
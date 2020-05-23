function [ Gamma , LogNorm , GP , ME , MGGP , AIC , AICc , BIC ] = fitIMF( IMF , isPlot )
%fitIMF Fit a Gamma, Generalized Pareto (GP), Mixed Exponential (ME),
%LogNormal (LN), hybrid Gamma and Generalized Pareto (FK08) and mixture
%model of Gamma and Generalized Pareto (MGGP) distributions for the 
%IMF. The HGGP (FK08) is not active, value can be improved if theta is subjectively
%predetrmined.
%   Inputs:
%   IMF - mean daily areal rainfall [mm h^-1]
%   Outputs:
%   Gamma , LogNorm , GP , ME , MGGP and FK08 - distribution objects or
%   parameters.
%   AIC , AICc , BIC - criterions parameters
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Fitting
Gamma=fitdist(IMF,'gamma');
LogNorm=fitdist(IMF,'lognormal');
GP=fitdist(IMF,'GeneralizedPareto');
pdf_expmixture = @(IMF,a,mu1,mu2) a*exppdf(IMF,mu1) + (1-a)*exppdf(IMF,mu2); % Mixed exponential
cdf_expmixture = @(x,a,mu1,mu2) a*expcdf(x,mu1) + (1-a)*expcdf(x,mu2);
muStart = quantile(IMF,[.25 .75]);
start = [0.5 muStart];
lb = [0 0 0];
ub = [1 Inf Inf];
ME = mle(IMF,'pdf',pdf_expmixture,'start',start,'lower',lb,'upper',ub);
pdfMGGP = @(IMF,W,A,B,K,S) (1-W).*gampdf(IMF,A,B) + W.*gppdf(IMF,K,S); % Mixed Gamma and GP
cdfMGGP = @(x,W,A,B,K,S) (1-W).*gamcdf(x,A,B) + W.*gpcdf(x,K,S);
MGGPstart = [0.5 1 1 1 1];
MGGPlb = [0 0 0 0 0];
MGGPub = [1 Inf Inf Inf Inf];
MGGPoptions = statset('MaxIter',10000,'MaxFunEvals',15000,'TolX',1e-30,'TolFun',1e-30);
MGGP = mle(IMF,'pdf',pdfMGGP,'cdf',cdfMGGP,'start',MGGPstart,'lower',MGGPlb,'upper',MGGPub, 'options',MGGPoptions);
% pdfFK08 = @(IMF,a,b,k,theta) gampdf(IMF,a,b).*Ibol(IMF,theta)+(1-gamcdf(IMF,a,b)).*gppdf(IMF,k,(1-gamcdf(theta,a,b))./gampdf(theta,a,b),theta).*Ibol2(IMF,theta);
% cdfFK08 = @(x,a,b,k,theta) gamcdf(x,a,b).*Ibol(x,theta)+(1-gamcdf(x,a,b)).*gpcdf(x,k,(1-gamcdf(theta,a,b))./gampdf(theta,a,b),theta).*Ibol2(x,theta)+gamcdf(x,a,b).*Ibol2(x,theta);
% startFK08 = [1 1 1 1];
% lbFK08 = [0 0 0 0];
% ubFK08 = [Inf Inf Inf Inf];
% optionsFK08 = statset('MaxIter',10000,'MaxFunEvals',15000,'TolX',1e-100,'TolFun',1e-100);
% FK08 = mle(IMF,'pdf',pdfFK08,'start',startFK08,'lower',lbFK08,'upper',ubFK08, 'options',optionsFK08);
%% Estimate
% BIC - Bayesian information criterion
% AIC - Akaike information criterion
% AICc - AIC with a correction for finite sample sizes
%Beta:
NLL=Gamma.NLogL; % -Log(L)
k=numel(Gamma.Params); % Number of parameters
n=numel(IMF); % Number of data points
BIC.gamma=-2*(-NLL)+k*log(n);
AIC.gamma=-2*(-NLL)+2*k;
AICc.gamma=AIC.gamma+((2*k*(k+1))/(n-k-1));
%LogNormal:
NLL=LogNorm.NLogL; % -Log(L)
k=numel(LogNorm.Params); % Number of parameters
BIC.LN=-2*(-NLL)+k*log(n);
AIC.LN=-2*(-NLL)+2*k;
AICc.LN=AIC.LN+((2*k*(k+1))/(n-k-1));
%GP:
NLL=GP.NLogL; % -Log(L)
k=numel(GP.Params); % Number of parameters
BIC.GP=-2*(-NLL)+k*log(n);
AIC.GP=-2*(-NLL)+2*k;
AICc.GP=AIC.GP+((2*k*(k+1))/(n-k-1));
%ME:
k=3;
NLL=-sum(log(pdf_expmixture(IMF,ME(1),ME(2),ME(3))));
BIC.ME=-2*(-NLL)+k*log(n);
AIC.ME=-2*(-NLL)+2*k;
AICc.ME=AIC.ME+((2*k*(k+1))/(n-k-1));
%MGGP:
k=5;
NLL=-sum(log(pdfMGGP(IMF,MGGP(1),MGGP(2),MGGP(3),MGGP(4),MGGP(5))));
BIC.MGGP=-2*(-NLL)+k*log(n);
AIC.MGGP=-2*(-NLL)+2*k;
AICc.MGGP=AIC.MGGP+((2*k*(k+1))/(n-k-1));
%HGGP:
% k=4;
% NLL=-sum(log(pdfFK08(IMF,FK08(1),FK08(2),FK08(3),FK08(4))));
% BIC.FK08=-2*(-NLL)+k*log(n);
% AIC.FK08=-2*(-NLL)+2*k;
% AICc.FK08=AIC.FK08+((2*k*(k+1))/(n-k-1));
%% Plot
if isPlot
    x=linspace(min(IMF),max(IMF),100);
    figure('units','normalized','outerposition',[0 0 1 0.5])
    %Gamma
    subplot(1,5,1)
    h1=plot(x,Gamma.cdf(x),'r');
    hold on
    ecdf(IMF);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','southeast')
    xlabel('Precipitation [mm h^{-1}]')
    ylabel('CDF')
    title(['Gamma']);
    axis('square')
    %Mixed Exponential
    subplot(1,5,2)
    h1=plot(x,cdf_expmixture(x,ME(1),ME(2),ME(3)),'r');
    hold on
    ecdf(IMF);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','southeast')
    xlabel('Precipitation [mm h^{-1}]')
    ylabel('CDF')
    title(['Mixed Exponential']);
    axis('square')
    %Mixed Gamma and Generelaized Pareto
    subplot(1,5,3)
    h1=plot(x,cdfMGGP(x,MGGP(1),MGGP(2),MGGP(3),MGGP(4),MGGP(5)),'r');
    hold on
    ecdf(IMF);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','southeast')
    xlabel('Precipitation [mm h^{-1}]')
    ylabel('CDF')
    title(['Mixed Gamma and GP']);
    axis('square')
    %Generelazied Pareto
    subplot(1,5,4)
    h1=plot(x,GP.cdf(x),'r');
    hold on
    ecdf(IMF);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','southeast')
    xlabel('Precipitation [mm h^{-1}]')
    ylabel('CDF')
    title(['Generalized Pareto']);
    axis('square')
    %Log Normal
    subplot(1,5,5)
    h1=plot(x,LogNorm.cdf(x),'r');
    hold on
    ecdf(IMF);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','southeast')
    xlabel('Precipitation [mm h^{-1}]')
    ylabel('CDF')
    title(['LogNormal']);
    axis('square')
    %Hybrid Gamma and Generelaized Pareto
%     subplot(2,3,6)
%     h1=plot(x,cdfFK08(x,FK08(1),FK08(2),FK08(3),FK08(4)),'r');
%     hold on
%     ecdf(IMF);
%     h1 = get(gca,'children');
%     legend('HGGP distribution','Empirical','location','southeast')
%     xlabel('Precipitation [mm h^{-1}]')
%     ylabel('CDF')
%     title(['Hybrid Gamma and Generalized Pareto distribution fit (AICc=',num2str(AICc.FK08),')']);
%     axis('square')
end
end
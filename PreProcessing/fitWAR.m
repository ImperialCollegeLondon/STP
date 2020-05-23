function [ BETAmodel , LogNorm , John , AIC , AICc , BIC ] = fitWAR( WAR , isPlot )
%fitWAR Fit Beta, Lognormal or Johnson distribution for the wet area ratio
%   Inputs:
%   WAR - wet area ratio [fraction]
%   Outputs:
%   BETAmodel , LogNorm and John - Beta or lognormal distribution objects;
%   Johnson parameters
%   AIC , AICc , BIC - criterions parameters
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Fitting
BETAmodel=fitdist(WAR,'beta');
LogNorm=fitdist(WAR,'lognormal');
John=f_johnson_fit(WAR,'M');
%% Estimate
% BIC - Bayesian information criterion
% AIC - Akaike information criterion
% AICc - AIC with a correction for finite sample sizes
%Beta:
NLL=BETAmodel.NLogL; % -Log(L)
k=numel(BETAmodel.Params); % Number of parameters
n=numel(WAR); % Number of data points
BIC.beta=-2*(-NLL)+k*log(n);
AIC.beta=-2*(-NLL)+2*k;
AICc.beta=AIC.beta+((2*k*(k+1))/(n-k-1));
%LogNormal:
NLL=LogNorm.NLogL; % -Log(L)
k=numel(LogNorm.Params); % Number of parameters
n=numel(WAR); % Number of data points
BIC.LN=-2*(-NLL)+k*log(n);
AIC.LN=-2*(-NLL)+2*k;
AICc.LN=AIC.LN+((2*k*(k+1))/(n-k-1));
%John
[AIC.John,AICc.John,BIC.John] = f_johnson_aic(WAR,John.coef,John.type);
%% Plot
if isPlot
    x=linspace(min(WAR),max(WAR),100);
    figure('units','normalized','outerposition',[0 0 1 0.5])
    %Beta
    subplot(1,5,2)
    h1=plot(x,BETAmodel.cdf(x),'r');
    hold on
    ecdf(WAR);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','northwest')
    xlabel('Wet area ratio [-]')
    ylabel('CDF')
    title(['Beta']);
    axis('square')
    %LogNormal
    subplot(1,5,3)
    h1=plot(x,LogNorm.cdf(x),'r');
    hold on
    ecdf(WAR);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','northwest')
    xlabel('Wet area ratio [-]')
    ylabel('CDF')
    title(['LogNormal']);
    axis('square')
    %Johnson
    subplot(1,5,4)
    try
        h1=plot(x,f_johnson_cdf(x',John.coef,John.type),'r');
    catch
        x=linspace(min(WAR)+.1,max(WAR)-.1,100);
        h1=plot(x,f_johnson_cdf(x',John.coef,John.type),'r');
    end
    hold on
    ecdf(WAR);
    h1 = get(gca,'children');
    legend(['Fit [',John.type,']'],'Empirical','location','northwest')
    xlabel('Wet area ratio [-]')
    ylabel('CDF')
    title(['Johnson']);
    axis('square')
end
end
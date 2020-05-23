function [ LogNorm , Normal , Weibull , John , AIC , AICc , BIC ] = fitCAR( CAR , isPlot )
%fitCAR Fit a Log Normal, Normal, Weibull and Johnson distributions for cloud areal ratio. 
%   Inputs:
%   CAR - areal cloud ratio [-]
%   Outputs:
%   LogNorm , Normal , Weibull and John - distribution objects or parameters.
%   AIC , AICc , BIC - criterions parameters
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Fitting
try
    LogNorm=fitdist(CAR,'lognormal');
catch
    CAR(CAR==0)=0.001;
    LogNorm=fitdist(CAR,'lognormal');
end
Normal=fitdist(CAR,'normal');
Weibull=fitdist(CAR,'weibull');
John=f_johnson_fit(CAR,'M');
%% Estimate
% BIC - Bayesian information criterion
% AIC - Akaike information criterion
% AICc - AIC with a correction for finite sample sizes
n=numel(CAR); % Number of data points
%LogNormal:
NLL=LogNorm.NLogL; % -Log(L)
k=numel(LogNorm.Params); % Number of parameters
BIC.LN=-2*(-NLL)+k*log(n);
AIC.LN=-2*(-NLL)+2*k;
AICc.LN=AIC.LN+((2*k*(k+1))/(n-k-1));
%Normal:
NLL=Normal.NLogL; % -Log(L)
k=numel(Normal.Params); % Number of parameters
BIC.Normal=-2*(-NLL)+k*log(n);
AIC.Normal=-2*(-NLL)+2*k;
AICc.Normal=AIC.Normal+((2*k*(k+1))/(n-k-1));
%Weibull:
NLL=Weibull.NLogL; % -Log(L)
k=numel(Weibull.Params); % Number of parameters
BIC.Weibull=-2*(-NLL)+k*log(n);
AIC.Weibull=-2*(-NLL)+2*k;
AICc.Weibull=AIC.Weibull+((2*k*(k+1))/(n-k-1));
%John:
[AIC.John,AICc.John,BIC.John] = f_johnson_aic(CAR,John.coef,John.type);
%% Plot
if isPlot
    x=linspace(min(CAR),max(CAR),100);
    figure('units','normalized','outerposition',[0 0 1 0.5])
    %Weibull
    subplot(1,5,1)
    h1=plot(x,Weibull.cdf(x),'r');
    hold on
    ecdf(CAR);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','northwest')
    xlabel('Areal cloud ratio [-]')
    ylabel('CDF')
    title(['Weibull']);
    axis('square')
    %Johnson
    subplot(1,5,2)
    try
        h1=plot(x,f_johnson_cdf(x',John.coef,John.type),'r');
    catch
        try
            x=linspace(min(CAR)+.1,max(CAR)-.1,100);
            h1=plot(x,f_johnson_cdf(x',John.coef,John.type),'r');
        catch
            x=linspace(min(CAR)+.2,max(CAR)-.2,100);
            h1=plot(x,f_johnson_cdf(x',John.coef,John.type),'r');
        end
    end
    hold on
    ecdf(CAR);
    h1 = get(gca,'children');
    legend(['Fit [',John.type,']'],'Empirical','location','northwest')
    xlabel('Areal cloud ratio [-]')
    ylabel('CDF')
    title(['Johnson']);
    axis('square')
    %Normal
    subplot(1,5,3)
    h1=plot(x,Normal.cdf(x),'r');
    hold on
    ecdf(CAR);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','northwest')
    xlabel('Areal cloud ratio [-]')
    ylabel('CDF')
    title(['Normal']);
    axis('square')
    %Log Normal
    subplot(1,5,4)
    h1=plot(x,LogNorm.cdf(x),'r');
    hold on
    ecdf(CAR);
    h1 = get(gca,'children');
    legend('Fit','Empirical','location','northwest')
    xlabel('Areal cloud ratio [-]')
    ylabel('CDF')
    title(['LogNormal']);
    axis('square')
end
end
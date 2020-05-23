function [ dry , wet , dryFit , wetFit ] = WetDryG( isWetEvent , isWetEventYear , GenParetoTheta , dt )
%eventPR Define time series of wet events by thresholds.
%   Inputs:
%   isWetEvent - boolean time series of the observed data, 1 for wet time step 0 for dry
%   GenParetoTheta - Theta value for the Generalized Pareto distribution
%   (i.e., the lowest statring point for the distribution, 1 - one time step)
%   dt - time step [min]
%   Output:
%   dry and wet - the length of each dry or wet event
%   dryFit and wetFit - a lognormal or generalized pareto distribution fit
%   for the dry and wet events series
%% Initilazing
dry(isWetEvent(1)==0)=1;
wet(isWetEvent(1)==1)=1;
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Wet / dry time series
for i=2:size(isWetEvent,1)
    if isWetEvent(i)==isWetEvent(i-1)
        switch isWetEvent(i)
            case 0
                if isWetEventYear(i)==isWetEventYear(i-1)
                    dry(size(dry,1),1)=dry(size(dry,1),1)+1;
                else
                    dry(size(dry,1)+1,1)=1;
                end
            case 1
                if isWetEventYear(i)==isWetEventYear(i-1)
                    wet(size(wet,1),1)=wet(size(wet,1),1)+1;
                else
                    wet(size(wet,1)+1,1)=1;
                end
        end
    else
        switch isWetEvent(i)
            case 0
                dry(size(dry,1)+1,1)=1;
            case 1
                wet(size(wet,1)+1,1)=1;
        end
    end
end
dry(dry==0)=[];
wet(wet==0)=[];
%% Fitting distributions
try
    wetFit1=fitdist(wet,'GeneralizedPareto','Theta',GenParetoTheta);
catch
    wetFit1=fitdist(wet,'GeneralizedPareto','Theta',0);
end
try
    dryFit1=fitdist(dry,'GeneralizedPareto','Theta',GenParetoTheta);
catch
    dryFit1=fitdist(dry,'GeneralizedPareto','Theta',0);
end
wetFit2=fitdist(wet,'LogNormal');
dryFit2=fitdist(dry,'LogNormal');
wetFit3=f_johnson_fit(wet,'M');
dryFit3=f_johnson_fit(dry,'M');
[fw xw]=ecdf(wet);
[fd xd]=ecdf(dry);
wetFit3.x=linspace(min(wet),max(wet),100);
dryFit3.x=linspace(min(dry),max(dry),100);
LL1=wetFit1.negloglik; LL2=wetFit2.negloglik; LL3=f_johnson_lik(wet,wetFit3.coef,wetFit3.type);
cW=find([LL1 LL2 LL3]==max([LL1 LL2 LL3]));
switch cW
    case 1
        wetFit=wetFit1;
    case 2
        wetFit=wetFit2;
    case 3
        wetFit=wetFit3;
end
LL1=dryFit1.negloglik; LL2=dryFit2.negloglik; LL3=f_johnson_lik(dry,dryFit3.coef,dryFit3.type);
cD=find([LL1 LL2 LL3]==max([LL1 LL2 LL3]));
switch cD
    case 1
        dryFit=dryFit1;
    case 2
        dryFit=dryFit2;
    case 3
        dryFit=dryFit3;
end
%% Plot
figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(2,1,1)
ecdf(wet);
hold on
switch cW
    case {1,2}
        h2=plot(wetFit3.x,wetFit.cdf(wetFit3.x),'r');
        legend(h2,[wetFit.DistributionName],'location','SouthEast')
    case 3
         try
            h2=plot(wetFit.x,f_johnson_cdf(wetFit.x',wetFit.coef,wetFit.type),'r');
        catch
            wetFit.x=linspace(min(wet),max(wet)*.9,100);
            h2=plot(wetFit.x,f_johnson_cdf(wetFit.x',wetFit.coef,wetFit.type),'r');
        end
        legend(h2,['Johnson [',wetFit.type,']'],'location','SouthEast')
end
title('Wet period')
xlabel(['Wet length [',num2str(dt),'-min intervals]'])
axis('square')
subplot(2,1,2)
ecdf(dry);
hold on
switch cD
    case {1,2}
        h2=plot(dryFit3.x,dryFit.cdf(dryFit3.x),'r');
        legend(h2,[dryFit.DistributionName],'location','SouthEast')
    case 3
        try
            h2=plot(dryFit.x,f_johnson_cdf(dryFit.x',dryFit.coef,dryFit.type),'r');
        catch
            dryFit.x=linspace(min(dry)+1,max(dry)*.9,100);
            h2=plot(dryFit.x,f_johnson_cdf(dryFit.x',dryFit.coef,dryFit.type),'r');
        end
        legend(h2,['Johnson [',dryFit.type,']'],'location','SouthEast')
end
title('Dry period')
xlabel(['Dry length [',num2str(dt),'-min intervals]'])
axis('square')
end
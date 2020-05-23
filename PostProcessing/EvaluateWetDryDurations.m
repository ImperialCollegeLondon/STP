%% Initilazing
rng('shuffle')
dt=10; % Time interval [min]
II=1000; % Number of realizations to generate
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Load
load('gaugeIsEvent.mat')
load('Precipitation.mat', 'data')
%% Number of events and total wet event length - Observed
yS=min(year(gaugeIsEvent(:,2)));
yE=max(year(gaugeIsEvent(:,2)));
for i=yS:yE
    for m=1:12
        wetSumObs(i-yS+1,m)=sum(gaugeIsEvent(month(gaugeIsEvent(:,2))==m&year(gaugeIsEvent(:,2))==i,1)>0)*dt;
        tmp=gaugeIsEvent(month(gaugeIsEvent(:,2))==m&year(gaugeIsEvent(:,2))==i,1);
        tmp=unique(tmp);
        eventSumObs(i-yS+1,m)=sum(tmp>0);
    end
end
%% Number of events and total wet event length - Simulation
clear wetSum eventSum wetS dryS
for m=1:12
    try
        dryPool{m}=f_johnson_rnd(data(m).dryFitG.coef,data(m).dryFitG.type,200*II,1);
    catch
        dryPool{m}=data(m).dryFitG.random(200*II,1);
    end
    try
        wetPool{m}=data(m).wetFitG.random(200*II,1);
    catch
        wetPool{m}=f_johnson_rnd(data(m).wetFitG.coef,data(m).wetFitG.type,200*II,1);
    end
    dryPool{m}(dryPool{m}<1)=1;
    wetPool{m}(wetPool{m}<1)=1;
    dryPool{m}=ceil(dryPool{m});
    wetPool{m}=floor(wetPool{m});
end
progressbar('Generating the ensemble') % Init single bar
for I=1:II
    [ General.simulation.tsMonth , General.simulation.tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , dt/(24*60) );
    [ tsMatrix ] = tsTable( General.simulation.tsWetDry , General.simulation.tsMonth );
    tsMatrix(:,5)=tsMatrix(:,4)-tsMatrix(:,3)+1;
    for m=1:12
        eventSum(I,m)=sum(tsMatrix(tsMatrix(:,2)==m,1));
        wetSum(I,m)=sum(tsMatrix(tsMatrix(:,1)==1&tsMatrix(:,2)==m,5))*dt;
        if I==1
            wetS{m}=(tsMatrix(tsMatrix(:,1)==1&tsMatrix(:,2)==m,5))*dt;
        else
            wetS{m}=cat(1,wetS{m},(tsMatrix(tsMatrix(:,1)==1&tsMatrix(:,2)==m,5))*dt);
        end
        if I==1
            dryS{m}=(tsMatrix(tsMatrix(:,1)==0&tsMatrix(:,2)==m,5))*dt;
        else
            dryS{m}=cat(1,dryS{m},(tsMatrix(tsMatrix(:,1)==0&tsMatrix(:,2)==m,5))*dt);
        end
    end
    progressbar(I/II) % Update progress bar
end
%% Plot
figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(2,1,1)
errorbar(1:12,mean(eventSumObs),std(eventSumObs),'.','MarkerSize',15)
hold on
errorbar(1:12,mean(eventSum),std(eventSum),'.r','MarkerSize',15)
axis('square')
axis([0.5 12.5 10 50])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('J');('F');('M');('A');('M');('J');('J');('A');('S');('O');('N');('D')] );
ylabel('Number of wet events')
legend('Obs.','Sim.')
subplot(2,1,2)
errorbar(1:12,mean(wetSumObs),std(wetSumObs),'.','MarkerSize',15)
hold on
errorbar(1:12,mean(wetSum),std(wetSum),'.r','MarkerSize',15)
axis('square')
axis([0.5 12.5 2000 1.8e4])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12] );
set(gca,'XTickLabel',[('J');('F');('M');('A');('M');('J');('J');('A');('S');('O');('N');('D')] );
ylabel('Total duration of wet events [min]')
export_fig Figs\tif\Evaluate_wet_events.tif -tif -transparent -r300
savefig(gcf,['Figs\Evaluate_wet_events.fig']); close all
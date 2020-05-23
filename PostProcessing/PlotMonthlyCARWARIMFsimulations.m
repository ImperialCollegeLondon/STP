function [ ] = PlotMonthlyCARWARIMFsimulations( MaternCov , Precipitation , Cloud , Copula , numR )
%PlotMonthlyCARWARIMFsimulations Plot n realiziations versus observed WAR,
%CAR and IMF data
figure('units','normalized','outerposition',[0 0 0.5 1])
figure('units','normalized','outerposition',[0 0 0.5 1])
figure('units','normalized','outerposition',[0 0 0.5 1])
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
for i=1:12
    clear war imf car tsWARFIMA
    %% Plot WAR
    parfor I=1:numR
        [ tsWARFIMA{I} ] = MultiTriVARFIMA( MaternCov , Precipitation , Cloud , Copula , i , size(Precipitation.WAR(i).WAR,1) , 1 );
    end
    figure(1)
    for I=1:numR
        war(:,I)=tsWARFIMA{I}.WARreal;
    end
    subplot(4,3,i)
    [n1,c1] = hist(Precipitation.WAR(i).WAR,0:0.1:1);
    n1=n1/sum(n1);
    h1=bar(c1,n1);
    hold on
    [n2,c2] = hist(war,0:0.1:1);
    n2=n2./size(war,1);
    h2=plot(c2,mean(n2')','.red','markers',16);
    h3=plot(c2,min(n2')','--black');
    h4=plot(c2,max(n2')','--black');
    axis([-0.05 1.05 0 max(max(n1,max(max(n2))))]);
    title([datestr(datenum(01,i,1),'mmm')])
    %% Plot IMF
    figure(2)
    for I=1:numR
        imf(:,I)=tsWARFIMA{I}.IMFreal;
    end
    subplot(4,3,i)
    [n1,c1] = hist(Precipitation.IMF(i).IMF,0:0.5:5);
    n1=n1/sum(n1);
    h1=bar(c1,n1);
    hold on
    [n2,c2] = hist(imf,0:0.5:5);
    n2=n2./size(imf,1);
    h2=plot(c2,mean(n2')','.red','markers',16);
    h3=plot(c2,min(n2')','--black');
    h4=plot(c2,max(n2')','--black');
    axis([-0.5 5.5 0 max(max(n1,max(max(n2))))]);
    title([datestr(datenum(01,i,1),'mmm')])
    %% Plot CAR
    % wet
    figure(3)
    for I=1:numR
        car(:,I)=tsWARFIMA{I}.CARreal;
    end
    subplot(4,3,i)
    [n1,c1] = hist(Cloud.Wet(i).CAR,0:0.1:1);
    n1=n1/sum(n1);
    h1=bar(c1,n1);
    hold on
    [n2,c2] = hist(car,0:0.1:1);
    n2=n2./size(car,1);
    h2=plot(c2,mean(n2')','.red','markers',16);
    h3=plot(c2,min(n2')','--black');
    h4=plot(c2,max(n2')','--black');
    axis([-0.05 1.05 0 max(max(n1,max(max(n2))))]);
    title([datestr(datenum(01,i,1),'mmm')])
end
figure(1)
export_fig Figs\tif\WAR_stats.tif -tif -transparent -r300
figure(2)
export_fig Figs\tif\IMF_stats.tif -tif -transparent -r300
figure(3)
export_fig Figs\tif\CAR_stats.tif -tif -transparent -r300
figure(1)
savefig(gcf,'Figs\WARsimulatedMonthly.fig');
figure(2)
savefig(gcf,'Figs\IMFsimulatedMonthly.fig');
figure(3)
savefig(gcf,'Figs\CARsimulatedMonthlyW.fig');
close all
end
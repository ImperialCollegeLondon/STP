% clear;clc;

% filename = 'AnimatedSTREAPRain.gif';
% load('H:\CODE_MATLAB\AWE-GEN-2D\Files\Output_Pr_Clt_500UV_dAR_linearAR.mat','simDATA','simIMF','simWAR')
% figurePath = 'H:\CODE_MATLAB\TestRainfall\IMAGES_linearAR';
% load('H:\CODE_MATLAB\AWE-GEN-2D\Files\Output_Pr_Clt_500UV_dAR_old.mat','simDATA','simIMF','simWAR')
% figurePath = 'H:\CODE_MATLAB\TestRainfall\IMAGES';

load('H:\CODE_MATLAB\AWE-GEN-2D\Files\Output_Pr_Sim_01.mat','simDATA','simIMF','simWAR')
figurePath = 'H:\CODE_MATLAB\TestRainfall\IMAGES_cpm';


mkdir(figurePath)
%% create a bunch of images and then use ffpmeg to create video
for MON = [2]
h = figure;
setFigureProperty('Paper');
XYWH = [50,50,240,432];
axis tight manual % this ensures that getframe() returns a consistent size

monthDATA = extractOnePeriod(simDATA,simDATA.Time.Month == MON);
monthIMF = simIMF(simDATA.Time.Month == MON);
monthWAR = simWAR(simDATA.Time.Month == MON);
Istart = 1;
for I = 1:2000%length(monthDATA.Time)
    R = originalData(extractOnePeriod(monthDATA,I));
    timn = monthDATA.Time(I);
    timn.Format = '20**/MM/dd HH:mm';
    R(R<=0) = NaN;
    if STATS(extractOnePeriod(monthDATA,I),'war')==0% ~any(~isnan(R(:)))
        Istart = [];
        continue;
    elseif isempty(Istart)
        Istart = I;
    end
    hp = uipanel('Title','CurrentTime(in 110Kmx110Km domain)','FontSize',12,...
        'BackgroundColor','white',...
        'Position',[0 0.44 1 .56]);
    axes(hp)
    pcolor(monthDATA.XX,monthDATA.YY,R);hold on;
    shading flat
    % cptcmap('GMT_no_green', 'mapping', 'direct');
    % plot(UKMap.borderE,UKMap.borderN,'k-','linewidth',1);
    % 'precip_meteoswiss'%'precip2_17lev'
    cptcmap('precip_meteoswiss', 'mapping', 'direct');
    % xlim([-300,800]*1000);
    % ylim([0,1200]*1000)
    % caxis([0,30])
    % caxis([0,15])
    axis equal
    axis off
    c = colorbar('Fontsize',8);%('location','Manual', 'position', [0.83 0.5 0.02 0.31]);
    c.TickLabels = strcat(c.TickLabels,'mm/h');
    % hold off
    text(390*1000,216*1000,sprintf('%s',timn),'color','k',...
        'fontweight', 'bold','fontsize',8 )
    
    hp = uipanel('Title','SimIMF','FontSize',12,...
        'BackgroundColor','white',...
        'Position',[0 0.22 1 .22]);
    axes(hp)
    plot(monthDATA.Time([Istart:I]),monthIMF(Istart:I));hold on;
    plot(monthDATA.Time(I),monthIMF(I),'ro')
    ax = gca;
    ax.XLim(2) = monthDATA.Time([I]);
    % ax.XLim(2).Day = ax.XLim(2).Day+1;
    ax.XLim(2).Hour = ceil(ax.XLim(2).Hour/3)*3+3;ax.XLim(2).Minute = 0;
    ax.YLim(1) = 0;
    ax.YLim(2) = ceil(ax.YLim(2));
    hp = uipanel('Title','SimWAR','FontSize',12,...
        'BackgroundColor','white',...
        'Position',[0 0 1 .22]);
    axes(hp)
    plot(monthDATA.Time([Istart:I]),monthWAR(Istart:I));hold on;
    plot(monthDATA.Time(I),monthWAR(I),'ro')
    ax = gca;
    ax.XLim(2) = monthDATA.Time([I]);
    % ax.XLim(2).Day = ax.XLim(2).Day+1;
    ax.XLim(2).Hour = ceil(ax.XLim(2).Hour/3)*3+3;ax.XLim(2).Minute = 0;
    ax.YLim(1) = 0;
    ax.YLim(2) = ceil(ax.YLim(2));
    drawnow
    savePlot([figurePath,filesep,sprintf('Mon%02d-%04d',MON,I)],...
        'units','points','XYWH',XYWH,'onlyPng',true,'needreply','N');
    clf(h)
end
close all
end
%% create GIF
UKMap = getUKMap();
h = figure;
setFigureProperty('Paper');
XYWH = [50,50,480,500];
set(gcf,'units','points','position',XYWH,'defaultTextFontSize',20);
axis tight manual % this ensures that getframe() returns a consistent size
ntag = 1;

for i = I:I+24*12
    timn = datetime(2016,3,1)+minutes(5*(i-1));
    timn.Format = 'yyyyMMddHHmm';
    R = squeeze(simRain(i,:,:));
    R(R<0.03) = NaN;
    % hold on
    pcolor(XX,YY,R);hold on;
    shading flat
    % plot(UKMap.borderE,UKMap.borderN,'k-','linewidth',1);
    % 'precip_meteoswiss'
    % cptcmap('precip2_17lev', 'mapping', 'direct');
    % xlim([-300,800]*1000);
    % ylim([0,1200]*1000)
    caxis([0,30])
    axis equal
    axis off
    c = colorbar;%('location','Manual', 'position', [0.83 0.5 0.02 0.31]);
    c.TickLabels = strcat(c.TickLabels,'mm/h');
    % hold off
    timn=datetime(2016,3,1)+i/24/12;
    % text(490*1000,-200*1000,sprintf('%s',timn),'color','k',...
    %     'fontweight', 'bold','fontsize',8 )
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ntag == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    savePath = 'H:\CODE_MATLAB\SpatialTemporalDATA\RadarImages';
    saveName = sprintf('%03d',i);
    savePlot([savePath,filesep,saveName],'units','points','XYWH',XYWH,...
        'needreply','N','onlyPng',true);
    clf(h)
    ntag = ntag+1;
    if ntag>24*7.
        continue;
    end
    pause(0.05);
    % axis off
    
end

filename = 'AnimatedSTREAPRain.gif';
XX = DATA.XX;
YY = DATA.YY;
% [XX,YY] = meshgrid([1:13]*1000,[1:13]*1000);
UKMap = getUKMap();
h = figure;
setFigureProperty('Paper');
XYWH = [50,50,480,500];
set(gcf,'units','points','position',XYWH,'defaultTextFontSize',20);
axis tight manual % this ensures that getframe() returns a consistent size
ntag = 1;

for I = 12*24*6*30:12*24*12*30
    R = squeeze(simRain(I,:,:));
    R(R<0.03) = 0;
    if nansum(R(:)) ~= 0
        break;
    end
end
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
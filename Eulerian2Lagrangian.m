% 
% This script is to Estimate ARMA parameters for temporal evolution of
% precipitation field.
%
% Input: 
%       A bunch of Radar Images;
% Output:
%       Transformed Lagrangian Coordinate System;
% 
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk
%

%% Configuration
clear;clc

[OBSDATA] = getData();

ObsSUMMARY = struct;
ObsSUMMARY.meanTimeseries = nanmean(OBSDATA,'timeseries');
ObsSUMMARY.meanSpatial = nanmean(OBSDATA,'Spatial');
ObsSUMMARY.WAR = STATS(OBSDATA,'war');

%% Split Storms and calculate velocity

EventSUMMARY = struct;

[ObsSUMMARY.radarIsEvent] = findRadarEvents(ObsSUMMARY.WAR,OBSDATA.dt);
[Events] = getStorms(OBSDATA,ObsSUMMARY.radarIsEvent);

EventSUMMARY.velocity = arrayfun(@(event)Advection(event,0.1),Events);
EventSUMMARY.velocity = arrayfun(@(vel)smooth(vel,6),EventSUMMARY.velocity);

% select Storms advected with a constant velocity (approximately)
EventSUMMARY.okTag = selectStorms(EventSUMMARY.velocity,OBSDATA.dt);
calEvents = Events(EventSUMMARY.okTag);
speed = getStormSpeed(EventSUMMARY.velocity(EventSUMMARY.okTag));


%% Transform to Lagrangian System
dim = 500;
[calEvents500_Eulerian] = expandEventsDomain(calEvents,dim);
[calEvents500_Lagrangian] = transLagrangian(calEvents500_Eulerian,speed);

dim = 110;
calEvents_Lagrangian = cutEventsDomain(calEvents500_Lagrangian,dim);

warning on
save('H:\CODE_MATLAB\calEvents500_Eulerian.mat','calEvents500_Eulerian','-v7.3');
save('H:\CODE_MATLAB\calEvents500_Lagrangian.mat','calEvents500_Lagrangian','-v7.3');
clear calEvents500_Eulerian calEvents500_Lagrangian

save('H:\CODE_MATLAB\calEvents.mat','calEvents','speed','EventSUMMARY','ObsSUMMARY');
save('H:\CODE_MATLAB\calEvents_Lagrangian.mat','calEvents_Lagrangian','-v7.3');



%% visualise
unit = 1000;
XX = calEvents500_Eulerian.XX;
YY = calEvents500_Eulerian.YY;
for evi = 1:numel(calEvents500_Lagrangian)
    rain = originalData(calEvents500_Eulerian(evi),'double');
    RE = originalData(calEvents500_Lagrangian(evi),'double');
    for i = round(size(RE,3)/2) % i = 1:size(RE,1)
        subplot(1,2,1)
        R = squeeze(rain(:,:,i)); R(R==0) = NaN;
        pcolor(XX/unit,YY/unit,R);shading flat;cptcmap('rain-mmh')
        title('Eulerian Coor')
        subplot(1,2,2)
        R = squeeze(RE(:,:,i)); R(R==0) = NaN;
        pcolor(XX/unit,YY/unit,R);shading flat;cptcmap('rain-mmh')
        title('Lagrangian Coor')
        pause(0.1)
    end
    pause(1)
    savePlot(['Figs\Lagrangian\event',num2str(evi),'-timestep',num2str(i)],...
        'units','centimeters','XYWH',[5,5,20,8],'onlyPng',true,...
        'needreply','N');
    close(gcf);
end


%% AUXILLARY FUNCTION

function [qFields] = transLagrangian(calEvents,calSpeed)

dt = 5;
dx = 1;
qFields = arrayfun(@(data,speed)toLagrangian(data,speed),calEvents,calSpeed);

    function newData = toLagrangian(data,speed)
        rain = originalData(data,'double');
        rain = single(rain);
        simU = speed.Vx*1000/3600;% m/s
        simV = speed.Vy*1000/3600;
        
        U=0; V=0;
        qField = [];
        for i = 1:size(rain,3)
            field = squeeze(rain(:,:,i));
            u=(simU(i)/1000)*60*dt/dx;
            v=(simV(i)/1000)*60*dt/dx;
            field=circshift(field,round([-v+V -u+U]));
            V=-v+V;
            U=-u+U;
            qField{i} = field;
        end
        qField = cat(3,qField{:});
        newData = data;
        newData.Val = qField;
        newData = squeezeData(newData,'int16');
    end
end

function [calEvents_central] = cutEventsDomain(calEvents,dim)

calEvents_central = arrayfun(@(data)cutDomain(data,dim),calEvents);

    function PRS = cutDomain(data,dim)
        x1 = round(dim/2);
        cendim = 250-x1:250+(dim-x1-1);
        isInArea = ismember(data.XX,data.XX(cendim,cendim)) ...
            & ismember(data.YY,data.YY(cendim,cendim));
        PRS = extractOneArea(data,isInArea,data.XX(cendim,cendim));
    end
end


function [calEvents_bigger] = expandEventsDomain(calEvents,dim)

calEvents_bigger = arrayfun(@(data)expandDomain(data,dim),calEvents);

    function PRS = expandDomain(data,dim)
        unit = 1000;
        x_yr = data.XX(1,1)-(dim/2-55)*unit:1000:data.XX(1,1)+(dim/2+54)*unit;
        y_yr = data.YY(1,1)-(dim/2-55)*unit:1000:data.YY(1,1)+(dim/2+54)*unit;
        [XX,YY] = meshgrid(x_yr,y_yr);
        [PRS,~] = importNIMROD_P(XX,YY,data.Time);
    end
end

function speed = getStormSpeed(velocity)

% speed = arrayfun(@(x)updateAdvection(x,...
%     nanmean(x.Vx),nanmean(x.Vy),nanmean(x.Vstd),nanmean(x.Vdir)),velocity);
speed = velocity;

end

function okTag = selectStorms(velocity,dt)

Vstd = arrayfun(@(x)nanmean(x.Vstd),velocity);
Vdir = arrayfun(@(x)nanmean(x.Vdir),velocity);
Dur = arrayfun(@(x)length(x.Vdir)*dt/60,velocity);%[hour]

thre_VEstd = nanmean(Vstd);
thre_VEdir = nanmean(Vdir);
okTag = (Vstd<thre_VEstd) & (Vdir<thre_VEdir) & (Dur>3);

end


function radarIsEvent = findRadarEvents(WAR,dt)

minDur = 60;% min;
sepThre = 2;% hour
% seperate & find all events

[lengthi,stormsi] = seperateRadarEvents( WAR, 0.02, sepThre, minDur, dt);
EvNo = [zeros(1,stormsi(1)-1),cell2mat(arrayfun(@(si,nsi,li,evi)...
    [evi*ones(1,li),0*ones(1,nsi-si-li)],stormsi,...
    [stormsi(2:end),length(WAR)+1],lengthi,1:length(stormsi),'UniformOutput', false))];

radarIsEvent = EvNo;% with Tag of eventNo

end

function [Event] = getStorms(obsDATA,isEvent)

Event = RainfallDataClass();
tag = 1;
for evi = reshape(unique(isEvent),1,[])
    if evi~=0
        Event(tag,1) = extractOnePeriod(obsDATA,isEvent==evi);
        tag = tag+1;
    end
end

end

function [DATA] = getData()

load('H:\CODE_MATLAB\PRS_Birm_2013.mat','DATA');

end



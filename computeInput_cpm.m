% ------------------------------------------------------------------- %
% This script is to compute all required input file for STREAP module.
% where calibration is based on Climate Model
% study area: Clee Hill
% 110 Km * 110 Km
% @ Yuting Chen
% ------------------------------------------------------------------- %

%% configuration
region = getfield(REGIONS_info(),'Cleehill');
inputInfo = getInputInfo(region,'cpm-historical');
[Time,inputInfo] = getTime(inputInfo);
inputInfo.time = Time;
catchName = 'Birmingham_cpm';
mkdir(catchName);
%% file 1: MeanArealStats
[WAR,IMF] = deal([]);
tic
[DATA] = getData(inputInfo,Time);
save([catchName,'/CPMData.mat'],'DATA','-v7.3')
%%
[WAR,IMF] = getStormArrivalStats(inputInfo,DATA);
toc

[U500,V500] = getUV(inputInfo);
[CAR] = NaN;
for ensNo = 1:length(WAR)
    MeanArealStats(ensNo) = struct('WAR',WAR{ensNo}','IMF',IMF{ensNo}',...
        'Time',Time','CAR',CAR','U500',U500','V500',V500');
end
save([catchName,'/MeanArealStats.mat'],'MeanArealStats');

%% file 1: 'CV.mat' 
% Rainfall coefficient of variation [-], monthly data
[CV] = getCV(inputInfo,DATA);
save([catchName,'/CV.mat'],'CV');

%% file 2: 'ARMA.mat' 
% ARMA coefficients
% # time consuming #
% Data has been saved in 'Birmingham/ARMA.mat'
% including: 'ARMA'; 'mean_spatial_correlation','Rains'
edit('Eulerian2Lagrangian.m')
edit('Lagrangian2ARMA.m')
%% file 3: 'expoSpatialCorrelation.mat'
% # time consuming #
% Exponential coefficient of the spatial correlation
h = 0.1;
[covariance,distance] = getSpatialCov(inputInfo,DATA,h);
save([catchName,'/SpatialCov_ensNoH1.mat'],'covariance','distance','-v7.3')

% pool all ensemble together into one file.
cov = cell(12,1);
for ensNo = 1:length(covariance)
    for mon = 1:12
        cov{mon} = cat(2,cov{mon},covariance{ensNo}{mon});
    end
end
covariance = cov;
distance = distance{1};
save([catchName,'/SpatialCov.mat'],'covariance','distance','-v7.3')

% find out suitable expSpatialCorrelation
load([catchName,'/SpatialCov.mat'],'covariance','distance')
[expoSpatialCorr,cb_theo] = getExpoSpatialCorr(inputInfo,covariance,distance,h);
expoSpatialCorr = -1./expoSpatialCorr;% transform into: @(expCorr)exp(s*expCorr);(s:unit:km);
save([catchName,'/expoSpatialCorrelation.mat'],'expoSpatialCorr')

%% file 4: 'NHRO.mat' 
% Non homogenic rainfall ocurrence; 
% in this example number of days above a certain threshold, 
% can be also hours or minutes
[NHRO] = getNHRO(inputInfo,DATA);
save([catchName,'/NHRO.mat'],'NHRO')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% file 5: 'ObsAnnualRain.mat'
% ObsAnnualRain - Gridded annual rainfall [mm] {Y}(i,j).
% 1980-2017 CEH-GEAR was used
[ObsAnnualRain] = getObsAnnualRain(inputInfo,DATA);
save([catchName,'/ObsAnnualRain.mat'],'ObsAnnualRain')

%% file 6: 'thresholdI.mat'
[ thresholdI ] = I_Threshold_estimation(catchName,inputInfo);
save([catchName,'/thresholdI.mat'],'thresholdI')

%% file 7: 'gpFitSO.mat' 
% Generalized Pareto parameters for the observed and simulated daily data
% 
% For Engelberger catchment- GP parameters were derived from daily MeteoSwiss gridded rainfall data,
% 2-km resolution, for the period 1981-2012.
% A 30 years simulation was generated and analyzed for the simulated GP parameres, 
% see 'Rainfall_accumulation_estimation' for example
[gpFitO,gpFitS] = Rainfall_accumulation_estimation(catchName,inputInfo);
save([catchName,'/gpFitSO.mat'],'gpFitO','gpFitS')

%% AUXILLARY Func

function [Time,inputInfo] = getTime(inputInfo)
Time = arrayfun(@(x)datetime(x,1,1):hours(1):datetime(x,1,1)+days(12*30)-hours(1),...
    inputInfo.year(2:end),'UniformOutput',false);
Time = cat(2,[],Time{:});
Time([end-(30*24-1):end]) = [];
PTime0 = datetime(inputInfo.year(1),12,1):hours(1):datetime(inputInfo.year(1),12,1)+hours(30*24-1);
Time = [PTime0,Time];
end

function [DATA] = getData(inputInfo,Time)
% load all years, each ensembles, all months
DATA = cell(12,1);

% added because of the strange saving period of ukcp cpm 2.2: starting from 1980.12
for mon = 1:12
    fileName = ['H:\CODE_MATLAB\Cleehill\',sprintf('Ensems_mon%02d.mat',mon)];
    load(fileName,'RainEnsembles');
    for ensNo = 1:length(RainEnsembles)
        R = RainEnsembles{ensNo};
        R = permute(R,[2,1,3]);%%%%%%%%%%%%%%%%%%%%%%%%
        if mon == 12
            R = cat(3,R,NaN([size(R,[1,2]),24*30]));
        else
            R = cat(3,NaN([size(R,[1,2]),24*30]),R);
        end
        R = reshape(R,[size(R,[1,2]),30*24,size(R,3)/30/24]);
        DATA{ensNo} = cat(3,DATA{ensNo},R);
    end
end
for ensNo = 1:length(RainEnsembles)
    DATA{ensNo} = reshape(DATA{ensNo},...
        size(DATA{ensNo},1),size(DATA{ensNo},2),[]);
    DATA{ensNo}(:,:,[1:30*24*11]) = [];
    DATA{ensNo}(:,:,[end-(30*24-1):end-0]) = [];
    DATA{ensNo} = RainfallDataClass(DATA{ensNo},-1,32,'mm/h',Time',...
        inputInfo.XX,inputInfo.YY,'');
end
end

function [WAR,IMF] = getStormArrivalStats(inputInfo,DATA)
[WAR,IMF] = deal(cell(length(DATA),1));
for ensNo = 1:length(DATA)
    WAR{ensNo} = STATS(DATA{ensNo},'war');
    IMF{ensNo} = STATS(DATA{ensNo},'imf');
    % exclude typical measuring errors
    % WAR{ensNo}(IMF{ensNo}>10) = NaN;
    % IMF{ensNo}(IMF{ensNo}>10) = NaN;
end
end
function [U,V] = getUV(inputInfo)

[era5info,E,N] = defineEra5DataBound(inputInfo.year,inputInfo.UV.pressure,inputInfo.Er,inputInfo.Nr);
[U,V,time,missval,dt] = getERA5(era5info);
[lia1,~] = ismember(time,inputInfo.time);
U = U(:,:,lia1);V = V(:,:,lia1);time = time(lia1);
[U,V] = formulateUV(inputInfo,U,V,time,missval,dt);
    function [era5info,E,N] = defineEra5DataBound(year,pressure,Er,Nr)
        
        spaceEN = [Er,Nr];
        
        era5info = struct;
        filePath = 'K:\DATA_FCS\ERA5_Birm\ERA5_pressure';
        fileName = sprintf('era5_%dhPa_Birmingham_all_1979_2018.nc',pressure);
        era5info.fileN = [filePath,filesep,fileName];
        % A = ncinfo(datainfo.fileN);
        
        % ADDITIONAL GRID
        % fileName = sprintf('era5_%dhPa_Birmingham_all_additional_1979_2018.nc',pressure);
        % era5info.fileN = [filePath,filesep,fileName];
        % Aadd = ncinfo(datainfo.fileN);
        
        fprintf(['Notice that for full datasets around $Birmingham$',...
            'might need to include $Add$\n'])
        
        LON = ncread(era5info.fileN,'longitude');
        LAT = ncread(era5info.fileN,'latitude');
        [E, N] = ll2os(LAT, LON);
        E = E/1000;
        N = N/1000;
        time = ncread(era5info.fileN,'time')/24+datenum(datetime(1900,1,1,0,0,0));
        era5info.timestart = find(datetime(datevec(time)).Year==year(1),1);
        era5info.timeend = find(datetime(datevec(time)).Year==year(end)+1,1);
        if isempty(era5info.timeend)
            era5info.t_len = Inf;
        else
            era5info.t_len = era5info.timeend-era5info.timestart+1;
        end
        era5info.loni = find(min(abs(E-spaceEN(1))) == abs(E-spaceEN(1)));
        era5info.lonlen = find(min(abs(E-spaceEN(2))) == abs(E-spaceEN(2)))-era5info.loni+1;
        era5info.lati = find(min(abs(N-spaceEN(4))) == abs(N-spaceEN(4)));
        era5info.latlen = find(min(abs(N-spaceEN(3))) == abs(N-spaceEN(3)))-era5info.lati+1;
        
        % Result Checking Module
        if any([era5info.loni,era5info.lonlen,era5info.lati,era5info.latlen]<=0)
            f = msgbox('Climate Data Format (N Axis) is changed! Please Contact Yuting to modify code', 'Error','error');
        end
        if any([E(:),N(:)] > 1e4)
            f = msgbox('Climate Data Format (E,N) is changed! Please Contact Yuting to modify code', 'Error','error');
        end
    end
    function [U,V,time,missval,dt] = getERA5(era5info)
        fileN = era5info.fileN;
        loni = era5info.loni;
        lati = era5info.lati;
        timei = era5info.timestart;
        lonlen = era5info.lonlen;
        latlen = era5info.latlen;
        timelen = era5info.t_len;

        getOneVar = @(fileN,var)ncread(fileN,var,[loni,lati,timei],[lonlen,latlen,timelen]);
        % wind speed: 500hPa: for motion of many tropical cyclones
        % wind speed: 700hPa: for shallower tropical cyclones
        U = getOneVar(fileN,'u');% 'U component of wind' 'eastward_wind' 'm s**-1'
        V = getOneVar(fileN,'v');% 'northward_wind' 'V component of wind' 'm s**-1'
        time = datetime(datevec(ncread(fileN,'time',timei,timelen)/24+datenum(datetime(1900,1,1,0,0,0))));
        missval = -32767;
        dt = hours(1);
    end
    function [U,V] = formulateUV(inputInfo,U,V,time,missval,dt)
        U(U == missval) = NaN;
        U = nanmean(reshape(U,[],size(U,3)),1);
        V(V == missval) = NaN;
        V = nanmean(reshape(V,[],size(V,3)),1);
        if inputInfo.dt<=dt
            U = interp1(time, U, inputInfo.time,'nearest');
            V = interp1(time, V, inputInfo.time,'nearest');
            U = fillmissing(U,'nearest');
            V = fillmissing(V,'nearest');
        else
            % # unfinished #
        end
    end
end
function [CV] = getCV(inputInfo,DATAs)
CV = cell(12,1);
for ensNo = 1:length(DATAs)
    DATA = DATAs{ensNo};
    for mon = 1:12
        if ~isempty(DATA)
            monthData = extractOnePeriod(DATA,DATA.Time.Month == mon);
            CV{ensNo}(mon,1) = STATS(monthData,'cvpos');
        else
            monthData = getData(inputInfo,inputInfo.year(1));
            monthData = extractOnePeriod(monthData,monthData.Time.Month == mon);
            for year = inputInfo.year(2:end)
                [DATA_temp] = getData(inputInfo,year);
                monthData = appendTime(monthData,...
                    extractOnePeriod(DATA_temp,DATA_temp.Time.Month == mon));
            end
            CV{ensNo}(mon,1) = STATS(monthData,'cvpos');
        end
    end
end
end
function [expoSpatialCorr,cb_theo] = getExpoSpatialCorr(inputInfo,covariance,distance,h)
[expoSpatialCorr,allSpatialCorr] = deal([]);
for mon = 1:12
    cb_obs = nanmedian(covariance{mon}(1:100,:),2);% use median not MEAN for skewed data!
    s = distance{1}(1:100);
    IniNo = 5;% increase for preventing local optima;
    opts = optimoptions(@fmincon,'Display','iter','TolX',10^-5,...
        'TolFun',10^-6,'MaxIter',10^5,'MaxFunEvals',10^5,'Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',...
        @(a)obj(a,cb_obs),'x0',10,'lb',1.1,'ub',500,'options',opts);
    ms = MultiStart('StartPointsToRun','bounds-ineqs','MaxTime',20*60);
    [expoSpatialCorr(mon,1),fval] = run(ms,problem,IniNo);
    cb_theo = getCB_THEO(expoSpatialCorr(mon,1));
    subplot(3,4,mon)
    plot(cb_obs,'ko');hold on;plot(cb_theo,'r-');
    title(getMonthName(mon))
    drawnow
end
    function res = obj(alpha_g,cb_obs)
        cb_theo = getCB_THEO(alpha_g);
        rmse = @(sim,obs)sqrt(nansum(((sim(:))-(obs(:))).^2));% sqrt(nanmean(((sim(:))-(obs(:))).^2));
        res = rmse(cb_theo,cb_obs);
    end
    function cb_theo = getCB_THEO(alpha_g)
        cb_theo = zeros(length(s),1);
        for i = 1:length(s)
            s_this = s(i);
            cb_theo(i,1) = func_cbs(alpha_g,s_this,h);
        end
    end
    function cb_s = func_cbs(alpha_g,s_this,h)
        % all input/output is scalar
        expFunc = @(s,alpha)exp(-s/alpha);% exp(-s/alpha_g);
        cg_s = expFunc(s_this,alpha_g);
        cg_0 = expFunc(0,alpha_g);
        f1 = @(a,b)exp(-(a.^2+b.^2-2*cg_s*a.*b)./(2*(1-cg_s.^2)))./...
            (2*pi*sqrt(1-cg_s^2));
        r = integral(@(y)exp(-y.^2/2/cg_0),h,Inf,'RelTol',1e-4,'AbsTol',1e-5)...
            /sqrt(2*pi*cg_0);
        ps = integral2(@(a,b)f1(a,b),h,Inf,h,Inf,'RelTol',1e-4,'AbsTol',1e-5);
        cb_s = ps-r^2;
    end
end
function [covariance,distance] = getSpatialCov(inputInfo,DATAs,h)
[r,disI,disJ,Ind_is,Ind_js] = deal([]);
[covariance,distance] = deal(cell(12,1));
if isempty(DATAs)
    % ...
else
    for mon = 1:12
        for ensNo = 1:length(DATAs)
            DATA = DATAs{ensNo};
            subplot(3,4,mon)
            monthData = extractOnePeriod(DATA,DATA.Time.Month == mon);
            [covariance{ensNo}{mon,1},distance{ensNo}{mon,1}] = getSpatialCov_mon(monthData);
            title(getMonthName(mon));
            drawnow
        end
        save('Birmingham/SpatialCov.mat','covariance','distance')
    end
end
    function [covariance,distance] = getSpatialCov_mon(monthData)
        covariance = [];
        binaryRain = originalData(monthData)>h;
        setStart = 1;
        thisWholeLen = size(binaryRain,3);
        setLen = thisWholeLen; % 4000; % depends on PC's RAM.
        for setNo = 1:ceil(thisWholeLen/setLen)
            setInd = setStart:1:min(setStart+setLen-1,thisWholeLen);
            unit = 1000;%m-km
            if ~isempty(setInd)
                bRain = binaryRain(:,:,setInd);
                typRainInd = STATS(extractOnePeriod(monthData,setInd),'WAR')>0.02;
                subData = extractOnePeriod(monthData,typRainInd);
                tic
                domainX = size(inputInfo.XX,1);
                [ii,jj] = meshgrid(-(domainX-1):(domainX-1));
                cod = [];
                for i = 1:length(subData.Time)
                    thisD = extractOnePeriod(subData,i);
                    data = originalData(thisD);
                    corBi = xcorr2(double(data>h))./(domainX-abs(ii))./(domainX-abs(jj));
                    cod = cat(3,cod,corBi);
                end
                dis = round(sqrt(ii.^2+jj.^2)*2.2);
                dis = repmat(reshape(dis,size(dis,1),size(dis,2),1),1,1,size(cod,3));
                ff = @(x)reshape(x,[prod(size(x,[1,2])),size(x,3)]);
                covariance_temp = grpstats(ff(cod),nanmean(ff(dis),2),{'mean'});
                toc
                distance = unique(reshape(squeeze(dis(:,:,1)),1,[]));
                % [covariance_temp,distance,r,disI,disJ,Ind_is,Ind_js]=spatialcov(...
                %     monthData.XX/unit,monthData.YY/unit,...
                %     bRain(:,:,typRainInd),monthData.dx,50,r,disI,disJ,...
                %     Ind_is,Ind_js);
                setStart = setInd(end)+1;
            end
        end
        covariance = cat(2,covariance,covariance_temp);
    end
end
function [NHRO] = getNHRO(inputInfo,DATAs)
% UKCP18 data was used (period:...)
% resolution: 1hour; 2.2Km
dailyThre = 1;
for ensNo = 1:length(DATAs)
    NHRO{ensNo} = struct('occurrence', cell(1, 12));
    DATA = DATAs{ensNo};
    for mon = 1:12
        for year = [unique(DATA.Time.Year)]'
            thisTime = DATA.Time.Year == year & DATA.Time.Month == mon;
            if any(thisTime)
                data0 = extractOnePeriod(DATA,thisTime);
                data0 = aggregate(data0,1,days(1)./minutes(DATA.dt),'mm/day');
                RAIN = originalData(data0);
                NHRO{ensNo}(mon).occurrence(end+1,:) = nansum(reshape(RAIN,[],size(RAIN,3))>dailyThre,2);
            end
        end
    end
end
end

function [ObsAnnualRain] = getObsAnnualRain(inputInfo,DATAs)
% UKCP18 data was used
% all ensemble were pooled together to represent the annula variability.
dailyThre = 1;
ObsAnnualRain = [];
yearTag = 1;
for ensNo = 1:length(DATAs)
    DATA = DATAs{ensNo};
    for year = [unique(DATA.Time.Year)]'
        thisTime = DATA.Time.Year == year;
        data0 = extractOnePeriod(DATA,thisTime);
        if length(data0.Time)>days(300)/minutes(data0.dt) % exclude the first year data, only december is available
            data0 = aggregate(data0,1,days(1)./minutes(DATA.dt),'mm/day');% unit:mm/day
            RAIN = originalData(data0);
            ObsAnnualRain{yearTag} = squeeze(nansum(RAIN,3));
            yearTag = yearTag + 1;
        end
    end
end
end






% ------------------------------------------------------------------- %
% This script is to compute all required input file for STREAP module.
% study area: Clee Hill
% 110 Km * 110 Km
% @ Yuting Chen
% ------------------------------------------------------------------- %

%% configuration
region = getfield(REGIONS_info(),'Westuk');
inputInfo = getInputInfo(region);

%% file 1: MeanArealStats
[Time,inputInfo] = getTime(inputInfo);
[DATA] = getData(inputInfo);
[WAR,IMF] = getStormArrivalStats(inputInfo,DATA);
[U500,V500] = getUV(inputInfo);
[CAR] = NaN;
MeanArealStats = struct('WAR',WAR','IMF',IMF','Time',Time','CAR',CAR',...
    'U500',U500','V500',V500');
save('Birmingham/MeanArealStats.mat','MeanArealStats');
%% file 1: 'CV.mat' 
% Rainfall coefficient of variation [-], monthly data
[CV] = getCV(inputInfo,DATA);
save('Birmingham/CV.mat','CV');
%% file 2: 'ARMA.mat' 
% ARMA coefficients
% # time consuming #
% Data has been saved in 'Birmingham/ARMA.mat'
% including: 'ARMA'; 'mean_spatial_correlation','Rains'
edit('Eulerian2Lagrangian.m')
edit('Lagrangian2ARMA.m')
%% file 3: 'expoSpatialCorrelation.mat') 
% Exponential coefficient of the spatial correlation
[covariance,distance] = getSpatialCov(inputInfo,DATA);
save('Birmingham/SpatialCov.mat','covariance','distance','-v7.3')
%%
[expoSpatialCorr] = getExpoSpatialCorr(inputInfo,covariance,distance);
%% file 4: 'NHRO.mat' 
% Non homogenic rainfall ocurrence; 
% in this example number of days above a certain threshold, 
% can be also hours or minutes


%% file 5: 'gpFitSO.mat' 
% Generalized Pareto parameters for the observed and simulated daily data


%% file 6: 'thresholdI.mat'










%% AUXILLARY Func



function inputInfo = getInputInfo(region)
inputInfo = struct;

inputInfo.year = 2013;
inputInfo.Er = [region.minE,region.minE+(region.dimE-1)*region.dx];
inputInfo.Nr = [region.minN,region.minN+(region.dimN-1)*region.dx];
inputInfo.dt = minutes(5);%min;
inputInfo.UV = struct;inputInfo.UV.pressure = 500;

inputInfo.datafile = 'H:\CODE_MATLAB\PRS_Birm_2013.mat';
end
function [Time,inputInfo] = getTime(inputInfo)
inputInfo.time = datetime(inputInfo.year,1,1):inputInfo.dt:...
    datetime(inputInfo.year+1,1,1)-inputInfo.dt;
Time = datenum(inputInfo.time);
end
function [DATA] = getData(inputInfo)
load(inputInfo.datafile,'DATA');
end
function [WAR,IMF] = getStormArrivalStats(inputInfo,DATA)
WAR = STATS(DATA,'war');
IMF = STATS(DATA,'imf');
end
function [U,V] = getUV(inputInfo)

[era5info,E,N] = defineEra5DataBound(inputInfo.year,inputInfo.UV.pressure,inputInfo.Er,inputInfo.Nr);
[U,V,time,missval,dt] = getERA5(era5info);
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
        era5info.timeend = find(datetime(datevec(time)).Year==year(end)+1,1)-1;
        era5info.t_len = era5info.timeend-era5info.timestart+1;
        
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
            U = interp1(time, U, inputInfo.time);
            V = interp1(time, V, inputInfo.time);
        else
            % # unfinished #
        end
    end
end
function [CV] = getCV(inputInfo,DATA)
CV = [];
for mon = 1:12
    monthData = extractOnePeriod(DATA,DATA.Time.Month == mon);
    CV(mon,1) = STATS(monthData,'cvpos');
end
end

function [covariance,distance] = getExpoSpatialCorr(inputInfo,covariance,distance)
[covariance,distance] = deal([]);
end

function [covariance,distance] = getSpatialCov(inputInfo,DATA)
[r,disI,disJ,Ind_is,Ind_js] = deal([]);
[covariance,distance] = deal([]);
for mon = 1:12
    subplot(3,4,mon)
    monthData = extractOnePeriod(DATA,DATA.Time.Month == mon);
    [covariance{mon,1},distance{mon,1}] = caliSpatialCoef(monthData);
    title(getMonthName(mon));
    drawnow
end
    function [covariance,distance] = caliSpatialCoef(monthData)
        covariance = [];
        h = 0.1;
        binaryRain = originalData(monthData)>h;
        setStart = 1;
        thisWholeLen = size(binaryRain,3);
        setLen = thisWholeLen; % 2000;
        for setNo = 1:ceil(thisWholeLen/setLen)
            tic
            setInd = setStart:1:min(setStart+setLen-1,thisWholeLen);
            unit = 1000;
            if ~isempty(setInd)
                [covariance_temp,distance,r,disI,disJ,Ind_is,Ind_js]=spatialcov(...
                    monthData.XX/unit,monthData.YY/unit,...
                    binaryRain(:,:,setInd),monthData.dx,50,r,disI,disJ,...
                    Ind_is,Ind_js);
                setStart = setInd(end)+1;
            end
            toc
        end
        covariance = cat(2,covariance,covariance_temp);
        drawnow
    end
end













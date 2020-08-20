function inputInfo = getInputInfo(region,datasource)
% input info for computing model input
arguments
    region
    datasource (1,:) char = 'radar'
end
inputInfo = struct;
inputInfo.Er = [region.minE,region.minE+(region.dimE-1)*region.dx];
inputInfo.Nr = [region.minN,region.minN+(region.dimN-1)*region.dx];
switch(datasource)
    case 'radar'
        inputInfo.year = 2007:2018;
        inputInfo.dt = minutes(5);%min;
        inputInfo.dx = 1;
        inputInfo.UV = struct;inputInfo.UV.pressure = 500;
        inputInfo.datafile = ['H:\DATA_RADAR\CLEEHILL_RADAR\',sprintf('PRS_CLEEHILL_%04d-%04d.mat',...
            inputInfo.year(1),inputInfo.year(end))];
        x_yr = inputInfo.Er(1):inputInfo.Er(end);
        y_yr = inputInfo.Nr(1):inputInfo.Nr(end);
    case 'cpm-historical'
        inputInfo.year = 1980:2000;
        inputInfo.dt = minutes(60);%min;
        inputInfo.dx = 2.2;
        inputInfo.UV = struct;inputInfo.UV.pressure = 500;
        inputInfo.datafile = ['H:\CODE_MATLAB\Cleehill\'];
        x_yr = inputInfo.Er(1):2.2:inputInfo.Er(end);
        y_yr = inputInfo.Nr(1):2.2:inputInfo.Nr(end);
    case 'cpm-future'
        % ###
    otherwise
end
[inputInfo.XX,inputInfo.YY] = meshgrid(x_yr*1000,y_yr*1000);
end
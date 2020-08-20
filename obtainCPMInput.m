
x_yr = inputInfo.Er(1):inputInfo.Er(end);
y_yr = inputInfo.Nr(1):inputInfo.Nr(end);
[XX,YY] = meshgrid(x_yr*1000,y_yr*1000);
region = getfield(REGIONS_info(),'Cleehill');
IntFac = 32; % to make the resolution = (Nimrod)

% [DATA,status] = importNIMROD_P(XX,YY,YEAR);
addpath(genpath('H:\CODE_MATLAB\SpatialTemporalDATA\UKCP18'));

ENSEMBLENO = getEnsNos();
ENSEMBLENO = ENSEMBLENO(1:12);

options = 'Save';
data = struct('Years',[1980,2000],...
    'fileGetPath','K:/UkCp18/',...
    'savePath',['H:\CODE_MATLAB\',region.Name]);
[RainEnsembles,IntFac] = readCPM_pr(region,ENSEMBLENO,1:12,IntFac,data,options);
% save(['H:\CODE_MATLAB\',sprintf('PRS_CLEEHILL_%04d.mat',YEAR)],'DATA','-v7.3');
% DATA = RainfallDataClass(PRS,-1,32,'mm/h',PTime',XX,YY,'');

rmpath(genpath('H:\CODE_MATLAB\SpatialTemporalDATA\UKCP18'));
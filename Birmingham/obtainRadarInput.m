x_yr = inputInfo.Er(1):inputInfo.Er(end);
y_yr = inputInfo.Nr(1):inputInfo.Nr(end);
[XX,YY] = meshgrid(x_yr*1000,y_yr*1000);
YEAR = 2013;
[DATA,status] = importNIMROD_P(XX,YY,YEAR);
save(inputInfo.datafile,'DATA','-v7.3');
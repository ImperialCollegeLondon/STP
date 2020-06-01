for YEAR = [2007:2012,2015:2018]
    x_yr = inputInfo.Er(1):inputInfo.Er(end);
    y_yr = inputInfo.Nr(1):inputInfo.Nr(end);
    [XX,YY] = meshgrid(x_yr*1000,y_yr*1000);
    [DATA,status] = importNIMROD_P(XX,YY,YEAR);
    save(['H:\CODE_MATLAB\',sprintf('PRS_CLEEHILL_%04d.mat',YEAR)],'DATA','-v7.3');
end

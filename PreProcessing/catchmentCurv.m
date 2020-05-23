function [ profc , planc ] = catchmentCurv( domainDTM )
%catchmentCurv Define curvature index for the catchment
%   Inputs:
%   domainDTM - catchment elevation [m] map
%   Outputs:
%   profc , planc - curvature parameters
%% Find curvature and plot
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
hsize=[12 12];
Sigma=1;
h = fspecial('gaussian', hsize, Sigma);
DTM2 = filter2(h,domainDTM,'same');
[profc,planc] = curvature(DTM2);
figure,
imagesc(planc);
axis('square')
set(gca,'YDir','normal');
caxis([-10 10])
title('Planform curvature'); colorbar;
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
end
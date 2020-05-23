function [ AccMatrix ] = catchmentAccumulation( domainDTM , AccThreshold )
%catchmentAccumulation A function to find the main channels in the
%catchment
%   Inputs:
%   domainDTM - elevation [m] map of the domain
%   AccThreshold - a lower threshold to define the main streams,
%   subjectively determined
%   Output:
%   AccMatrix - accumulation field index
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Find and plot stream accumulation
[ ~ , ~ , ~ , AccMatrix ] = FlowAcc( domainDTM , 100 );
AccMatrix=smoothn(AccMatrix);
AccMatrix=smoothn(AccMatrix);
figure,
h(1)=imagesc(AccMatrix);
set(gca,'YDir','normal');
axis('square')
caxis([0 AccThreshold]) 
colorbar
title('Flow accomulation [-]')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
figure,
h(2)=imagesc(domainDTM.*(AccMatrix<AccThreshold)); 
set(gca,'YDir','normal');
axis('square')
Acc_comp=colormap;
Acc_comp(1,:)=1;
colormap(Acc_comp);
colorbar
title('Catchment elevation [m] - Main streams in white')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
end
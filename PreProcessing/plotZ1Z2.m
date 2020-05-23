function [ ] = plotZ1Z2( Z1 , Z2 )
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Plot
figure
imagesc(Z1)
axis('square')
set(gca,'YDir','normal');
colorbar
title('Z1 - elevation surface [m]')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
figure
imagesc(Z2)
axis('square')
set(gca,'YDir','normal');
colorbar
title('Z2 - elevation surface [m]')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
end
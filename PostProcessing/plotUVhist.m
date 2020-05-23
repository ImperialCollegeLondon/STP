function [ h3 ] = plotUVhist( U , V , Uobs , Vobs )
%% Plot
h1=figure('units','normalized','outerposition',[0 0 1 1]);
scatterhist(Uobs,Vobs,'Marker','.');
axis([-50 50 -50 50])
title('Observed')
xlabel('U component [m s^{-1}]')
ylabel('V component [m s^{-1}]')

h2=figure('units','normalized','outerposition',[0 0 1 1]);
scatterhist(U,V,'Marker','.');
axis([-50 50 -50 50])
title('Simulated (downscaled)')
xlabel('U component [m s^{-1}]')
ylabel('V component [m s^{-1}]')

h3=figure('units','normalized','outerposition',[0 0 1 1]);
u1=uipanel('position',[0 0 0.5 1]);
u2=uipanel('position',[0.5 0 0.5 1]);

set(get(h1,'Children'),'parent',u1);
set(get(h2,'Children'),'parent',u2);
close(h1,h2)
end
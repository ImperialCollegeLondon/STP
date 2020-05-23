function [ h3 ] = plotUVrose( U , V , Uobs , Vobs )
%% Plot
[ObsD,ObsV] = cart2pol(Uobs,Vobs);
ObsD=wrapTo2Pi(ObsD);
ObsD=radtodeg(ObsD);
[SimD,SimV] = cart2pol(U(:,1),V(:,1));
SimD=wrapTo2Pi(SimD);
SimD=radtodeg(SimD);
%% Plot
h1=WindRose(ObsD,ObsV,'maxfrequency',8,'ndirections',72,'nfreq',8,'nspeeds',7,'freqlabelangle',45,'titlestring',{'Observed';' '},'lablegend','Advection velocity [m s^{-1}]','anglenorth',270,'angleeast',180,'speedround',10);
set(h1,'units','normalized','outerposition',[0 0 1 1]);

h2=WindRose(SimD,SimV,'maxfrequency',8,'ndirections',72,'nfreq',8,'nspeeds',7,'freqlabelangle',45,'titlestring',{'Simulated';' '},'lablegend','Advection velocity [m s^{-1}]','anglenorth',270,'angleeast',180,'speedround',10);
set(h2,'units','normalized','outerposition',[0 0 1 1]);

h3=figure('units','normalized','outerposition',[0 0 1 1]);
u1=uipanel('position',[0 0 0.5 1]);
u2=uipanel('position',[0.5 0 0.5 1]);

set(get(h1,'Children'),'parent',u1);
set(get(h2,'Children'),'parent',u2);
close(h1,h2)
end
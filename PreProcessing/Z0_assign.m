%% Setting roughness length [m] to th Switzerland land use map
% based on the 2004 land use, 17 classes partions (see table in the GIS
% folder)
z0(domainLU==1)=0.7;
z0(domainLU==2)=1.5;
z0(domainLU==3)=0.7;
z0(domainLU==4)=0.7;
z0(domainLU==5)=0.03;
z0(domainLU==6)=0.18;
z0(domainLU==7)=0.12;
z0(domainLU==8)=0.06;
z0(domainLU==9)=0.03;
z0(domainLU==10)=2.3;
z0(domainLU==11)=0.45;
z0(domainLU==12)=1.7;
z0(domainLU==13)=0.0002;
z0(domainLU==14)=0.0002;
z0(domainLU==15)=0.008;
z0(domainLU==16)=0.004;
z0(domainLU==17)=0.001;
function [ HZ , Z ] = HorizonAnglePolar( DTM , cellsize , dthe , maxDist )
% HorizonAnglePolar This function calculates the horizon angle for different directions
%   Inputs:
%   DTM - matrix digital elevation model [m]
%   cellsize - dimension cell [m]
%   dthe - directions intervals [deg]
%   maxDist - max search distance [m]
%	Outputs:
%	HZ - horizon angle array [angular degree] per direction
%	Z - azimuth directions [angualr degree] from N
%% Initilazing
Z=0:dthe:360-dthe; %% Search angle [angular degree]
iDTM=zeros(size(DTM));
for i=1:size(DTM,1)*size(DTM,2)
    iDTM(i)=i;
end
[dx,dy]=meshgrid((1:size(DTM,1))*cellsize,(1:size(DTM,2))*cellsize);
HZ=zeros(size(DTM,1)*size(DTM,2),length(Z));
%% Compute Horizon angle f(azimuth)
progressbar('Compute Horizon angle f(azimuth)') % Init single bar
for k=1:length(Z)
    rDTM=imrotate(iDTM,-Z(k)-180);
    for i=1:size(rDTM,1)
        for j=1:size(rDTM,2)
            if rDTM(i,j)>0
                idx=rDTM(1:i,j);
                idx(idx==rDTM(i,j)|idx==0)=[];
                idx=unique(idx);
                if ~isempty(idx)
                    dE=(DTM(idx)-DTM(rDTM(i,j)));
                    dL=sqrt((dx(rDTM(i,j))-dx(idx)).^2+(dy(rDTM(i,j))-dy(idx)).^2);
                    dE(dL>maxDist)=[];
                    dL(dL>maxDist)=[];
                    mang=atand(dE./dL);
                    mang(mang<0)=0;
                    HZ(rDTM(i,j),k)=min(90-mang);
                end
            end
        end
    end
    progressbar(k/length(Z)) % Update progress bar
end
HZ=reshape(HZ,size(DTM,1),size(DTM,2),length(Z));
HZ(HZ==0)=90;
end
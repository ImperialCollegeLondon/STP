function [ betaS , xiS ] = findSlope( domainDTM , dxdy )
%findSlope Find the maximum slope for each pixel in the downslope direction
%   Inputs:
%   domainDTM - elevation DTM [m]
%   dxdy - pixel length [m]
%   Output:
%   betaS - maximum slope for each pixel [deg]
%   xiS - terrain slope azimuth [deg from north]
%% Initilazing
betaS=zeros(size(domainDTM));
xiS=zeros(size(domainDTM));
%% Find maximum slope
for i=2:size(domainDTM,1)-1
    for j=2:size(domainDTM,2)-1
        for ik=[-1 0 1]
            for jk=[-1 0 1]
                if ik==0 && jk==0
                    continue
                end
                if ik==-1 && jk==-1
                    xi=45;
                elseif ik==1 && jk==1
                    xi=225;
                elseif ik==-1 && jk==1
                    xi=315;
                elseif ik==1 && jk==-1
                    xi=135;
                elseif ik==1 && jk==0
                    xi=180;
                elseif ik==-1 && jk==0
                    xi=0;
                elseif ik==0 && jk==1
                    xi=270;
                elseif ik==0 && jk==-1
                    xi=90;
                end
                s=sqrt(ik^2+jk^2)*dxdy;
                beta=atand((domainDTM(i,j)-domainDTM(i+ik,j+jk))/s);
                if beta<=betaS(i,j)
                    betaS(i,j)=beta;
                    xiS(i,j)=xi;
                end
            end
        end
    end
end
betaS=abs(betaS);
%% Boundary
betaS(1,:)=betaS(2,:);
betaS(end,:)=betaS(end-1,:);
betaS(:,1)=betaS(:,2);
betaS(:,end)=betaS(:,end-1);
xiS(1,:)=xiS(2,:);
xiS(end,:)=xiS(end-1,:);
xiS(:,1)=xiS(:,2);
xiS(:,end)=xiS(:,end-1);
% 
% betaS(1,:)=nan;
% betaS(end,:)=nan;
% betaS(:,1)=nan;
% betaS(:,end)=nan;
% betaSb=smoothn(betaS);
% xiSb=smoothn(xiS);
% betaS(1,:)=betaSb(1,:);
% betaS(end,:)=betaSb(end,:);
% betaS(:,1)=betaSb(:,1);
% betaS(:,end)=betaSb(:,end);
% xiS(1,:)=xiSb(1,:);
% xiS(end,:)=xiSb(end,:);
% xiS(:,1)=xiSb(:,1);
% xiS(:,end)=xiSb(:,end);
% xiS(xiS<0)=0;
% betaS(betaS<0)=0;
end
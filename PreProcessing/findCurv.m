function [ omegaC ] = findCurv( domainDTM , dxdy , eta )
%findCurv Find the mean curve for each pixel in the DTM. Following Liston
%and Elder (2006).
%   Inputs:
%   domainDTM - elevation DTM [m]
%   dxdy - pixel length [m]
%   eta - the curvature length scale [m]
%   Output:
%   omegaC - curvature
%% Initilazing
omegaC=zeros(size(domainDTM));
nP=round(eta/dxdy); % Number of pixels for the search
nP(nP==0)=1; % Set minimum search radius
%% Find curvature
for i=1+nP:size(domainDTM,1)-nP
    for j=1+nP:size(domainDTM,2)-nP
        a=(domainDTM(i,j)-0.5*(domainDTM(i,j-nP)+domainDTM(i,j+nP)))/(2*eta);
        b=(domainDTM(i,j)-0.5*(domainDTM(i-nP,j)+domainDTM(i+nP,j)))/(2*eta);
        c=(domainDTM(i,j)-0.5*(domainDTM(i+nP,j+nP)+domainDTM(i-nP,j-nP)))/(2*eta*sqrt(2));
        d=(domainDTM(i,j)-0.5*(domainDTM(i+nP,j-nP)+domainDTM(i-nP,j+nP)))/(2*eta*sqrt(2));
        omegaC(i,j)=0.25*(a+b+c+d);
    end
end
%% Boundaries
omegaC(1:1+nP,:)=nan;
omegaC(end-nP:end,:)=nan;
omegaC(:,1:1+nP)=nan;
omegaC(:,end-nP:end)=nan;
%% Scale to [-0.5 0.5]
tmp=reshape(omegaC,[],1);
tmp(isnan(tmp))=[];
[muhat,sigmahat] = normfit(tmp);
tmp=normcdf(omegaC,muhat,sigmahat);
omegaC=-0.5+tmp;
tmp=smoothn(omegaC);
tmp(tmp<-0.5)=-0.5;
tmp(tmp>0.5)=0.5;
omegaC(1:1+nP,:)=tmp(1:1+nP,:);
omegaC(end-nP:end,:)=tmp(end-nP:end,:);
omegaC(:,1:1+nP)=tmp(:,1:1+nP);
omegaC(:,end-nP:end)=tmp(:,end-nP:end);
end
function [ Z1 ] = defineZ1( domainDTM , AccMatrix , AccThreshold , planc , CurvThreshold , CurveAccDist )
%defineZ1 Compute Z1 field
%% Combine flow accumulation and curvature
a=AccMatrix>=AccThreshold;
b=planc<=CurvThreshold(2)&planc>=CurvThreshold(1);
c=single(a)+single(b);
c(c>0)=1;
c=bwmorph(c,'clean');
c=bwmorph(c,'open');
c=bwmorph(c,'spur');
c=bwmorph(c,'clean');
[x,y]=meshgrid(0:size(domainDTM,1)-1,0:size(domainDTM,2)-1);
parfor i=1:size(c,1)*size(c,2)
    if c(i)>0
        d=sqrt((x(i)-x).^2+(y(i)-y).^2).*a;
        d=reshape(d,[],1);
        d(d==0)=[];
        if min(d)>CurveAccDist
            c(i)=0;
        end
    end
end
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Plot
figure,
imagesc(domainDTM.*(c==0))
set(gca,'YDir','normal');
axis('square')
d=colormap;
d(1,:)=1;
colormap(d);
colorbar
title('Z1 mask surface')
set(gca,'XTick',[50 100 150 200 250] );
set(gca,'XTickLabel',[5 10 15 20 25] );
set(gca,'YTick',[50 100 150 200 250] );
set(gca,'YTickLabel',[5 10 15 20 25] );
xlabel('Distance [km]')
ylabel('Distance [km]')
%% Z1 field
Z1=domainDTM.*single(c);
Z1(Z1==0)=nan;
z1=Z1(Z1>0);
Z1(x==0&isnan(Z1))=mean(z1);
Z1(y==0&isnan(Z1))=mean(z1);
Z1(x==max(max(x))&isnan(Z1))=mean(z1);
Z1(x==max(max(y))&isnan(Z1))=mean(z1);
X=reshape(x,[],1);
Y=reshape(y,[],1);
z1=double(reshape(Z1,[],1));
X(isnan(z1))=[];
Y(isnan(z1))=[];
z1(isnan(z1))=[];
Z1=griddata(X,Y,double(z1),x,y);
end
function [ corrAlpha ] = fftAutoCorr( dim , alpha , dxdy , aTarget)
%fftAutoCorr Function finds the fft simulated parameters that correnspond to the spatial correlation exponential decay. The function is calibrated
%toward an inputed exponential coefficient decay
%   Inputs:
%   dim - domain dimension [dim1 , dim2] [pixels]
%   alpha - fft stochastic Gaussian correlation parameter range for fitting
%   [min max]
%   dxdy - pixels` length (must be orthogonal) [km]
%   aTarget - exponential coefficient of the spatial correlation decay 
%   Output:
%   corrAlpha - exponential 'a' variable and it alpha value
%% Initilazing
ft = fittype( 'exp(a*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.399346009457777;
corrAlpha=[]; ij=[];
%% Simulating the corrField
for ALPHA=alpha(1):1:alpha(2)
    D=[]; C=[]; 
    for i=1:100
        [ corrField ] = fieldGen( dim , ALPHA );
        a(:,i)=reshape(corrField,dim(1)*dim(2),1);
    end
    %% Preapering the data
    [x,y]=meshgrid(linspace(1,dxdy*size(corrField,1),size(corrField,1)),linspace(1,dxdy*size(corrField,1),size(corrField,1)));
    if ALPHA==alpha(1)
        for c=1:1000000
            if c==size(corrField,1)^2-1 || c==1000
                break
            end
            i=round(random('uniform',1,size(corrField,1)^2)); j=round(random('uniform',1,size(corrField,1)^2));
            if sqrt((x(i)-x(j))^2+(y(i)-y(j))^2)<=dim(1)*dxdy/2
                C(size(C,1)+1,1)=corr(a(i,:)',a(j,:)');
                D(size(D,1)+1,1)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
                ij(size(ij,1)+1,1)=i; ij(size(ij,1),2)=j;
            end
        end
    end
    if ALPHA>alpha(1)
        C=zeros(size(ij,1),1); D=zeros(size(ij,1),1); 
        for c=1:size(ij,1)
            C(c,1)=corr(a(ij(c,1),:)',a(ij(c,2),:)');
            D(c,1)=sqrt((x(ij(c,1))-x(ij(c,2)))^2+(y(ij(c,1))-y(ij(c,2)))^2);
        end
    end
    %% Exponential fit
    [xData, yData] = prepareCurveData( D(D<mean(D)), C(D<mean(D)) );
    corrFieldFit = fit( xData, yData, ft, opts );
    corrAlpha(size(corrAlpha,1)+1,1)=corrFieldFit.a;
    corrAlpha(size(corrAlpha,1),2)=ALPHA;
end
%% Finding the best match
corrAlpha=corrAlpha(find(min(sqrt((corrAlpha(:,1)-aTarget).^2))==sqrt((corrAlpha(:,1)-aTarget).^2)),:);
end
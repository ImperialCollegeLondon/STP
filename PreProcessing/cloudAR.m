function [ EstMdl , invS ] = cloudAR( clt , cltDis , disMiu , disStd , tThreshold , cldFit , numLags , plotPACF )
%% cloudAR A function that calculates the ARMA proces of the cloud time series.
%The time series is first normalized for the dry period and than estimated.
%% Initilazing
cltDis(cltDis>tThreshold)=tThreshold;
%% Find cloudiness for distance from a wet event
t=cltDis(cltDis>0);
c=clt(cltDis>0);
%% Normalized
clt_norm=(c-disMiu(t))./disStd;
%% ARMA
for i=0:4
    for j=0:2
        Mdl=arima(i,0,j);
        [~,~,LLF] = estimate(Mdl,clt_norm);
        NumActive = i+j+1;
        AIC(i+1,j+1) = aicbic(LLF,NumActive);
    end
end
%% Choose the model to use
[i,j]=find(min(min(AIC))==AIC);
Mdl=arima(i-1,0,j-1);
EstMdl = estimate(Mdl,clt_norm);
%% Plot Partial Autocorrelation Function (PACF)
if plotPACF
    [ S ] = simulate(EstMdl,size(c,1));
    normClt=S.*disStd+disMiu(t);
    normClt=normcdf(normClt,mean(normClt),std(normClt));
    try
        invS=cldFit.CAR.fit.icdf(normClt);
    catch
        jType=cldFit.type; % Johnson dist.
        invS=f_johnson_inv(normClt,cldFit.coef,jType);
    end
    invS(invS>1)=1;
    invS(invS<0)=0;
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    parcorr(clt_norm,numLags);
    title('Observed');
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    subplot(2,2,2)
    parcorr(S(:,1),numLags);
    title(['Simulated ARMA(',num2str(i-1),',',num2str(j-1),')']);
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    subplot(2,2,3)
    hist(c,0:0.01:1);
    axis([0 1 0 10000])
    axis square
    title('Observed cloud histogram');
    subplot(2,2,4)
    hist(invS,0:0.01:1);
    title('Simulated cloud histogram');
    axis([0 1 0 10000])
    axis square
end
end
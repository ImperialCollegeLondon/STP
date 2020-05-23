function [ EstMdl , unorm , vnorm , h ] = UVadvectionVAR( U , V , numLags , plotPACF )
%UVadvectionVAR Function compute vector autoregressive model for the
%advection quantiles time series.
%   Input:
%   U and V - observed u and v component time series [m s^-1]
%   numLags - the number of lags of the PACF to plot
%   Output:
%   EstMdl - vgxvarx object model
%   unorm and vnorm - mean and standard deviation of the normalized u and v
%   components
%   h - figure handle
%% Initilazing
[unorm(1) , unorm(2)]=normfit(U);
[vnorm(1) , vnorm(2)]=normfit(V);
%% Set
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 20);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
%% Quantile and normalize
u=normcdf(U,unorm(1),unorm(2));
v=normcdf(V,vnorm(1),vnorm(2));
u=norminv(u);
v=norminv(v);
%% VAR - AR 1-10
Mdl=vgxset('n',2,'nAR',1);
[EstMdl1,~,LLF1] = vgxvarx(Mdl,[u v]);
[~,NumActive1] = vgxcount(EstMdl1);
Mdl=vgxset('n',2,'nAR',2);
[EstMdl2,~,LLF2] = vgxvarx(Mdl,[u v]);
[~,NumActive2] = vgxcount(EstMdl2);
Mdl=vgxset('n',2,'nAR',3);
[EstMdl3,~,LLF3] = vgxvarx(Mdl,[u v]);
[~,NumActive3] = vgxcount(EstMdl3);
Mdl=vgxset('n',2,'nAR',4);
[EstMdl4,~,LLF4] = vgxvarx(Mdl,[u v]);
[~,NumActive4] = vgxcount(EstMdl4);
AIC = aicbic([LLF1 LLF2 LLF3 LLF4],[NumActive1 NumActive2 NumActive3 NumActive4]);
%% Choose the model to use
switch find(min(AIC)==AIC)
    case 1
        EstMdl=EstMdl1;
    case 2
        EstMdl=EstMdl2;
    case 3
        EstMdl=EstMdl3;
    case 4
        EstMdl=EstMdl4;
end
%% Plot Partial Autocorrelation Function (PACF)
if plotPACF
    [ S ] = vgxsim(EstMdl,size(u,1));
    h(1)=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    parcorr(u,numLags);
    title('U component - observed');
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    subplot(2,2,2)
    parcorr(S(:,1),numLags);
    title(['U component - VAR(',num2str(find(min(AIC)==AIC)),') simulated']);
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    subplot(2,2,3)
    parcorr(v,numLags);
    title('V component - observed');
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    subplot(2,2,4)
    parcorr(S(:,2),numLags);
    title(['V component - VAR(',num2str(find(min(AIC)==AIC)),') simulated']);
    ylabel('PACF')
    axis([0  numLags -1 1])
    axis square
    %% Plot CDFs of CDFs (obs vs. simulated)
    S=normcdf(S);
    h(2)=figure('units','normalized','outerposition',[0 0 0.5 1]);
    subplot(2,1,1)
    h1=cdfplot(U);
    hold on
    S(:,1)=norminv(S(:,1),unorm(1),unorm(2));
    h2=cdfplot(S(:,1));
    title('U component');
    set(h2,'color','red','LineStyle','--')
    legend([h1 , h2],'Obs.','Sim.','Location','Southeast')
    axis('square')
    subplot(2,1,2)
    h1=cdfplot(V);
    hold on
    S(:,2)=norminv(S(:,2),vnorm(1),vnorm(2));
    h2=cdfplot(S(:,2));
    title('V component');
    set(h2,'color','r','LineStyle','--')
    axis('square')
end
end
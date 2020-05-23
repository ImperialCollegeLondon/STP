function [ ] = PlotCloudDis( clt , cltDis , dt , tThreshold , Mdl , cldFit )
%% CloudWetEventDis1
t1=cltDis(cltDis>0);
c1=clt(cltDis>0);
x=dt:dt:max(t1);
for i=1:length(x)
    M1(i)=mean(c1(t1==x(i)));
    S1(i)=std(c1(t1==x(i)));
end
ul=find(tThreshold==x);
[xData, yData] = prepareCurveData( x(1:ul), M1(1:ul) );
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [-0.00463936380272764 0.000161408091677754 0.479245770646933 -5.89286370999931e-05];
cltMiu1 = fit( xData, yData, ft, opts );
cltStd1=mean(S1(1:ul));
%% cloudAR
cltDis2=cltDis;
cltDis2(cltDis2>tThreshold)=tThreshold;
[ S ] = simulate(Mdl,size(clt,1));
normClt=S.*cltStd1+cltMiu1(cltDis2);
normClt=normcdf(normClt,mean(normClt),std(normClt));
try
    invS=cldFit.CAR.fit.icdf(normClt);
catch
    jType=cldFit.fit.type;
    invS=f_johnson_inv(normClt,cldFit.fit.coef,jType);
end
invS(invS>1)=1;
invS(invS<0)=0;
%% CloudWetEventDis2
t2=cltDis(cltDis>0);
c2=invS(cltDis>0);
for i=1:length(x)
    M2(i)=mean(c2(t2==x(i)));
    S2(i)=std(c2(t2==x(i)));
end
[xData, yData] = prepareCurveData( x(1:ul), M2(1:ul) );
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [-0.00463936380272764 0.000161408091677754 0.479245770646933 -5.89286370999931e-05];
[cltMiu2] = fit( xData, yData, ft, opts );
cltStd2=mean(S2(1:ul));
%% Plot
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
h1=plot(x,M1);
hold on
h2=plot(cltMiu1(x(1):x(ul)),'r');
h3=plot(x,M2,'black');
h4=plot(cltMiu2(x(1):x(ul)),'--g','LineWidth',2);
title('Mean cloud cover');
ylabel('[-]')
xlabel('Time [min] to the closest rain event')
legend([h1,h2,h3,h4],'Observed',[type(cltMiu1),' obs. fit'],'Simulated',[type(cltMiu2),' sim. fit']);
patch([0;0;x(ul);x(ul)],[0;1;1;0],[0;0;0;0],'FaceColor','b','FaceAlpha',0.2);
axis([0  x(end) 0 1])
axis square
subplot(1,2,2)
h1=plot(x,S1);
hold on
h3=plot(x,S2,'black');
x=x(1:ul); y(1:ul)=cltStd1;
h2=plot(x,y,'r');
y(1:ul)=cltStd2;
h4=plot(x,y,'--g','LineWidth',2);
title('Standard deviation of cloud cover');
ylabel('[-]')
xlabel('Time [min] to the closest rain event')
legend([h1,h2,h3,h4],'Observed','Obs. Std.','Simulated','Sim. Std.');
x=dt:dt:max(t1);
patch([0;0;x(ul);x(ul)],[0;1;1;0],[0;0;0;0],'FaceColor','b','FaceAlpha',0.2);
axis([0  x(end) 0 1])
axis square
end
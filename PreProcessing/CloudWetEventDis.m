function [ cltMiu , cltStd ] = CloudWetEventDis( clt , cltDis , dt , tThreshold , isplot )
%CloudWetEventDis A function that calculates and plot the mean cloud and
%standard deviation as a function of distance (in time) from the closest
%wet event.
%   Inputs:
%   clt - cloud time series [fraction]
%   cltDis - time series of the time [min] for the closest wet event, 0 indicate wet event time step
%   dt - time interval [min]
%   tThreshold - time threshold [min] of the upper limit of the exponential
%   decay (see AWE-GEN "fair weather")
%   Outputs:
%   cltMiu and cltStd - a decay function of the mean and standard deviation
%   as a function of distnace [min] from the closest rain event
%% Find cloudiness for distance from a wet event
t=cltDis(cltDis>0);
c=clt(cltDis>0);
%% Mean ans std
x=dt:dt:max(t);
for i=1:length(x)
    M(i)=mean(c(t==x(i)));
    S(i)=std(c(t==x(i)));
end
%% Fits
ul=find(tThreshold==x);
[xData, yData] = prepareCurveData( x(1:ul), M(1:ul) );
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [-0.00463936380272764 0.000161408091677754 0.479245770646933 -5.89286370999931e-05];
% Fit model to data.
[cltMiu] = fit( xData, yData, ft, opts );
cltStd=mean(S(1:ul));
%% Plot
if isplot
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    h1=plot(x,M);
    hold on
    h2=plot(cltMiu(x(1):x(ul)),'r');
    title('Mean cloud cover');
    ylabel('[-]')
    xlabel('Time [min] to the closest rain event')
    legend([h1,h2],'Mean',[type(cltMiu),' fit']);
    patch([0;0;x(ul);x(ul)],[0;1;1;0],[0;0;0;0],'FaceColor','b','FaceAlpha',0.2);
    axis([0  x(end) 0 1])
    axis square
    subplot(1,2,2)
    h1=plot(x,S);
    hold on
    h2=plot(x(1:ul),repmat(cltStd,ul),'r');
    title('Standard deviation of cloud cover');
    ylabel('[-]')
    xlabel('Time [min] to the closest rain event')
    legend('Std.','Constant Std.');
    patch([0;0;x(ul);x(ul)],[0;1;1;0],[0;0;0;0],'FaceColor','b','FaceAlpha',0.2);
    axis([0  x(end) 0 1])
    axis square
end
end
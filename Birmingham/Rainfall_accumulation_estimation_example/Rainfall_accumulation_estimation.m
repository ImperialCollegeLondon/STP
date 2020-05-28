function [ gpFitO , gpFitS ] = Rainfall_accumulation_estimation( )
%%
%Rainfall_accumulation_estimation This is an example for estimating the GP
%distribution parameters for the daily observed and simulated data
%This requires the model Process section to be set up to the point of the
%non-homogenous probability of precipitation cccurrence. Obsreved daily
%rainfall for this example were derived from the Meteo Swiss gridded daily
%products. If better temporal resolution exists, the model should be
%modified accordingly.
%% Load data
load('CV.mat') % Rainfall coefficient of variation [-], monthly data
load('ARMA.mat') % ARMA coefficients
load('expoSpatialCorrelation.mat') % Exponential coefficient of the spatial correlation
load('NHRO.mat') % Non homogenic rainfall ocurrence; in this example nu,ber of days above a certain threshold, can be also hours or minutes
load('Precipitation.mat');
Precipitation.data=data;
Precipitation.WAR=WAR;
Precipitation.IMF=IMF;
load('triNormTrans.mat')
load('Cloud.mat')
Cloud=load('Cloud.mat','ar');
Cloud.cloudDis=cloudDis;
Cloud.Dry=Dry;
Cloud.Wet=Wet;
load('obsOCC.mat')
load('thresholdI.mat')
%% Running a 30 years simulation for the domain
PreProc_30y_run
%% Estimating Generalized Pareto parameters for the 30 years simulated data
[ gpFitS ] = gpFitR( dR , thresholdI );
%% Load observed daily rainfall data
load('dRo.mat')
%% Estimating Generalized Pareto parameters for the observed data
[ gpFitO ] = gpFitR( dRo , 1 ); % Daily threshold of 1 mm per day
end
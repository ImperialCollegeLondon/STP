% compare output with cpm


% load data
load([catchName,'/CPMData.mat'],'DATA')
obsDATA = DATA{1};clear DATA
load(['Files\',sprintf('Output_Pr_sim_%02d.mat',1)],'simDATA');
for I = 2:5
    D = load(['Files\',sprintf('Output_Pr_sim_%02d.mat',I)],'simDATA');
    simDATA = appendTime(simDATA,D.simDATA);
end

%% compute statistics
DATA = obsDATA;
DATASTATS_obs = computeSTATS(DATA);
DATA = simDATA;
DATASTATS_sim = computeSTATS(DATA);

%% plot
figure;
ha = subplot(1,2,1);
plotSTATS_annulMean(DATASTATS_obs,DATA,ha);
ha = subplot(1,2,2);
plotSTATS_annulMean(DATASTATS_sim,DATA,ha);

ha = figure;hold on
plotSTATS_cv(DATASTATS_sim,ha);
plotSTATS_cv(DATASTATS_obs,ha);
legend('sim','obs')

ha = figure;hold on
plotSTATS_monthMean(DATASTATS_sim,ha);
plotSTATS_monthMean(DATASTATS_obs,ha);
legend('sim','obs')

ha = figure;hold on;
plotSTATS_imf(DATASTATS_sim,ha)
plotSTATS_imf(DATASTATS_obs,ha)
legend('sim','obs')

ha = figure;hold on;
plotSTATS_war(DATASTATS_sim,ha)
plotSTATS_war(DATASTATS_obs,ha)
legend('sim','obs')

ha = figure;hold on;
plotSTATS_acf(DATASTATS_sim,ha)
plotSTATS_acf(DATASTATS_obs,ha)
legend('sim','obs')

function plotSTATS_annulMean(DATASTATS,DATA,ha)
% annual mean
pcolor(ha,DATA.XX,DATA.YY,DATASTATS.spatialMean);shading flat
cptcmap('precip_annualMeanUK','mapping','direct');
axis('square')
colorbar
end

function plotSTATS_cv(DATASTATS,ha);
% seasonalCV
plot(1:12,DATASTATS.cvpos);
axis('square')
end

function plotSTATS_monthMean(DATASTATS,ha);
% seasonalCV
plot(1:12,DATASTATS.monthMean);
axis('square')
end

function DATASTATS = computeSTATS(DATA)
DATASTATS = struct;
% spatial annual mean
DATASTATS.spatialMean = nanmean(DATA,'spatial')*365*24;
DATASTATS.arealMean = nanmean(DATA,'timeseries');
% seasonality
DATASTATS.monthMean = nanmean(DATA,'month')*24;
DATASTATS.monthCV = nanmean(DATA,'monthCV');
DATASTATS.cvpos = nanmean(DATA,'monthcvpos');
% war,imf
DATASTATS.war = STATS(DATA,'war');
DATASTATS.imf = STATS(DATA,'imf');
% dry wet period

% autocorrelation
DATASTATS.acf_imf = autocorr(DATASTATS.imf,100-1);
DATASTATS.acf_war = autocorr(DATASTATS.war,100-1);
end

function plotSTATS_imf(DATASTATS,ha)
ha = cdfplot(DATASTATS.imf);
ha.Color(4) = 0.8;ha.LineWidth = 1;
set(gca,'XScale','log')
xlim([0.1,10])
xlabel('IMF')
axis('square')
end

function plotSTATS_war(DATASTATS,ha)
ha = cdfplot(DATASTATS.war);
ha.Color(4) = 0.8;ha.LineWidth = 1;
xlim([0.01,1])
xlabel('WAR')
axis('square')
end

function plotSTATS_acf(DATASTATS,ha)
subplot(121)
plot(0:99,DATASTATS.acf_imf);hold on
title('ACF imf[-]')
axis('square')
subplot(122)
plot(0:99,DATASTATS.acf_war);hold on
title('ACF war[-]')
axis('square')
end


% This script is to check several differences between CPM input and Radar input
% 
% Including:
% CV: 
%     check if the monthly CV (1h) from radar and cpm is the same
%     differences between 5min CV and 1h CV
% IMF and WAR:
%     check
% ExpoCorrelation:
%     check !
% obsAnnualRain:
%     check differences between 12 members
% 
% @ Yuting
% Imperial College London


%% Compare CV (notice here the CV is the CV for posistive rainfall)
cpm = load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham_cpm\CV.mat');
rad = load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham_obs\CV_1h2d2Km.mat');
rad5min = load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham_obs\CV.mat');

figure
hold on;
plot(1:12,rad.CV,'k-','linewidth',2);
% plot(1:12,rad5min.CV,'r-','linewidth',2,'color',[0 0 0 0.2]);
plot(cell2mat(cpm.CV'),'k-','color',[0 0 0 0.2])
axis('square')
xlabel('Mon');
ylabel('CV');
xlim([1,12]);
xticks([1:12])
xticklabels(getMonthName(1:12,1,1))
export_fig Figs\preCheck\cpm_obs_compare_CV.tif -tif -transparent -r300

%% Check IMF and WAR
cpm = load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham_cpm\MeanArealStats.mat');
rad = load('H:\CODE_MATLAB\AWE-GEN-2D\Birmingham_obs\MeanArealStats.mat');
figure
subplot(121)%IMF
hold on;
for ensNo = 1:length(cpm.MeanArealStats)
    H = cdfplot(cpm.MeanArealStats(ensNo).IMF);H.Color = [0.5 0.5 0.5 0.2];H.LineWidth = 1;
end
H1 = cdfplot(aggregate(rad.MeanArealStats.IMF, hours(1)/minutes(5),'mean'));
H1.Color(4) = 0.8;H1.LineWidth = 1;
H2 = cdfplot(aggregate(rad.MeanArealStats.IMF, 1,'mean'));
H2.Color(4) = 0.8;H2.LineWidth = 1;
set(gca,'XScale','log')
xlim([0.1,10])
xlabel('IMF')
axis('square')
subplot(122)
hold on;
for ensNo = 1:length(cpm.MeanArealStats)
    H = cdfplot(cpm.MeanArealStats(ensNo).WAR);H.Color = [0.5 0.5 0.5 0.2];H.LineWidth = 1;
end
H1 = cdfplot(aggregate(rad.MeanArealStats.WAR, 60/5,'mean'));
H1.Color(4) = 0.8;H1.LineWidth = 1;
H2 = cdfplot(aggregate(rad.MeanArealStats.WAR, 1,'mean'));
H2.Color(4) = 0.8;H2.LineWidth = 1;
xlim([0.01,1])
xlabel('WAR')
axis('square')
legend([H,H1,H2],{'cpm','rad','rad(5min1km)'},'Location','SouthEast')
export_fig Figs\preCheck\cpm_obs_compare_imfwar.tif -tif -transparent -r300
%%
figure;hold on;
subplot(121)
qqplot(aggregate(rad.MeanArealStats.WAR, 60/5,'mean'),rad.MeanArealStats.WAR);
axis('square')
xlim([0,1]);ylim([0,1]);xlabel('WAR(1h)');ylabel('WAR(5min)')
subplot(122)
qqplot(aggregate(rad.MeanArealStats.IMF, 60/5,'mean'),rad.MeanArealStats.IMF);
axis('square');xlabel('IMF(1h)');ylabel('IMF(5min)')
export_fig Figs\preCheck\cpm_obs_compare_imfwar_5min_1h.tif -tif -transparent -r300

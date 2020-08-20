clear;clc
close all
pack

dt = 5;

severalSites_sim = getSites('streap');
severalSites_obs = getSites('radar');

%%
for index = 1:numel(severalSites_sim.XX)
    
    oneSite_sim = extractOneArea(severalSites_sim,index,severalSites_sim.XX(index));
    oneSite_obs = extractOneArea(severalSites_obs,index,severalSites_obs.XX(index));
    
    [intensity_sim{index},returnT,duration] = getIDF(oneSite_sim,minutes(dt));
    [intensity_obs{index},returnT,duration] = getIDF(oneSite_obs,minutes(dt));

end

%%
figure;hold on;
for index = 1:numel(severalSites_sim.XX)
    subplot(3,3,index)
    plot(hours(duration),intensity_obs{index}','-')
    hold on;
    plot(hours(duration),intensity_sim{index}','--')
    set(gca,'YScale','log')
    ylim([2,2e2]);xlim([0,25])
    grid on
end
set(gca,'YScale','log')

%% AUXILLARY FUNCTION

function severalSites = getSites(datasource)
severalSites = [];
for I = 1:12
    
    if strcmp(datasource,'streap')
        load(['H:\DATA_STREAP\',sprintf('Output_Pr_sim_%02d.mat',I)],'simDATA','simIMF','simWAR','simU','simV');
        
        region = getfield(REGIONS_info(),'Cleehill');
        inputInfo = getInputInfo(region,'radar');
        simDATA.XX = inputInfo.XX;
        simDATA.YY = inputInfo.YY;
        
    elseif strcmp(datasource,'radar')
        load(['H:\DATA_RADAR\CLEEHILL_RADAR\',sprintf('PRS_CLEEHILL_%04d.mat',2007-1+I)],'DATA');
        simDATA = DATA;
    end
    XX = simDATA.XX;
    
    simDATA.Time.Year = 2007-1+I;
    isInArea = logical(zeros(size(XX)));
    isInArea([1,55,110],[1,55,110]) = 1;
    thisYear = extractOneArea(simDATA,isInArea,XX(isInArea));
    if isempty(severalSites)
        severalSites = thisYear;
    else
        severalSites = appendTime(severalSites,thisYear);
    end
    
end
end

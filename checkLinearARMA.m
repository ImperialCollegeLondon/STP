%% check linear ARMA

DECAY_TEST = -1./linspace(10,200,20);

% plot_OneStorm(DECAY_TEST);

[simACF,obsACF] = computeACF(DECAY_TEST);

bias = checkBias(simACF,obsACF,DECAY_TEST);

save('K:\DATA_ARMA_TRIALANDERROR\outputBias.mat','bias','simACF',...
    'obsACF','DECAY_TEST','-v7.3')
%%
load('K:\DATA_ARMA_TRIALANDERROR\outputBias.mat','bias','DECAY_TEST')
plotBias(bias,DECAY_TEST)
function plotBias(Ys,X)
Y = nanmean(Ys,2);
plot(X,Y,'linewidth',2);
hold on
plot(X(Y ==  min(Y)),min(Y),'ro','linewidth',2);
set(gca,'XScale','log')
grid minor
xlim([-0.1,-0.005])
ax = gca;
ax.XTick = X([1,3,6,12,20]);
ax.XTickLabels = ax.XTick;
xlabel('A');
ylabel('ObjFunc');
end
%% AUXILLARY FUNCTION
function bias = checkBias(simACF,obsACF,DECAY_TEST)

bias = [];
for expNo = 1:length(DECAY_TEST)
    bias(expNo,:) = fun1(simACF{expNo},obsACF);
end
    function bias = fun1(simACF,obsACF)
        rmse = @(sim,obs)sqrt(nanmean((sim-obs).^2));
        for stormi = 1:numel(simACF)
            bias(1,stormi) = rmse(nanmean(simACF{stormi},1),nanmean(obsACF{stormi},1));
        end
    end
end

function [simACF,obsACF] = computeACF(DECAY_TEST);

load('K:\DATA_ARMA_TRIALANDERROR\inputObs.mat','tsMatrix')
simACF = getSim(tsMatrix);
obsACF = getObs(tsMatrix);

    function acf = getObs(tsMatrix)
        D = load('H:\CODE_MATLAB\PRS_CLEEHILL_2015.mat');
        simDATA = D.DATA;
        acf = computeACF(simDATA,tsMatrix);
    end

    function acf = getSim(tsMatrix)
        acf = [];
        for expNo = 1:20
            load(['K:\DATA_ARMA_TRIALANDERROR\expNo',sprintf('%03d.mat',expNo)],...
                'simDATA');
            acf{expNo} = computeACF(simDATA,tsMatrix);
        end
    end

    function acf = computeACF(simDATA,tsMatrix)
        simStorm = [];
        acf = [];
        for i = 1:size(tsMatrix,1)
            stormInd = tsMatrix(i,3):tsMatrix(i,4);
            simStorm{i} = originalData(extractOnePeriod(simDATA,stormInd));
            % squeeze(simRain(stormInd,:,:));
            simStorm{i} = reshape(simStorm{i},[],size(simStorm{i},3));
            simStorm{i} = fillmissing(simStorm{i}, 'previous');
            acf{i} = plotACF(simStorm{i},[],false);
        end
    end

end

function plot_OneStorm(DECAY_TEST)
load('K:\DATA_ARMA_TRIALANDERROR\inputObs.mat')
h = figure;
cmap = cptcmap('blue_continuous','ncol',20);
for expNo = 1:20
    load(['K:\DATA_ARMA_TRIALANDERROR\expNo',sprintf('%03d.mat',expNo)],...
        'simDATA');
    simStorm = [];
    for i = 1:size(tsMatrix,1)
        stormInd = tsMatrix(i,3):tsMatrix(i,4);
        simStorm{i} = originalData(extractOnePeriod(simDATA,stormInd));
        % squeeze(simRain(stormInd,:,:));
        simStorm{i} = reshape(simStorm{i},[],size(simStorm{i},3));
        simStorm{i} = fillmissing(simStorm{i}, 'previous');
        subplot(2,1,1)
        plotACF(simStorm{i},cmap(expNo,:),true);hold on;
        title('Autocorrelation');
        xlabel('Lags');ylabel('Acorr')
        ylim([-0.1,1])
        subplot(2,1,2)
        plot(STATS(extractOnePeriod(simDATA,stormInd),'imf'));hold on;
        title('Areal Mean Precipitation');
        xlabel('Time steps[5min]');ylabel('IMF[mm/h]')
        drawnow
        % clf(h)
    end
end
end
function acf = plotACF(rain,colorM,pl)
acf = [];
for sitei = 1:size(rain,1)
    try
        if size(rain,2)>20
            [acf(sitei,:),lags,bounds] = autocorr(rain(sitei,:),20);
        else
            [acf_aux,lags,bounds] = autocorr(rain(sitei,:),size(rain,2)-1);
            acf(sitei,1:length(acf_aux)) = acf_aux;
        end
    catch
        acf(sitei,:) = NaN(1,21);
    end
end
if pl
    plot(lags,nanmean(acf,1),'color',colorM);
end
end

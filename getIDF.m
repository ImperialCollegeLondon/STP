function [intensity,returnT,duration] = getIDF(oneSite,dt)
% dt: <duration> 
% oneSite: <double> or <RainfallDataClass>

if strcmp(class(oneSite),'RainfallDataClass')
     %time = oneSite.Time;
    % oneSite = originalData(oneSite);
end

returnT = [5;20;100];
duration = [minutes([5,30]),hours([1:24])];

returnT = reshape(returnT,[],1);
index = 1;
intensity = [];

for dur = duration
    rainfallTS = aggregate(oneSite,1,dur/dt,'mm/h');
    time = rainfallTS.Time;
    % rainfallTS = squeeze(originalData(rainfallTS));
    rainfallTS = aggregate(squeeze(originalData(oneSite)),dur/dt,'mean');% use this one to prevent issue of imresize3 in aggregation
    am = getAM(rainfallTS,time);
    intensity(:,index) = getReturnLevels(am,returnT);
    index = index+1;
end

end

function am = getAM(rainfallTS,time)
am = grpstats(rainfallTS,time.Year,'max');
end
function intensity = getReturnLevels(am,rt)
options = statset('gevfit');
options.MaxIter = 500;
P = 1-(1./rt);
warning off
[parmhat,parmci] = gevfit(am,0.05,options);
intensity = gevinv(P,parmhat(1),parmhat(2),parmhat(3));
[warnMsg,warnId] = lastwarn;
if ~isempty(warnMsg)
    warnMsg = '';
    [parmhat,parmci] = evfit(am);% 
    intensity = evinv(P,parmhat(1),parmhat(2));
end
end
function [ dryPool , wetPool ] = drywetPool( II , data )
%drywetPool A pool of dry and wet events to form the yearly periods
%   Inputs:
%   II - number of realizations
%   data - Precipitation data package
%   Outputs:
%   dryPool and wetPool - the pool of wet and dry durations
for m=1:12
    try
        dryPool{m}=f_johnson_rnd(data(m).dryFitG.coef,data(m).dryFitG.type,200*II,1);
    catch
        dryPool{m}=data(m).dryFitG.random(200*II,1);
    end
    try
        wetPool{m}=data(m).wetFitG.random(200*II,1);
    catch
        wetPool{m}=f_johnson_rnd(data(m).wetFitG.coef,data(m).wetFitG.type,200*II,1);
    end
    dryPool{m}(dryPool{m}<1)=1;
    wetPool{m}(wetPool{m}<1)=1;
    dryPool{m}=ceil(dryPool{m});
    wetPool{m}=floor(wetPool{m});
end
end
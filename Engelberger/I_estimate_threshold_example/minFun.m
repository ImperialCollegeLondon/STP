function [ minX ] = minFun( thresholdI , RFIELD , objFun , idx )
[ simNday ] = countDays( thresholdI , RFIELD , idx );
minX=mean2(sqrt((simNday-objFun(idx)).^2));
end
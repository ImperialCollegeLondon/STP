function [ simNday ] = countDays( thresholdI , RFIELD , idx )
simNday=0;
for i=1:30
    tmp2=RFIELD{i};
    for d=1:365
        tmp3(d,:,:)=sum(squeeze(tmp2(1+288*(d-1):288*d,:,:)).*5/60);
    end
    tmp3(tmp3(:,idx)<thresholdI,idx)=0;
    simNday=simNday+squeeze(sum(tmp3(:,idx)>0));
end
simNday=simNday./i;
end
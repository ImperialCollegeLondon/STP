function [ simNdays ] = countDays2( thresholdI , RFIELD )
simNdays=zeros(13,13);
for i=1:30
    tmp2=RFIELD{i};
    tmp=zeros(13,13);
    for d=1:365
        tmp3(d,:,:)=sum(squeeze(tmp2(1+288*(d-1):288*d,:,:)).*5/60);
    end
    for j=1:size(thresholdI,1)*size(thresholdI,2)
        tmp3(tmp3(:,j)<thresholdI(j),j)=0;
    end
    simNdays=simNdays+squeeze(sum(tmp3>0));
end
simNdays=simNdays./i;
end
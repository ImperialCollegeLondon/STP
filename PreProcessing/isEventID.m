function [ IsEvent ] = isEventID( IsEvent )
c=1;
for i=2:length(IsEvent)
    if IsEvent(i)==1
        IsEvent(i)=c;
    end
    if IsEvent(i)==0 && IsEvent(i-1)>0
        c=c+1;
    end
end
end
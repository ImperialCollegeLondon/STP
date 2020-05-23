function [ tsMonth , tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , dt )
%% Initilazing
lastDay=eomday(2001,1:12); % The last day in a month without leap years
N=sum(lastDay)/dt; % Length of ts
tsMonth=zeros(1,N); % For each time step record it relevant month
tsWetDry=zeros(1,N); % For each time step mark 0 for dryG and 1 for wetG
for i=2:12
    lastDay(i)=lastDay(i)+lastDay(i-1);
end
lastDay=lastDay/dt;
for i=12:-1:1
    tsMonth(1:lastDay(i))=i;
end
dryEvents=[];
wetEvents=[];
%% Build annualy events time series (Starts - January 1, not a leap year)
%% First Event
% Start with a dry event
a=dryPool{1}(1);
dryPool{1}(1)=[];
tsWetDry(1:a)=0;
I=a+1;
dw=1; % Boolean value to determine if the next event will be dryG or wetG
M=1;
%% Constructing the ts for the rest of the year
while I<N
    switch dw
        case 0
            idx=dryPool{M}+I<=lastDay(M);
            idx=find(idx==1,1);
            if isscalar(idx)
                a=dryPool{M}(idx);
                dryPool{M}(idx)=[];
                tsWetDry(I:I+a-1)=0;
                I=I+a;
            else
                a=lastDay(M)-I;
                tsWetDry(I:I+a)=0;
                I=I+a+1;
            end
            I(I>N)=N;
            M=tsMonth(I);
            dw=1; 
        case 1
            idx=wetPool{M}+I<=lastDay(M);
            idx=find(idx==1,1);
            if isscalar(idx)
                a=wetPool{M}(idx);
                wetPool{M}(idx)=[];
                tsWetDry(I:I+a-1)=1;
                I=I+a;
            else
                a=lastDay(M)-I;
                tsWetDry(I:I+a)=1;
                I=I+a+1;
            end
            I(I>N)=N;
            M=tsMonth(I);
            dw=0;
    end
end
end
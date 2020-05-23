%% Processing
%% Initilazing
rng('shuffle');
t=(datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,23,55,0))'; % Dummy - One year dates based on the year 2001
I=30; % Number of years to simulate
%% Storm arrival process
[ dryPool , wetPool ] = drywetPool( I , Precipitation.data );
for I=1:30
    [ tsMonth , tsWetDry , wetPool , dryPool ] = tsWetDryMinG( wetPool , dryPool , 10/(24*60) );
    for i=1:length(tsWetDry) % For this example the 5-min temporal resolution is interpolated from the 10-min resolution
        tsWetDry5Min(i*2-1)=tsWetDry(i);
        tsWetDry5Min(i*2)=tsWetDry(i);
        tsMonth5Min(i*2-1)=tsMonth(i);
        tsMonth5Min(i*2)=tsMonth(i);
    end
    tsWetDry=tsWetDry5Min;
    tsMonth=tsMonth5Min;
    clear tsMonth5Min tsWetDry5Min
    [ tsMatrix ] = tsTable( tsWetDry , tsMonth );
    %% Advection - constant 5 m s^-1 toward the East
    simU=ones(tsMatrix(end),1).*5;
    simV=zeros(tsMatrix(end),1);
    simU=single(simU);
    simV=single(simV);
    %% Varfima - simulating WAR, CAR and IMF
    simWAR=zeros(tsMatrix(end),1);
    simIMF=zeros(tsMatrix(end),1);
    simCAR=zeros(tsMatrix(end),1);
    [ simDis ] = WetEventDis( simCAR , tsWetDry , 5 );
    for i=1:size(tsMatrix,1)
        m=tsMatrix(i,2);
        if tsMatrix(i,1)==1
            [ tsWARFIMA ] = MultiTriVARFIMA( MaternCov , Precipitation , Cloud , Copula , m , tsMatrix(i,4)-tsMatrix(i,3)+1 , 5 );
            E=abs(0-tsWARFIMA.WARreal(1,:))+abs(0-tsWARFIMA.WARreal(end,:));
            E=find(min(E)==E);
            if length(E)>1
                E=E(1);
            end
            simWAR(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA.WARreal(:,E);
            simCAR(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA.CARreal(:,E);
            simIMF(tsMatrix(i,3):tsMatrix(i,4))=tsWARFIMA.IMFreal(:,E);
            simCAR(simCAR<simWAR)=simWAR(simCAR<simWAR);
        end
    end
    %% Generating rain fields
    rField=zeros(105120,13,13);
    for m=1:12
        mmin(m,1)=min(tsMatrix(tsMatrix(:,2)==m,3));
        mmax(m,1)=max(tsMatrix(tsMatrix(:,2)==m,4));
        [ QField(mmin(m,1):mmax(m,1),:,:) ] = quantileFieldGen( [26 26] , ARMA , 5 , Precipitation.data(m).SpatialAlpha , [mmin(m,1) mmax(m,1)] , simU , simV , 2 , 5 );
        %% Non-homogeneous probability of precipitation cccurrence
        for i=mmin(m,1):mmax(m,1)
            [ rField(i,:,:) ] = NHPO( squeeze(QField(i,:,:)) , simWAR(i) , obsOCC{m} );
        end
        %% Assigning rainfall intensity for the Gaussian field
        for i=1:size(tsMatrix,1)
            M=tsMatrix(i,2);
            if tsMatrix(i,1)==1 && m==M
                [ rField(tsMatrix(i,3):tsMatrix(i,4),:,:) ] = invLN2( rField(tsMatrix(i,3):tsMatrix(i,4),:,:) , simWAR(tsMatrix(i,3):tsMatrix(i,4)) , simIMF(tsMatrix(i,3):tsMatrix(i,4)) , CV(m));
            end
        end
    end
    %% Saving the data
    SIMIMF{I}=simIMF;
    RFIELD{I}=rField;
    disp([num2str(I),'/30'])
end
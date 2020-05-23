function [ Alpha , Beta ] = fitAODpar( m , lambda , AOD )
%fitAODpar Fit the Angstorm turbidity parameters
%   Inputs:
%   m - month
%   lambda - wavelength
%   AOD - aerosol optical depth data
%% Initilazing
Alpha=[]; Beta=[];
f = fittype('b*x^(-a)');
%% Fit two parameters
for i=1:size(AOD.Data,2)
    if month(AOD.Time(i))==m
        if any(isnan(AOD.Data(:,i)))
            clear DATA LAMBDA
            c=1;
            for j=1:7
                if ~isnan(AOD.Data(j,i))
                    DATA(c,1)=AOD.Data(j,i);
                    LAMBDA(c,1)=lambda(j,1);
                    c=c+1;
                end
            end
            if length(DATA)>=6
                fit1 = fit(LAMBDA,DATA,f,'StartPoint',[1.3 0.05]);
                if fit1.a>=0 && fit1.a<=2.5 && fit1.b>=0 && fit1.b<=1
                    Alpha(size(Alpha,1)+1,1)=fit1.a;
                    Beta(size(Beta,1)+1,1)=fit1.b;
                end
            end
        else
            fit1 = fit(lambda,AOD.Data(:,i),f,'StartPoint',[1.3 0.05]);
            if fit1.a>=0 && fit1.a<=2.5 && fit1.b>=0 && fit1.b<=1
                Alpha(size(Alpha,1)+1,1)=fit1.a;
                Beta(size(Beta,1)+1,1)=fit1.b;
            end
        end
    end
end
end
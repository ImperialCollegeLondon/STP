function [ ARMA ] = ARMA_estimate( R, timeDim, pl )
%%
%ARMA_estimate This is an example for estimating the ARMA parameters for the Lagrangian system of coordinates
%This example is for a SINGLE storm and given only as an example and code suggestion.
%For further details see Paschalis (2013):
%% 
% Input: R: [time, loc1, loc2]
% Output: ARMA
%%
% <http://e-collection.library.ethz.ch/view/eth:7209 MODELLING THE SPACE-TIME STRUCTURE OF PRECIPITATION AND ITS IMPACT ON BASIN RESPONSE, DISS. ETH No. 21112>.

%% input checking
if timeDim == 1
    % ok
elseif timeDim == 3
    R = permute(R,[3,1,2]);
else
    error('Please Check Dimension order of R.');
end
%% Initilazing
lag=15; % Number of lags
BS=1; % Block size (increase for a larger domain)
%% Find temporal correlation
temp_autocorr(1:lag+1,1:size(R,1))=NaN;
for i=1:size(R,1)
    % Estimation of the Langrangian correlation
    for tlag=0:lag
        if tlag+i<size(R,1)
            % Identify corelation on the space-time distrortion
            im1=BlockMean(squeeze(R(i,:,:)),BS); % Convert to boolean value (1-rain, 0-no rain)
            % 
            % # Update on func 'BlockMean' #
            % In original func blockmean, im1/im2 is double when R is double,
            % is boolean when R is not double.
            % However, isa(R,'single') is taken into account and output became
            % to boolean.
            %
            % I think for <single> rainfall data, output should also be 
            % <single>/<double>.
            % (because in the model, ARMA coef is used for estimating value, 
            % instead of for estimating rain/no-rain).
            %
            % 2020.05.25
            % @ Yuting
            % 
            im2=BlockMean(squeeze(R(i+tlag,:,:)),BS);
            war1=mean2(R(i,:,:)>0);
            war2=mean2(R(i+tlag,:,:)>0);
            if war1>.1 && war2>.1 % For larger domains the treshold should increase
                try
                    taut=normxcorr2(im1,im2); % Normalized 2D cross-correlation
                    temp_autocorr(tlag+1,i)=nanmax(taut(:));
                catch
                    temp_autocorr(tlag+1,i)=NaN;
                end
            else
                temp_autocorr(tlag+1,i)=NaN;
            end
        end
    end
end
av_lang_corr=nanmean(temp_autocorr(:,1:size(R,1)),2);
av_lang_corr(1)=1;
mean_spatial_correlation(1,:)=av_lang_corr; % For multiply storms analysis, each storm is stored in this array
mean_spatial_correlation(2,:)=mean_spatial_correlation(1,:); % The storm is repeated for the sake of example
%% ARMA estimation
c=1;
for i=2:4 % Auto-regressive order
    for j=0:4 % Moving average order
        [ ar_coefs , th ] = Estimate_ARMA_coef_corr( nanmean(mean_spatial_correlation) , i , j , length(nanmean(mean_spatial_correlation)));
        ARMA(c).ar=ar_coefs{1};
        ARMA(c).ma=ar_coefs{2};
        ARMA(c).th=th;
        ARMA(c).rmse=sum((nanmean(mean_spatial_correlation)-th(1:end-1)').^2)/length(nanmean(mean_spatial_correlation));
        c=c+1;
    end
end
%% Plot
if pl
    figure('units','normalized','outerposition',[0 0 1 1]);
    c=1;
    tempx = (0:lag)*5;
    for i=2:4
        for j=0:4
            subplot(3,5,c)
            plot(tempx,nanmean(mean_spatial_correlation),'marker','.','linestyle','none','markerfacecolor','black');
            hold on;
            tempxx = (0:(length(ARMA(c).th)-1))*5;
            plot(tempxx,ARMA(c).th,'linewidth',1,'color','r');
            box on
            xlabel('Lag [min]');
            ylabel('ACF [-]')
            legend({'Observed','Fit'});
            title(['ARMA(',num2str(i),',',num2str(j),') RMSE=',num2str(ARMA(c).rmse)]);
            c=c+1;
            ylim([0.6,1])
        end
    end
    c=find(min([ARMA.rmse])==[ARMA.rmse]);
    set(subplot(3,5,c),'Color',[135/255 206/255 250/255]);
end

%% Nested functions
%% BlockMean
    function Y = BlockMean(X,V)
        % 2D block mean over 1st and 2nd dim [MEX]
        % Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de
        W=V;
        S = size(X);
        M = S(1) - mod(S(1), V);
        N = S(2) - mod(S(2), W);
        if M * N == 0
            Y = X([]);
            return;
        end
        MV = M / V;
        NW = N / W;
        XM = reshape(X(1:M, 1:N, :), V, MV, W, NW, []);
        if isa(X, 'double') || isa(X, 'single')
            Y = sum(sum(XM, 1), 3) ./ (V * W);
        elseif uint8(0.8) == 1
            Y = uint8(sum(sum(XM, 1), 3) ./ (V * W));
        else
            Y = uint8(round(sum(sum(XM, 1), 3) ./ (V * W)));
        end
        S(1) = MV;
        S(2) = NW;
        Y    = reshape(Y, S);
    end
%% Estimate ARMA coefficient correlation
    function [ coeffs , th ] = Estimate_ARMA_coef_corr( corr_fun , p , q , lag_num )
        corr_fun = corr_fun(:)';
        cl = length(corr_fun);
        x0 = [zeros(p,1);zeros(q,1)];
        options = optimset('Display','off');
        lb = [];
        ub = [];
        arma_coeffs = fmincon(@(x)obj(x,corr_fun,cl,p,q),x0,[],[],[],[],lb,ub,@(x)nonlcon(x,p,q),options );
        if p>0
            coeffs{1} = arma_coeffs(1:p);
        else
            coeffs{1} = [];
        end
        if q>0
            coeffs{2} = arma_coeffs(p+1:end);
        else
            coeffs{2} = [];
        end
        th = acf(coeffs{1},coeffs{2},lag_num,1);
        function err = obj(x,corr_fun,cl,p,q)
            try
                if p>0 && q>0
                    temp = acf(x(1:p),x(p+1:end),cl-1,1);
                    err = sum((corr_fun(1:cl)'-temp).^2);
                elseif p>0 && q==0
                    temp = acf(x,[],cl-1,1);
                    err = sum((corr_fun(1:cl)'-temp).^2);
                elseif q>0 && p==0
                    temp = acf([],x,cl-1,1);
                    err = sum((corr_fun(1:cl)'-temp).^2);
                end
            catch
                err = 2*cl;
            end
        end
        function [c,ceq] = nonlcon(yy,p,q)
            ceq = [];
            if p>0
                zz = yy(1:p);
                r = roots([1;-zz(:)]);
                c1 = nanmax(abs(r))-.99;
            else
                c1=-1;
            end
            if q>0
                zz = yy(p+1:end);
                r = roots([1;-zz(:)]);
                c2 = nanmax(abs(r))-.99;
            else
                c2=-1;
            end
            c=[c1;c2];
        end
    end
%% Computes the theoretical autocorrelations and long-run variance of an ARMA(p,q) process
    function [autocorr, sigma2_y] = acf(phi,theta,N,sigma2_e)
        % Copyright: Kevin Sheppard
        % kevin.sheppard@economics.ox.ac.uk
        % Revision: 3    Date: 1/1/2007
        if nargin<3 || nargin>4
            error('Requires 3 or 4 arguments.');
        end
        if length(N)>1 || any(N<0) || N~=floor(N)
            error('N must be a non-negative scalar integer.')
        end
        if nargin<4
            sigma2_e=1;
        elseif length(sigma2_e)>1 || any(sigma2_e<=0)
            error('sigma2_e must be a positive scalar');
        end
        p=length(phi);
        if min(size(phi))>1
            error('phi must be a vector');
        end
        if size(phi,1)~=p
            phi=phi';
        end
        q=length(theta);
        if min(size(theta))>1
            error('theta must be a vector');
        end
        if size(theta,1)~=q
            theta=theta';
        end
        if p>0
            [rho,stationary]=inverse_ar_roots(phi); %#ok<ASGLU>
            if ~stationary
                error('Autoregressive roots (phi) do not correspond to a stationary process');
            end
        end
        phi_original=phi;
        phi = [1 ; -phi; zeros(q-p,1) ];
        theta = [1 ; theta; zeros(p-q,1)];
        m=max(p,q)+1;
        phi_transformed  = zeros(m,m);
        T = toeplitz(1:m);
        for iii=1:m
            for jjj=1:m
                phi_transformed(iii,T(iii,jjj)) = phi_transformed(iii,T(iii,jjj))+phi(jjj);
            end
        end
        delta=tril(toeplitz(phi))^(-1)*theta;
        theta_transformed=flipud(tril(toeplitz(flipud(theta))));
        autocov=phi_transformed^(-1)*(theta_transformed*delta);
        if N+1<m
            autocov=autocov(1:N+1);
        elseif N+1>m
            autocov=[autocov;zeros(N-m+1,1)];
            if p>0
                for iii=m:N
                    autocov(iii+1)=phi_original'*autocov(iii:-1:(iii-p+1));
                end
            end
        end
        sigma2_y=sigma2_e*autocov(1);
        autocorr=autocov./autocov(1);
    end
%% Computes the inverted roots of the characteristic equation of an AR(P)
    function [rho, stationary]=inverse_ar_roots(phi)
        % Copyright: Kevin Sheppard
        % kevin.sheppard@economics.ox.ac.uk
        % Revision: 3    Date: 1/1/2007
        if size(phi,1)>=size(phi,2) && min(size(phi))==1
            phi=phi';
        else
            error('Phi should be a column vector.')
        end
        rho=roots([-fliplr(phi) 1]);
        stationary=min(abs(rho))>1;
    end
end
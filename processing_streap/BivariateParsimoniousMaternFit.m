function [ N1 , N2 , A , N12 , RHO12 ] = BivariateParsimoniousMaternFit( ser1 , ser2 , lagNum )
%    ( ser1 , ser2 , ser3 , lagNum )
%BivariateParsimoniousMaternFit estimates the parameters of the trivariate
% Matern parsimonious model.
% The estimation procedure is done by minimizing the square error of the
% estiamted autocovariances and cross covariances from the data. The
% estimation prcedure is a non linear least square constrained
% optimizatation.
%   Assumptions:
%   d=1 - time series
%   a1=a2=a12=a
%   beta12=1 - strong correlation between series
%   Model after:
%   Tilmann Gneiting, William Kleiber & Martin Schlather (2010) Matérn
%   Cross-Covariance Functions for Multivariate Random Fields, Journal of
%   the American Statistical Association, 105:491, 1167-1177,
%   DOI: 10.1198/jasa.2010.tm09420
%% Initilazing
covSer1 = xcov(ser1,ser1,lagNum,'unbiased');
covSer2 = xcov(ser2,ser2,lagNum,'unbiased');
% covSer3 = xcov(ser3,ser3,lagNum,'unbiased');
[crossCov12,dist] = xcov(ser1,ser2,lagNum,'unbiased');
% [crossCov13] = xcov(ser1,ser3,lagNum,'unbiased');
% [crossCov23] = xcov(ser2,ser3,lagNum,'unbiased');
var1 = nanvar(ser1(:));
var2 = nanvar(ser2(:));
% var3 = nanvar(ser3(:));
% Assumption - Fully symmetric cross covariance function
covSer1(dist<=0) = []; covSer2(dist<=0) = []; % covSer3(dist<=0) = []; 
crossCov12(dist<=0) = []; % crossCov13(dist<=0) = []; crossCov23(dist<=0) = []; 
dist(dist<=0) = [];
x0 = [.2 .2 0.01 .9 .3];
lb = [10^-5 10^-5 10^-5 10^-5 -1];
ub = [2 2 1 Inf 1];
A = [1 1 0 -1 0];
b = -10^-10;
options = optimset('Display','iter','Tolx',10^-15,'Tolfun',10^-15,'maxiter',10^10,'MaxFunEvals',10^10);
%% Fitting
res = fmincon(@(x) obj(x(1),x(2),x(3),x(4),x(5),covSer1,covSer2,crossCov12,dist,var1,var2),x0,A,b,[],[],lb,ub,@nonlcon,options);
N1 = res(1); N2 = res(2); A = res(3); N12 = res(4); RHO12 = res(5);
%% Nested functions
    function err = obj(n1,n2,a,n12,rho12,cov_ser1,cov_ser2,cross_cov12,dist,var1,var2)
        % func err = obj(n1,n2,n3,a,n12,n13,n23,rho12,rho13,rho23,cov_ser1,cov_ser2,cov_ser3,cross_cov12,cross_cov13,cross_cov23,dist,var1,var2,var3)
        cov1th = var1*(2^(1-n1))/gamma(n1)*(a*dist).^n1.*besselk(n1,a*dist);
        cov2th = var2*(2^(1-n2))/gamma(n2)*(a*dist).^n2.*besselk(n2,a*dist);
        % cov3th = var3*(2^(1-n3))/gamma(n3)*(a*dist).^n3.*besselk(n3,a*dist);
        cross_th12 = rho12.*(var1^.5)*(var2^.5)*(2^(1-n12))/gamma(n12)*(a*dist).^n12.*besselk(n12,a*dist);
        % cross_th13 = rho13.*(var1^.5)*(var3^.5)*(2^(1-n13))/gamma(n13)*(a*dist).^n13.*besselk(n13,a*dist);
        % cross_th23 = rho23.*(var2^.5)*(var3^.5)*(2^(1-n23))/gamma(n23)*(a*dist).^n23.*besselk(n23,a*dist);
        % err = sum((cov_ser1(:)-cov1th(:)).^2)+sum((cov_ser2(:)-cov2th(:)).^2)+sum((cov_ser3(:)-cov3th(:)).^2)+sum((cross_cov12(:)-cross_th12(:)).^2)+sum((cross_cov13(:)-cross_th13(:)).^2)+sum((cross_cov23(:)-cross_th23(:)).^2);
        err = sum((cov_ser1(:)-cov1th(:)).^2)+sum((cov_ser2(:)-cov2th(:)).^2)+sum((cross_cov12(:)-cross_th12(:)).^2);
        if ~isfinite(err)
            err=10^10;
        end
    end
    function [c,ceq] = nonlcon(yy)
        % Has been changed into Bivariate case.
        n1 = yy(1); n2 = yy(2); 
        n12 = yy(4);
        rho12 = yy(5);
        if n12==.5*(n1+n2)
            c = abs(rho12)-((gamma(n1+0.5)^.5)/gamma(n1)^.5)*((gamma(n2+0.5)^.5)/gamma(n2)^.5)*((gamma(.5*(n1+n2)))/gamma(.5*(n1+n2)+0.5));
            ceq = [];
        else
            ceq = [];
            c = [];
        end        
    end
end
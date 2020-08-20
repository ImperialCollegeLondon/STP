function A()
load('Birmingham/SpatialCov.mat','covariance','distance')
load('Birmingham/expoSpatialCorrelation.mat','expoSpatialCorr')
h = 1;
for mon = 1:12
    cb_obs = nanmean(covariance{mon},2);
    s = distance{1};
    cb_theo = getCB_THEO(-1./expoSpatialCorr(mon,1));
    subplot(3,4,mon)
    plot(cb_obs,'k.');hold on;plot(cb_theo,'r-');
    title(getMonthName(mon))
    drawnow
end
    function cb_theo = getCB_THEO(alpha_g)
        cb_theo = zeros(length(distance),1);
        for i = 1:length(s)
            s_this = s(i);
            cb_theo(i,1) = func_cbs(alpha_g,s_this,h);
        end
    end
    function cb_s = func_cbs(alpha_g,s_this,h)
        % all input/output is scalar
        expFunc = @(s,alpha)exp(-s/alpha);% exp(-s/alpha_g);
        cg_s = expFunc(s_this,alpha_g);
        cg_0 = expFunc(0,alpha_g);
        f1 = @(a,b)exp(-(a.^2+b.^2-2*cg_s*a.*b)./(2*(1-cg_s.^2)))./...
            (2*pi*sqrt(1-cg_s^2));
        r = integral(@(y)exp(-y.^2/2/cg_0),h,Inf,'RelTol',1e-8,'AbsTol',1e-13)...
            /sqrt(2*pi*cg_0);
        ps = integral2(@(a,b)f1(a,b),h,Inf,h,Inf,'RelTol',1e-8,'AbsTol',1e-13);
        cb_s = ps-r^2;
    end
end
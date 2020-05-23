function [ E0 ] = findE0( t )
%E0 Calculate correction factor E0
days = [31 28 31 30 31 30 31 31 30 31 30 31];    
nowYR=year(t); nowMO=month(t); nowDA=day(t);
if(nowMO==1)
    jDay = nowDA;
elseif(nowMO==2)
    jDay = days(1) + nowDA;
else
    jDay = sum(days(1:(nowMO-1))) + nowDA;
    if(mod(nowYR,4)==0)
        if(mod(nowYR,400)==0)
            jDay = jDay + 1;
        elseif(mod(nowYR,100)~=0)
            jDay = jDay + 1;
        end
    end
end
GA=2*pi*(jDay-1)/365;
E0=1.00011+0.034221*cos(GA)+0.00128*sin(GA)+0.000719*cos(2*GA)+0.000077*sin(2*GA);
end
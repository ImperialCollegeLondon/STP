function [ err ] = findGPparam( X , M , V )
[m,v] = gpstat(X(1),X(2),1);
err=sqrt((M-m).^2+(V-v).^2)/2;
end
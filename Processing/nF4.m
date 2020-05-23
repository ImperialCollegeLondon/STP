function [ x4 ] = nF4(fv,vk,G,AU,BU,cF,z0)
x4=(vk.*G)./sqrt((log(0.3.*(fv./abs(cF))./z0)-BU).^2+AU.^2)-fv;
end
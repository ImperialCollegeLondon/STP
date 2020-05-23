function [ x6 ] = nF6( fv,vk,G,cF,z0,MOL )
x6 = (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1./(1+0.882.*(sqrt((vk.*fv)./(abs(cF).*MOL))))))./z0)-0.5+2.55.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2+(4.5+1.765.*sqrt((vk.*fv)./(abs(cF).*MOL))).^2)-fv;
end
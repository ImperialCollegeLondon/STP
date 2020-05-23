function [ X1 ] = nF1( fv,vk,G,cF,z0,MOL )
X1 = (vk.*G)./sqrt((log((0.3.*(fv./abs(cF)).*(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))./z0)-10+(9.5./(1+0.027.*sqrt(-(vk.*fv)./(abs(cF).*MOL))))).^2+(4.5./(1+1.581.*sqrt(-(vk.*fv)./(abs(cF).*MOL)))).^2)-fv;
end
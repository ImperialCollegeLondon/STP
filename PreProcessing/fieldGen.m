function [ corrField , InvcorrField ] = fieldGen( dim , alpha )
%% fieldGen Generates a random correlated field, based on Athanasios Paschalis PhD thesis (2013).
%   Inputs:
%   dim - dimension of the domain [dim1 , dim2] [pixels]
%   alpha - correlation factor (Correlation parameter of the Gaussian field
%   divided by the pixel length)
%   Outputs:
%   field - random correlated field
%   InvcorrField - inverse normalized Gaussian field
%% Initilazing
[kx, ky] = meshgrid(-dim(1)/2:dim(1)/2-1,-dim(1)/2:dim(1)/2-1);
kx = kx/dim(1)*2*pi;
ky = ky/dim(2)*2*pi;
%% Exponential filter
FILTER=sqrt(2*pi*alpha^2./(sqrt(1^2+alpha^2*kx.^2+alpha^2*ky.^2).^3));
FILTER(~isfinite(FILTER)) = 0;
%% Generate normalized field N(miu,sigam)
corrField = real(ifft2(ifftshift(random('normal',0,1,[dim(1) dim(2)]).*FILTER),'symmetric'));
%% inverse the field U(N)
num_el = numel(corrField);
R = zeros(num_el,1);
[~, index] = sort(corrField(:));
R(index) = 1:num_el;
R = reshape(R, size(corrField));
InvcorrField = R/(num_el+1);
end
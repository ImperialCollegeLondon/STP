function [ qField ] = quantileFieldGen( dim , ARMA , Iter , alpha , t , simU , simV , dx , dt )
%quantileFieldGen A function that generates stochastic Gaussian field
%using FFT, taking into account the correlation ARMA proccess and the
%advection.
%   Inputs:

%% Initilazing
AR_coeffs = ARMA.ar;
ar_order = size(ARMA.ar,1);
MA_coeffs = ARMA.ma;
ma_order = size(ARMA.ma,1);
AR_std = sqrt(estimate_ARMA_noise_var(AR_coeffs,MA_coeffs,1));
U=0; V=0;
qField=cell(t(2)-t(1)+1,1);
%% Indices of repetition fourier coeffs - TO BE FIXED(only even numbers)
mid_ind = sub2ind([dim(1) dim(2)],dim(1)/2+1,dim(2)/2+1);
new_ind = 1:mid_ind;
rep_ind = mid_ind+1:(dim(1)*dim(2));
sam_ind = mid_ind-1:-1:mid_ind-length(rep_ind);
sam_l = length(new_ind);
buffer_re = randn(sam_l,ar_order);
buffer_im = randn(sam_l,ar_order);
buffer_re_ma = randn(sam_l,ma_order);
buffer_im_ma = randn(sam_l,ma_order);
%% run first iterations to stabilize buffer
temp_re = randn(dim(1),dim(2));
temp_im = randn(dim(1),dim(2));
for i = 1:Iter % Iterations
    err_noise_re = randn(sam_l,1)*AR_std;
    err_noise_im = randn(sam_l,1)*AR_std;
    new_re = buffer_re*AR_coeffs(end:-1:1)+buffer_re_ma*MA_coeffs(end:-1:1)+err_noise_re;
    new_im = buffer_im*AR_coeffs(end:-1:1)+buffer_im_ma*MA_coeffs(end:-1:1)+err_noise_im;
    temp_re(new_ind) = new_re;
    temp_im(new_ind) = new_im;
    temp_re(rep_ind) = temp_re(sam_ind);
    temp_im(rep_ind) = -temp_im(sam_ind); % complex conjugate
    buffer_re(:,1:end-1) = buffer_re(:,2:end);
    buffer_im(:,1:end-1) = buffer_im(:,2:end);
    buffer_re(:,end) = new_re;
    buffer_im(:,end) = new_im;
    buffer_re_ma(:,1:end-1) = buffer_re_ma(:,2:end);
    buffer_im_ma(:,1:end-1) = buffer_im_ma(:,2:end);
    buffer_re_ma(:,end) = err_noise_re;
    buffer_im_ma(:,end) = err_noise_im;
end
%% simulate images in the storm
[kx, ky] = meshgrid(-dim(1)/2:dim(1)/2-1,-dim(1)/2:dim(1)/2-1);
kx = kx/dim(1)*2*pi;
ky = ky/dim(2)*2*pi;
FILTER=sqrt(2*pi*alpha^2./(sqrt(1^2+alpha^2*kx.^2+alpha^2*ky.^2).^3));
FILTER(~isfinite(FILTER)) = 0;
for i = t(1)*2:t(2)*2
    %% Simulation of the new autoregressive map
    % new autoregressive random map
    err_noise_re = randn(sam_l,1)*AR_std;
    err_noise_im = randn(sam_l,1)*AR_std;
    new_re = buffer_re*AR_coeffs(end:-1:1)+buffer_re_ma*MA_coeffs(end:-1:1)+err_noise_re;
    new_im = buffer_im*AR_coeffs(end:-1:1)+buffer_im_ma*MA_coeffs(end:-1:1)+err_noise_im;
    temp_re(new_ind) = new_re;
    temp_im(new_ind) = new_im;
    temp_re(rep_ind) = temp_re(sam_ind);
    temp_im(rep_ind) = -temp_im(sam_ind); 
    noise_map = temp_re + 1i*temp_im; 
    if mod(i,2)==0 % Even numbers
        field = real(ifft2(ifftshift(noise_map.*FILTER),'symmetric'));
        u=round((simU(i/2)/1000)*60*dt/dx);
        v=round((simV(i/2)/1000)*60*dt/dx);
        field=circshift(field,[v+V u+U]);
        V=v+V;
        U=u+U;
        num_el = numel(field);
        R = zeros(num_el,1);
        [~, index] = sort(field(:));
        R(index) = 1:num_el;
        R = reshape(R, size(field));
        tmp = R/(num_el+1);
        qField{i/2} = tmp(1:dim(1)/2,1:dim(2)/2);
    end
    buffer_re(:,1:end-1) = buffer_re(:,2:end);
    buffer_im(:,1:end-1) = buffer_im(:,2:end);
    buffer_re(:,end) = new_re;
    buffer_im(:,end) = new_im;
    buffer_re_ma(:,1:end-1) = buffer_re_ma(:,2:end);
    buffer_im_ma(:,1:end-1) = buffer_im_ma(:,2:end);
    buffer_re_ma(:,end) = err_noise_re;
    buffer_im_ma(:,end) = err_noise_im;
end
qField=permute(cat(3,qField{:}),[1 2 3]);
end
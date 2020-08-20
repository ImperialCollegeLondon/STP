function [ qField ] = quantileFieldGen_ARMA( dim , decayA , Iter , alpha , t , simU , simV , dx , dt )
% quantileFieldGen A function that generates stochastic Gaussian field
% using FFT, taking into account the correlation ARMA proccess and the
% advection.
%    Inputs:

%% Initilazing

ar_order = 1;
% MA_coeffs = ARMA.ma;
% AR_coeffs = 0.80;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ARMA.ar;
% ma_order = size(ARMA.ma,1);
% AR_std = sqrt(estimate_ARMA_noise_var(AR_coeffs,MA_coeffs,1));

wavenumerBoundary = [linspace(0,500,500),Inf];
AR_coeffs = reshape(decayA*wavenumerBoundary(1:end-1)+0.999,[],1);
AR_coeffs(AR_coeffs<0.2) = 0.2;

% wavenumerBoundary =[0,Inf];
% AR_coeffs = [0.965];
% AR_coeffs(AR_coeffs<0.2) = 0.2;

% scaleBoundary = [Inf,  128,   64,   32,  16,  8,  4, -Inf];
% AR_coeffs =  [   .99805;   .9925;.9776;.9297;.796;.482;.482];
% scaleBoundary = [Inf,-Inf];%[Inf,  -Inf];
% AR_coeffs =  [0.998];%[0.9665];
AR_stds = [];
for leveli = 1:size(AR_coeffs,1)
    MA_coeffs = [];% currently, MA is not considered.
    AR_stds(leveli,1) = sqrt(estimate_ARMA_noise_var(...
        AR_coeffs(leveli,:),MA_coeffs,1));
end
[kx, ky] = meshgrid(-dim(1)/2:dim(1)/2-1,-dim(1)/2:dim(1)/2-1);
K = sqrt(kx.^2+ky.^2);

mid_ind = sub2ind([dim(1) dim(2)],dim(1)/2+1,dim(2)/2+1);
new_ind = 1:mid_ind;
wavelength = dim(1)*dx./K;%[Km]
scale = wavelength/2;%[Km]
[AR_coeffs,AR_stds] = assignARMAcoef(K(new_ind),wavenumerBoundary,AR_coeffs,AR_stds);

U=0; V=0;

qField=cell(t(2)-t(1)+1,1);
if ~isempty(MA_coeffs)
    % see original code 'quantileFieldGen'
else
    % AR model
    % added by yuting for calibrating AR(p) model
    %% Indices of repetition fourier coeffs - TO BE FIXED(only even numbers)
    rep_ind = mid_ind+1:(dim(1)*dim(2));
    sam_ind = mid_ind-1:-1:mid_ind-length(rep_ind);
    sam_l = length(new_ind);
    buffer_re = randn(sam_l,ar_order);
    buffer_im = randn(sam_l,ar_order);
    %% run first iteration to stabilize buffer
    temp_re = randn(dim(1),dim(2));
    temp_im = randn(dim(1),dim(2));
    for i = 1:Iter % Iterations
        err_noise_re = randn(sam_l,1).*AR_stds;
        err_noise_im = randn(sam_l,1).*AR_stds;
        new_re = buffer_re.*AR_coeffs(:,end:-1:1)...
            +err_noise_re;
        new_im = buffer_im.*AR_coeffs(:,end:-1:1)...
            +err_noise_im;
        temp_re(new_ind) = new_re;
        temp_im(new_ind) = new_im;
        temp_re(rep_ind) = temp_re(sam_ind);
        temp_im(rep_ind) = -temp_im(sam_ind); % complex conjugate
        buffer_re(:,1:end-1) = buffer_re(:,2:end);
        buffer_im(:,1:end-1) = buffer_im(:,2:end);
        buffer_re(:,end) = new_re;
        buffer_im(:,end) = new_im;
    end
    
    %% simulate images in the storm
    kx = kx/dim(1)*2*pi;
    ky = ky/dim(2)*2*pi;
    FILTER=sqrt(2*pi*alpha^2./(sqrt(1^2+alpha^2*kx.^2+alpha^2*ky.^2).^3));
    FILTER(~isfinite(FILTER)) = 0;
    for i = t(1)*2:t(2)*2
        %% Simulation of the new autoregressive map
        % new autoregressive random map
        err_noise_re = randn(sam_l,1).*AR_stds;
        err_noise_im = randn(sam_l,1).*AR_stds;
        new_re = buffer_re.*AR_coeffs(:,end:-1:1)+err_noise_re;
        new_im = buffer_im.*AR_coeffs(:,end:-1:1)+err_noise_im;
        
        temp_re(new_ind) = new_re;
        temp_im(new_ind) = new_im;
        temp_re(rep_ind) = temp_re(sam_ind);
        temp_im(rep_ind) = -temp_im(sam_ind);
        noise_map = temp_re + 1i*temp_im;
    
        if mod(i,2)==0 % Even numbers
            if isnan(simU(i/2))
                qField{i/2} = NaN(dim(1)/2,dim(2)/2);
            else
                
                field = real(ifft2(ifftshift(noise_map.*FILTER),'symmetric'));
                
                % Yt: In previous code, when speed < 12 Km/h (when dx = 2), rainfall will not move becasue of 'round'.
                % u=round((simU(i/2)/1000)*60*dt/dx);% deleted by YT
                % v=round((simV(i/2)/1000)*60*dt/dx);% deleted by YT
                % field=circshift(field,[v+V u+U]);% deleted by YT
                u=(simU(i/2)/1000)*60*dt/dx;% addeded by YT
                v=(simV(i/2)/1000)*60*dt/dx;% addeded by YT
                field=circshift(field,round([v+V u+U])); % addeded by YT
                V=v+V;
                U=u+U;
                num_el = numel(field);
                R = zeros(num_el,1);
                [~, index] = sort(field(:));
                R(index) = 1:num_el;
                R = reshape(R, size(field));
                tmp = R/(num_el+1);
                % qField{i/2} = tmp(1:dim(1)/2,1:dim(2)/2);% deleted by YT
                centInd = round(dim(1)/4):round(dim(1)/4)+dim(1)/2-1;% addeded by YT
                qField{i/2} = tmp(centInd,centInd);% addeded by YT
            end
        end
        buffer_re(:,1:end-1) = buffer_re(:,2:end);
        buffer_im(:,1:end-1) = buffer_im(:,2:end);
        buffer_re(:,end) = new_re;
        buffer_im(:,end) = new_im;

    end
end
qField=permute(cat(3,qField{:}),[3 1 2]);%%%%%%[1 2 3]);%%%% yt:Order changed to [T,loc1,loc2]
end

function [AR_coeff_, AR_std_] = assignARMAcoef(wavenumbers,wnBoundary,AR_coeffs,AR_stds);
[AR_coeff_, AR_std_] = deal(NaN(length(wavenumbers),size(AR_coeffs,2)));

for b = 1:length(wnBoundary)-1
    index = wavenumbers<=wnBoundary(b+1) & wavenumbers>=wnBoundary(b);
    AR_coeff_(index,:) = AR_coeffs(b,:);
    AR_std_(index,:) = AR_stds(b,:);
end

end





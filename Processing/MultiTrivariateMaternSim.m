function [ s1 , s2 , s3 ] = MultiTrivariateMaternSim( N1 , N2 , N3 , A , N12 , N13 , N23 , RHO12 , RHO13 , RHO23 , N , K )
%TrivariateMaternSim: Simulates a gaussian vector process with specific
%covariance structrure beloging to the Matern class of Covariances
% The methodology uses the fft approach.
% Ruan, F., and D. McLaughlin (1998), An efficient multivariate random
% field generator using the fast Fourier transform, Advances in Water
% Resources, 21(5), 385-399, doi:10.1016/S0309-1708(96)00064-4.
% ------------------------------------------------------------- %
% Assumtion both processes are zero mean unit variance gaussian %
% ------------------------------------------------------------- %
%% Initilazing
if mod(N,2)
    cut = 1;
    N=N+1;
else
    cut = 0;
end
%% Wavenumbers of the Fourier transform to be inverted
wn = (-N/2:1:N/2-1)/N*2*pi;
sp11 = (gamma(N1+0.5)*A^(2*N1)/gamma(N1)/pi^(0.5))./(A^2+abs(wn).^2).^(N1+0.5);
sp22 = (gamma(N2+0.5)*A^(2*N2)/gamma(N2)/pi^(0.5))./(A^2+abs(wn).^2).^(N2+0.5);
sp33 = (gamma(N3+0.5)*A^(2*N3)/gamma(N3)/pi^(0.5))./(A^2+abs(wn).^2).^(N3+0.5);
sp12 = RHO12*(gamma(N12+0.5)*A^(2*N12)/gamma(N12)/pi^(0.5))./(A^2+abs(wn).^2).^(N12+0.5);
sp13 = RHO13*(gamma(N13+0.5)*A^(2*N13)/gamma(N13)/pi^(0.5))./(A^2+abs(wn).^2).^(N13+0.5);
sp23 = RHO23*(gamma(N23+0.5)*A^(2*N23)/gamma(N23)/pi^(0.5))./(A^2+abs(wn).^2).^(N23+0.5);
for i = 1:N/2+1
    try
        temp = [sp11(i) sp12(i) sp13(i); sp12(i) sp22(i) sp23(i); sp13(i) sp23(i) sp33(i)];
        if i<N/2+1 && i~=1
            series{i,1} = mvnrnd([0 0 0],temp,K)'+1i*mvnrnd([0 0 0],temp,K)';
        else
            series{i,1} = mvnrnd([0 0 0],temp,K)';
        end
    catch
        checkEig=eig(temp+temp')/2;
        while checkEig(1)<=1*10^-10 % SIGMA not a symmetric positive semi-definite matrix
            temp=temp+min(min([temp(1,1) temp(2,2) temp(3,3)]))*.01*eye(3);
            checkEig=eig(temp+temp')/2;
        end
        if i<N/2+1 && i~=1
            series{i,1} = mvnrnd([0 0 0],temp,K)'+1i*mvnrnd([0 0 0],temp,K)';
        else
            series{i,1} = mvnrnd([0 0 0],temp,K)';
        end
    end
end
for i=1:size(series,1)
    s1(i,:)=series{i}(1,:);
    s2(i,:)=series{i}(2,:);
    s3(i,:)=series{i}(3,:);
end
s1(N,K)=0; s2(N,K)=0; s3(N,K)=0;
s1(N/2+2:end,:) = s1(N/2+1:-1:3,:);
s2(N/2+2:end,:) = s2(N/2+1:-1:3,:);
s3(N/2+2:end,:) = s3(N/2+1:-1:3,:);
s1 = real(ifft(ifftshift(s1,1),[],1));
s2 = real(ifft(ifftshift(s2,1),[],1));
s3 = real(ifft(ifftshift(s3,1),[],1));
s1 = bsxfun(@rdivide,s1,nanstd(s1,0,1));
s2 = bsxfun(@rdivide,s2,nanstd(s2,0,1));
s3 = bsxfun(@rdivide,s3,nanstd(s3,0,1));
if cut
    s1 = s1(1:end-1,:);
    s2 = s2(1:end-1,:);
    s3 = s3(1:end-1,:);
end
end
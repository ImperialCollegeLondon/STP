function [ gpFit ] = gpFitR( data , theta )
%gpFitR Fit GP distribution to data
%   Input:
%   data - observed or simulated daily rainfall [mm] series for each pixel in
%   the domain. data{i}{m} - i is the pixel order and m is the month.
%   Output:
%   gpFit - cell array storing the k and sigma parameters for each pixel
%% Compute
for m=1:12
    gpFit{1}{m}=zeros(sqrt(size(data,2)),sqrt(size(data,2)));
    gpFit{2}{m}=zeros(sqrt(size(data,2)),sqrt(size(data,2)));
end
for i=1:size(data,2)
    for m=1:12
        tmp=data{i}{m};
        tmp(tmp<theta(i))=[];
        aa=gpfit(tmp-theta(i));
        gpFit{1}{m}(i)=aa(1);
        gpFit{2}{m}(i)=aa(2);
    end
end
end
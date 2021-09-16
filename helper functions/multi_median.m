function out = multi_median(volumes,radius,passes)
% function multi_median()
%
% reproduces the multi_median function from Dipy
%
% default values, as recommended in dipy:
% passes = 4
% radius = 4

switch nargin
    case 1
        radius = 4;
        passes = 4;
    case 2
        passes = 4;
end

% if volumes is 4D, sums along 4th direction
if length(size(volumes)) == 4
    out = sum(volumes,4);
else
    out = volumes;
end

% create kernel for median filter
kernel_size = radius*2+1;
kern = [kernel_size kernel_size kernel_size];

for i = 1:passes
    out = medfilt3(out,kern);
end

end
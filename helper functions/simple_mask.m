function out = simple_mask(volumes,radius,passes)
% function simple_mask()
%
% returns a brain mask 
% based on the median_otsu method from Dipy
%

% check inputs for multi_median
switch nargin
    case 1
        radius = 2;
        passes = 2;
    case 2
        passes = 2;
end

% sum all volumes along 4th direction, pass it with a median filter and normalize
Vsum = sum(volumes,4);
Vsum = multi_median(Vsum,radius,passes);
Vnorm = Vsum/max(Vsum(:));

% compute gradient & use Otzu method to create a binary image
out = zeros(size(Vnorm));
for i = 1 : size(Vnorm,3)
    th = graythresh(Vnorm(:,:,i));
    out(:,:,i) = imbinarize(Vnorm(:,:,i),th);
end

% fill holes slice by slice
for i = 1:size(out,3)
    tmp = out(:,:,i);
    out(:,:,i) = imfill(tmp, 'holes');
end

end

function [data,vox_add,R,Q] = reshape_data_for_multivoxfit(data,nvox)
% function reshape_data_for_multivox()
%
% reshape the input data for multivoxel SDP computation
% data must be in the shape [#meas, number_voxels]

if nargin < 2
    nvox = 50;
end

% get quotient and reminder
siz = size(data);
R = rem(siz(2),nvox);
Q = fix(siz(2)/nvox);

% if necessary, complete the data by adding extra voxels at the end
vox_add = 0;
if R > 0
vox_add = nvox - R;
data = cat(2,data,data(:,1:vox_add));
end

% reshape to do multivox fit
data = reshape(data, siz(1) * nvox,[]);

end
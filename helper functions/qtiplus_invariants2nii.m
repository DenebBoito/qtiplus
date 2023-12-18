function qtiplus_invariants2nii(dps, nii_h)
% function qtiplus_invariants2nii(dps)
%
% writes the invariants in the dps structure into nifti file and attaches
% the specified nifti header

invs = fieldnames(dps);
n = length(invs);
nii_h_mod = nii_h;
nii_h_mod.ImageSize = horzcat(nii_h.ImageSize, 3);
nii_h_mod.PixelDimensions = horzcat(nii_h.PixelDimensions, 1);


for i = 1 : n
    fieldvalue = getfield(dps, invs{i});
    
    % adjust for 4D images
    if numel(size(fieldvalue)) > 3
        niftiwrite(fieldvalue, invs{i}, nii_h_mod);
        %continue
    else
        niftiwrite(fieldvalue, invs{i}, nii_h);
    end
end
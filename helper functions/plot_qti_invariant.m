function plot_qti_invariant(s,inv_name,disp_range)
% function to plot one invariant for different fits/protocols
%

n_fits = max(size(s));

% by default, we take the middle slice to be shown
tmp_inv = getfield(s{1},inv_name);
slice = round(size(tmp_inv,3)/2);

% open a figure
figure()

% for each fit (invariant structure), we plot the QTI main invariants
tmp = [];
for i = 1:n_fits
    % for each fit (invariant structure) extract the invariant
    i_tmp = getfield(s{i},inv_name);
    i_tmp = (fliplr(i_tmp(:,:,slice)))';
    tmp = cat(3,tmp,i_tmp);
end

if license('test', 'image_toolbox')
    
    montage(tmp, 'size', [1, n_fits], 'DisplayRange', disp_range)
    
else % use subplot
    
    for i = 1:n_fits
        subplot(1,n_fits,i)
        imshow(tmp(:,:,i), disp_range)
        colormap gray
        axis square
        axis off
    end
    
end
end
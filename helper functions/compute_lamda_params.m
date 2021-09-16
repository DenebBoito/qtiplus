function [FA, AD, RD, ax_dir, FA_col] = compute_lamda_params(d,siz)
% function compute_lambda_params()
%
% helper function to compute parameters from a second order tensor's
% eigenvalues
%

FA = zeros(size(d,1),1);
AD = zeros(size(d,1),1);
RD = zeros(size(d,1),1);
ax_dir = zeros(size(d,1),3);
FA_col = zeros(size(d,1),3);


for i = 1:size(d,1)
   
    % avoid off mask voxels
    if sum(d(i,:)) ~= 0
    
    [evec,eval] = eigs(convert_1x6_to_3x3(d(i,:)));
    
    [l,indx] = sort(diag(eval), 'descend');
    
    FA(i) = sqrt(0.5) *...
                      sqrt( (l(1) - l(2))^2 + (l(2) - l(3))^2 + (l(1) - l(3))^2) / ...
                      max(sqrt( l' * l) ,eps);
    
    
    AD(i) =  l(1);
    RD(i) = (l(2) + l(3)) * 0.5;
    ax_dir(i,:) = evec(:,indx(1));
    FA_col(i,:) = FA(i) * abs(ax_dir(i,:));
    
    
    end
    
end

% reshape
FA = reshape(FA, siz(1), siz(2), siz(3));
AD = reshape(AD, siz(1), siz(2), siz(3));
RD = reshape(RD, siz(1), siz(2), siz(3));
ax_dir = reshape(ax_dir, siz(1), siz(2), siz(3),3);
FA_col = reshape(FA_col, siz(1), siz(2), siz(3),3);


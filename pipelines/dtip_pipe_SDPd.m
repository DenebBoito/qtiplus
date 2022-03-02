function model = dtip_pipe_SDPd(data,btensors,mask,nvox,ind,parallel,cvxsolver)
% function model = qtip_pipe_SDPdc()
%

% reshape data and retain voxels within the mask
siz = size(data);
data = reshape(data, prod(siz(1:3)), siz(4))';
mask = reshape(mask, prod(siz(1:3)),1)';
si = find((mask>0).*(~all(data==0,1))); 
data = double(data(:,si));
% data(data <= 0) = eps;

% reshape data for multivox fit
[data_fit,vox_add,~,~] = reshape_data_for_multivoxfit(data,nvox);
[~,y_fit] = size(data_fit);
m_fit = zeros(7*nvox,y_fit);

% do fit
if parallel
    parfor k = 1:y_fit
        
        m_fit(:,k) = dtip_data2fit_SDPd(data_fit(:,k), btensors, ...
                                       'nvox', nvox,...
                                       'ind', ind,...
                                       'outflag', 1,...
                                       'cvxsolver',cvxsolver);
        
    end
else
    for k = 1:y_fit
        
        m_fit(:,k) = dtip_data2fit_SDPd(data_fit(:,k), btensors, ...
                                       'nvox', nvox,...
                                       'ind', ind,...
                                       'outflag', 1,...
                                       'cvxsolver',cvxsolver);
        
    end
end

% revert shape & compute outputs
m_fit = reshape(m_fit,7,[]);
if vox_add > 0
    m_fit(:,end - vox_add + 1:end) = [];
end
model = zeros(7, prod(siz(1:3)));
model(:,si) = m_fit;
model = reshape(model', siz(1), siz(2), siz(3), 7);

end
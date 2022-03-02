function model = dtip_pipe_NLLSd(model,data,btensors,mask,ind,parallel)
% function model = dtip_pipe_NLLSd()
%

% reshape data and retain voxels within the mask
siz = size(data);
sizm = size(model);
data = reshape(data, prod(siz(1:3)), siz(4))';
model = reshape(model, prod(sizm(1:3)), sizm(4))';
mask = reshape(mask, prod(siz(1:3)),1)';
si = find((mask>0).*(~all(data==0,1))); 
data = double(data(:,si));
model = double(model(:,si));
%data(data <= 0) = eps;
[~,y_fit] = size(data);
m_fit = zeros(7,y_fit);

% do fit
if parallel
    parfor k = 1:y_fit
        
        m_fit(:,k) = qtip_data2fit_NLLSdc(model(:,k), data(:,k), btensors, ind);
        
    end
else
    for k = 1:y_fit
        
        m_fit(:,k) = qtip_data2fit_NLLSdc(model(:,k), data(:,k), btensors, ind);
        
    end
end

% revert shape & compute outputs
m_fit = reshape(m_fit,7,[]);
model = zeros(7, prod(siz(1:3)));
model(:,si) = m_fit;
model = reshape(model', siz(1), siz(2), siz(3), 28);

end
function model = qtip_pipe_SDPdcm(model,data,btensors,mask,nvox,mcheckflag,ind,parallel,cvxsolver)
% function model = qtip_pipe_SDPdcm()
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
% data(data <= 0) = eps;

% reshape data for multivox fit
[data_fit,vox_add,~,~] = reshape_data_for_multivoxfit(data,nvox);
model = reshape_data_for_multivoxfit(model,nvox);
[~,y_fit] = size(data_fit);
m_fit = zeros(28*nvox,y_fit);


% do fit, check whether mcheck required
if parallel
    if mcheckflag
        
        mcheck = zeros(1,y_fit);
        
        parfor k = 1:y_fit
            
            % do m-check
            mcheck(k) = qtip_mcheck(reshape(model(:,k),28,[]));
            if mcheck(k)
                m_fit(:,k) = qtip_data2fit_SDPdcm(model(:,k), data_fit(:,k), btensors, ...
                                                 'nvox', nvox, ...
                                                 'ind', ind, ...
                                                 'outflag', 1, ...
                                                 'cvxsolver', cvxsolver);
            else
                m_fit(:,k) = model(:,k);
            end
        end
        
    else
        
        parfor k = 1:y_fit
            
            m_fit(:,k) = qtip_data2fit_SDPdcm(model(:,k), data_fit(:,k), btensors, ...
                                                 'nvox', nvox, ...
                                                 'ind', ind, ...
                                                 'outflag', 1, ...
                                                 'cvxsolver', cvxsolver);
            
        end
    end
else
    if mcheckflag
        
        mcheck = zeros(1,y_fit);
        
        for k = 1:y_fit
            
            % do m-check
            mcheck(k) = qtip_mcheck(reshape(model(:,k),28,[]));
            if mcheck(k)
                m_fit(:,k) = qtip_data2fit_SDPdcm(model(:,k), data_fit(:,k), btensors, ...
                                                 'nvox', nvox, ...
                                                 'ind', ind, ...
                                                 'outflag', 1, ...
                                                 'cvxsolver', cvxsolver);
            else
                m_fit(:,k) = model(:,k);
            end
        end
        
    else
        
        for k = 1:y_fit
            
            m_fit(:,k) = qtip_data2fit_SDPdcm(model(:,k), data_fit(:,k), btensors, ...
                                                 'nvox', nvox, ...
                                                 'ind', ind, ...
                                                 'outflag', 1, ...
                                                 'cvxsolver', cvxsolver);
            
        end
    end
end

% revert shape & compute outputs
m_fit = reshape(m_fit,28,[]);
if vox_add > 0
    m_fit(:,end - vox_add + 1:end) = [];
end
model = zeros(28, prod(siz(1:3)));
model(:,si) = m_fit;
model = reshape(model', siz(1), siz(2), siz(3), 28);

end
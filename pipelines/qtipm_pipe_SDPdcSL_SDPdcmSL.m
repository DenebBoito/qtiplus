function [model, model_SDPdcSL] = qtipm_pipe_SDPdcSL_SDPdcmSL(data,btensors,D0,mask,nvox,mcheckflag,ind,parallel,cvxsolver)
% function qtip_pipe_SDPdc_SDPdcm()
%
%
% utility to not have to reshape the data between fits


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
m_fit_SDPdcSL = zeros(28*nvox,y_fit);
m_fit = zeros(28*nvox,y_fit);
if mcheckflag
    mcheckp = zeros(1,y_fit);
    mcheckm = zeros(1,y_fit);
end

% do fit
if parallel
    parfor k = 1:y_fit
        
        % SDPdc
        m_fit_SDPdcSL(:,k) = qtipm_data2fit_SDPdcSL(data_fit(:,k), btensors, ...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
                                               'cvxsolver', cvxsolver);
        
        % SDPdcm
        if mcheckflag
            mcheckp(k) = qtip_mcheck(reshape(m_fit_SDPdcSL(:,k),28,[]));
            mcheckm(k) = qtipm_mcheck(reshape(m_fit_SDPdcSL(:,k),28,[]), D0);
            if mcheckp(k) || mcheckm(k)
                m_fit(:,k) = qtipm_data2fit_SDPdcmSL(m_fit_SDPdcSL(:,k), data_fit(:,k), btensors,...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
                                               'cvxsolver', cvxsolver);
            else
                m_fit(:,k) = m_fit_SDPdcSL(:,k);
            end
            
        else
            m_fit(:,k) = qtipm_data2fit_SDPdcm(m_fit_SDPdcSL(:,k),data_fit(:,k), btensors, ...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
                                               'cvxsolver', cvxsolver);
        end
    end
else
    for k = 1:y_fit
        
        % SDPdc
        m_fit_SDPdcSL(:,k) = qtipm_data2fit_SDPdcSL(data_fit(:,k), btensors, ...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
                                               'cvxsolver', cvxsolver);
        
        % SDPdcm
        if mcheckflag
            mcheckp(k) = qtip_mcheck(reshape(m_fit_SDPdcSL(:,k),28,[]));
            mcheckm(k) = qtipm_mcheck(reshape(m_fit_SDPdcSL(:,k),28,[]), D0);
            if mcheckp(k) || mcheckm(k)
                m_fit(:,k) = qtipm_data2fit_SDPdcmSL(m_fit_SDPdcSL(:,k), data_fit(:,k), btensors, ...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
                                               'cvxsolver', cvxsolver);
            else
                m_fit(:,k) = m_fit_SDPdcSL(:,k);
            end
            
        else
            m_fit(:,k) = qtipm_data2fit_SDPdcmSL(m_fit_SDPdcSL(:,k),data_fit(:,k), btensors, ...
                                               'nvox', nvox,...
                                               'ind', ind, ...
                                               'D0', D0, ...
                                               'outflag', 1,...
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

% revert also initial step
m_fit_init = reshape(m_fit_SDPdcSL,28,[]);
if vox_add > 0
    m_fit_init(:,end - vox_add + 1:end) = [];
end
model_SDPdcSL = zeros(28, prod(siz(1:3)));
model_SDPdcSL(:,si) = m_fit_init;
model_SDPdcSL = reshape(model_SDPdcSL', siz(1), siz(2), siz(3), 28);
end
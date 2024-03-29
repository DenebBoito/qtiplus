function model = dtipm_data2fit_SDPdSL(signal,bten,varargin)
% function model = dtipm_data2fit_SDPdSL()
%
% out_flag says which format to use for m : 
%                                           out_flag = 1 -> [28*nvox,1]
%                                           out_flag = 0 -> [28, nvox]

nmeas = size(bten,1);
% retrieve optional inputs
p = inputParser;
p.addParameter('nvox', 1)
p.addParameter('ind', ones(nmeas,1)>0)
p.addParameter('outflag',1)
p.addParameter('cvxsolver', 'sdpt3')
p.addParameter('D0', 3.1)
p.parse(varargin{:})
nvox        = p.Results.nvox;
ind         = p.Results.ind;
outflag     = p.Results.outflag;
cvxsolver   = p.Results.cvxsolver;
D0          = p.Results.D0;



% check signal size, determine which format has been used
flag = 0;
if size(signal,1) == nmeas * nvox
    % rearrange data
    signal = reshape(signal, [nmeas nvox]);
    flag = 1;
else
    % signal should be already [nmeas,nvox]  
    if size(signal,1) ~= nmeas
        signal = signal';
    end
end
nmeas = min([sum(ind), size(bten,1)]);


%% Setup variables 
% remove non ind signals
indx = repmat(ind, [1, nvox]);
signal = reshape(signal(indx),[sum(ind), nvox]);

% regressors, create a persistent variable 
b2 = bten(ind,:) * 1e-9 ;
regressors = get_regressors(nvox,b2,nmeas);

% normalize signals
signal_scales = repmat(max(signal,[],1),[nmeas,1]);
signal = signal ./signal_scales;
signallog = real(log(signal)); 
signallog(~isfinite(signallog)) = 0;

% least squares variables
b = signal .* signallog;
s_reshaped = repmat(reshape(signal,[],1,nvox), [1 7 1]);
A = s_reshaped .* regressors;

% compute P and c
% check if pagemtimes available
if verLessThan('Matlab','9.9')
    P = arrayfun(@(i)sqrtm(A(:,:,i)' * A(:,:,i)), 1:nvox, 'UniformOutput', false);
    P = cat(3, P{:});
    c = arrayfun(@(i)-2*A(:,:,i)'*b(:,i), 1:nvox, 'UniformOutput', false);
    c = cat(2, c{:});
    
else % pagemtimes should be available
    Q = pagemtimes(A,'transpose',A, 'none');
    c = -2 * pagemtimes(A,'transpose',permute(b,[1 3 2]), 'none');
    c = squeeze(c);
    P = zeros(7,7,nvox);
    for i = 1:nvox
        P(:,:,i) = sqrtm(Q(:,:,i));
    end
end

% variable for speed limit
D0eye3 = D0 * eye(3);
 
%% cvx 
cvx_begin sdp quiet

        % select solver
        if strcmpi(cvxsolver,'mosek')
            cvx_solver mosek
        else
            cvx_solver sdpt3
        end

        % define unknowns
        variable t(nvox)
        variable x(7,nvox) 

        % define objective function
        minimize(sum(t))

        % constraints
        for i = 1:nvox

            % problem
            [eye(7)                             P(:,:,i) * x(:,i); ...
             (P(:,:,i) * x(:,i))'               t(i) - c(:,i)'* x(:,i)]  >= 0;

            % (d) constraint
            convert_1x6_to_3x3(x(2:7,i)) >= 0;

            % (d-) speed limit
            D0eye3 - convert_1x6_to_3x3(x(2:7,i)) >= 0;

        end

cvx_end

%% output
% need to check how data were input
if flag && outflag
    x(1,:) = exp(x(1,:));
    model = x .* [signal_scales(1,:); 1e-9 * ones(6,nvox)] ;
    model = model(:);
else % switch to [7,nvox] format
    x(1,:) = exp(x(1,:));
    model = x .* [signal_scales(1,:); 1e-9 * ones(6,nvox)] ;
end

end

function regressors = get_regressors(nvox,b2,nmeas)

persistent reg;

if  ~isempty(reg) && (numel(reg.b2(:)) == numel(b2(:))) && reg.nvox == nvox && all(reg.b2(:) == b2(:))
    regressors = reg.regressors;
    return
end

b0 = ones(nmeas, 1);
regressors = repmat([b0 -b2], [1 1 nvox]);
reg.regressors = regressors;
reg.b2 = b2;
reg.nvox = nvox;

end
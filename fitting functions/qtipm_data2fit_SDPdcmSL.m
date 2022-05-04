function model = qtipm_data2fit_SDPdcmSL(model,signal,bten,varargin)
% function qtipm_data2fit_SDPdcmSL()
%
% flag is for define output size:
%                                   flag = 1 -> m is [nvox*28,1] (default);
%                                   flag = 0 -> m is [28,nvox]

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
if size(signal,1) == nmeas * nvox
    % rearrange data & m
    signal = reshape(signal, [nmeas nvox]);
    model = reshape(model, 28,[]);
else
    % signal & m should be already [nmeas,nvox]  
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
s_reshaped = repmat(reshape(signal,[],1,nvox), [1 28 1]);
A = s_reshaped .* regressors;

d = model(1:7,:);
d(1,:) = log(d(1,:)./signal_scales(1,:)); % if scaling signal, remember to scale this as well
d(2:7,:) = d(2:7,:) * 1e9; % change units

% compute Q,P,c2,cy
% check if pagemtimes available
if verLessThan('Matlab','9.9')
    Q = arrayfun(@(i)(A(:,:,i)' * A(:,:,i)), 1:nvox, 'UniformOutput', false);
    Q = cat(3,Q{:});
    Q12 = Q(1:7,8:28,:);
    Q22 = Q(8:28, 8:28,:);

    % compute P and c
    P = arrayfun(@(i)sqrtm(Q22(:,:,i)), 1:nvox, 'UniformOutput', false);
    P = cat(3,P{:});
    c = arrayfun(@(i)-2*A(:,:,i)'*b(:,i), 1:nvox, 'UniformOutput', false);
    c = cat(2,c{:});

    % compute c2,cy
    d = model(1:7,:);
    d(1,:) = log(d(1,:)./signal_scales(1,:)); % if scaling signal, remember to scale this as well
    d(2:7,:) = d(2:7,:) * 1e9;

    c2 = c(8:28,:);
    cy = arrayfun(@(i) (2 * d(:,i)' * Q12(:,:,i) + c2(:,i)')', 1:nvox, 'UniformOutput', false); % [21*1]
    cy = cat(2,cy{:});
    
else % pagemtimes should be available
    Q = pagemtimes(A,'transpose',A, 'none');
    c = -2 * pagemtimes(A,'transpose',permute(b,[1 3 2]), 'none');
    c = squeeze(c);
    Q12 = Q(1:7,8:28,:);
    c2 = c(8:28,:);
    dtimesQ12 = pagemtimes(permute(d,[1 3 2]),'transpose',Q12,'none');
    cy = 2 * squeeze(dtimesQ12) + c2;
    P = zeros(21,21,nvox);
    for i = 1:nvox
        P(:,:,i) = sqrtm(Q(8:28,8:28,i));
    end
end

% define variables for speed limit
mask_c2_SL = 3/4 * D0^2 * eye(6);
D024 = D0^2 / 4;
 
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
        variable y(21,nvox) 
        variables lambda(9,nvox) lambda_sl(9,nvox) lambda_gp(9,nvox) lambda_gm(9,nvox)

        % define objective function
        minimize(sum(t))

        % constraints
        for i = 1:nvox

            % problem
            [eye(21)                             P(:,:,i) * y(:,i); ...
            (P(:,:,i) * y(:,i))'               t(i) - cy(:,i)' * y(:,i)]  >= 0;

            % C & M constraints
            convert_1x21_to_6x6(y(:,i)) >= 0;
            get_M_9x9_from_d_c(convert_1x6_to_3x3(d(2:7,i)),convert_1x21_to_3x3x3x3(y(:,i)),lambda(:,i)) >= 0

            % (c-) diagonal components speed  limit
            D024 - y(1,i) >= 0
            D024 - y(2,i) >= 0
            D024 - y(3,i) >= 0

            % (c-) off-diagonal components speed limit, beware scaling due to voigt
            % format
            D024 + y(4,i)/sqrt(2) >= 0
            D024 + y(5,i)/sqrt(2) >= 0
            D024 + y(6,i)/sqrt(2) >= 0

            D024 - y(4,i)/sqrt(2) >= 0
            D024 - y(5,i)/sqrt(2) >= 0
            D024 - y(6,i)/sqrt(2) >= 0

            % (c-) eigenvalues speed limit
            mask_c2_SL - convert_1x21_to_6x6(y(:,i)) >= 0;

            % Gamma plus
<<<<<<< HEAD
            % qtipm_get_gammap(tm_1x21_to_3x3x3x3_cvx(y(:,i)),lambda_gp(:,i), D0) >= 0
=======
            % qtipm_get_gammap(convert_1x21_to_3x3x3x3(y(:,i)),lambda_gp(:,i), D0) >= 0
>>>>>>> dtiplus

            % Gamma minus
            qtipm_get_gammam(convert_1x21_to_3x3x3x3(y(:,i)),lambda_gm(:,i), D0) >= 0

            % (m-) speed limit
            get_M_9x9_from_d_c_SL(convert_1x6_to_3x3(d(2:7,i)),convert_1x21_to_3x3x3x3(y(:,i)),lambda_sl(:,i), D0) >= 0

        end

cvx_end

%% output
% check how data were input
if outflag
    model(8:28,:) = y * 1e-18;
    model = model(:);
else
    model(8:28,:) = y * 1e-18;
end

end

function regressors = get_regressors(nvox,b2,nmeas)

persistent reg;

if  ~isempty(reg) && (numel(reg.b2(:)) == numel(b2(:))) && reg.nvox == nvox && all(reg.b2(:) == b2(:))
    regressors = reg.regressors;
    return
end

b0 = ones(nmeas, 1);
b4 = convert_1x6_to_1x21(b2);
regressors = repmat([b0 -b2 1/2 * b4], [1 1 nvox]);
reg.regressors = regressors;
reg.b2 = b2;
reg.nvox = nvox;

end
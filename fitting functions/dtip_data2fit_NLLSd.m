function model = dtip_data2fit_NLLSd(model,signal,bten,ind)
% function m = dtip_data2fit_NLLSd(model,signal,bten,ind)
%
% fits the DTI parameters using NLLSd

if (nargin < 4), ind = ones(size(bten,1),1)>0; end

signal = signal(ind);
bten = bten(ind,:);
model = double(model);

% convert between fit and SI units
SI_to_fit = 1./[max(signal) [1 1 1 1 1 1] * 1e-9];
fit_to_SI = [max(signal) [1 1 1 1 1 1] * 1e-9];

% check if column or row
if size(model,1) > 1
    model = model';
end
model = model .* SI_to_fit;

% get the cholesky factorization of D
Dchol = cholesky_factorization(convert_1x6_to_3x3(model(2:7)));

% helper functions
function model = tmp_to_model(tmp)
    model(1) = tmp(1); 
    Dc = triu(convert_1x6_to_3x3(tmp(2:7)));
    model(2:7) = convert_3x3_to_1x6(Dc' * Dc);
    model = model .* fit_to_SI;
end

function s = dti_fit2data(t,varargin)
    model = tmp_to_model(t);
    s = model(1) * exp(-bten * model(2:7)');
end

% initial guess, fit opts, and fit attempt
init_guess = [model(1) convert_3x3_to_1x6(Dchol)];

opts = optimoptions(@lsqcurvefit,...
    'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',4e3,...
    'MaxIterations',2e3,...          
    'FunctionTolerance', 1e-8,...
    'StepTolerance', 1e-8,...
    'display', 'off');

try
    out_guess = lsqcurvefit(@dti_fit2data, init_guess, [], signal, [], [], opts);
catch
    out_guess = init_guess;
end

% check residuals: if worse than input, use input
if sum((signal -  dti_fit2data(out_guess, bten)).^2) > sum((signal -  dti_fit2data(init_guess, bten)).^2)
    out_guess = init_guess; 
end

% return model
model = tmp_to_model(out_guess);

end
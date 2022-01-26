function [model,invariants,varargout] = qtiplus_fit(data,btensors,varargin)
% function [model,invariants,varargout] = qtiplus_fit(data, btensors)
%
% This function computes the non-diffusion weighted signal S0, the mean
% diffusivity tensor D, and the covariance tensor C as implied by the QTI
% framework. During the fit, relevant positivity constraints are enforced
% on the estimated D and C tensors. The function returns the QTI model
% parameters and the invariants computed on them.
%
% For details, see 
%                  Herberthson et al (2021) NeuroImage.
%                  Westin et al (2016) NeuroImage.
%
% 
%
% INPUT:
%       required
%           - data:       a 4D matrix as [nx,ny,nz,ndiffusion_volumes]
%
%           - btensors:   a 2D matrix [#diffusion_volumes,6] where each row is a 1x6 vector
%                         identifying one of the b-tensors used in the experiment. The following
%                         convention for reducing a 3x3 Btensor to a 1x6 vector is used:
%                                       Bij = [Bxx Byy Bzz sqrt(2)*Bxy sqrt(2)*Bxz sqrt(2)*Byz]
%                         The tensors are to be input in SI units (s/m^2)
%
%       optional
%           - mask:       a mask for avoiding computations on voxels
%                         outside the brain. If not provided, one is computed.
%
%           - pipeline:   an integer specifying which steps to use (default: 0 for SDPdc)
%                         options:
%                                   - 0: SDPdc
%                                   - 1: SDPdc & NLLSdc
%                                   - 2: SDPdc & NLLSdc & m-check & SDPdcm
%                                   - 3: SDPdc & NLLSdc & SDPdcm
%                                   - 4: SPDdc & m-check & SDPdcm
%                                   - 5: SDPdc & SDPdcm
%
%           - nvox:       integer indicating how many voxel to process at
%                         once (default: 50). a warning: from experience,
%                         processing a nvox >= 100 slows the fit rather
%                         than speeding it up.
%
%           - parallel:   0 or 1 depending on whether to perform
%                         computations in parallel (default: 1)
%
%
%           - ind:        vector with 1s on measurements to be kept during
%                         the analysis, 0s on measurements to be excluded
%
%           - solver:     string containing the name of the wanted SDP
%                         solver (default : 'mosek')
%
% OUTPUT:
%           - model:      the final model parameters, in a [nx,ny,nz,28] matrix where
%                         the 28 parameters are ordered as follows:
%                                                    - model(x,y,z,1)   = S0
%                                                    - model(x,y,z,2:7) = D
%                                                    - model(x,y,z,8:28)= C
%                         The D and C tensor our output in SI units (m^2/s and m^4/s^2 respectively)
%
%           - invariants: a structure containing the invariants derived
%                         from the model parameters. The nomeclature
%                         follows the one reported in Westin et al (2016) NeuroImage
%
%           - varargout:  if a pipeline with multiple steps is selected,
%                         the intermediate model parameters are contained in this
%                         variable, ordered according to the chosen pipeline.
%
%
% Author: Deneb Boito 
% email:  deneb.boito@liu.se


% retrieve optional inputs
p = inputParser;
p.addParameter('mask', 0)
p.addParameter('pipeline', 0)
p.addParameter('nvox', 50)
p.addParameter('parallel', 1)
p.addParameter('ind', 0)
p.addParameter('solver', 'mosek')
p.parse(varargin{:})
mask        = p.Results.mask;
pipeline    = p.Results.pipeline;
nvox        = p.Results.nvox;
parallel    = p.Results.parallel;
ind         = p.Results.ind;
solver      = p.Results.solver;

% set solver for SDP problem
% it will save this preference until the current Matlab session runs
clear global
if ~strcmpi(solver, 'mosek') && ~strcmpi(solver, 'sdpt3')
    error('The entered solver name is wrong. Options:''mosek'', ''sdpt3''')
end
if strcmpi(solver, 'mosek')
    try
        cvx_solver mosek
        fprintf('Using Mosek as SDP solver \n')
        %this is apparently necessary for parfor
        cvxsolver = 'mosek';
        
    catch
        cvx_solver SDPT3
        fprintf('Mosek is not available, using SDPT3 instead \n')
        cvxsolver = 'sdpt3';
    end
else
    cvx_solver SDPT3
    fprintf('Using SDPT3 as SDP solver \n')
    cvxsolver = 'sdpt3';
end



% check if a mask has been provided
if mask == 0
    mask = simple_mask(data);
end

% if parallel computations desired, start a parpool if one is not active
% also check to see if parallel toolbox is installed
if parallel
    if license('test', 'Distrib_Computing_Toolbox')
        pp = gcp('nocreate');
        if isempty(pp)
            parpool
        end
    else
        fprintf('Parallel computing toolbox not available. Proceeding with one worker \n')
    end
end

% if ind vector has not been provided, use all measurements
if ind == 0
    ind = ones(size(btensors,1),1) > 0;
else
    ind = ind > 0;
end

% set mcheckflag to default to 1, each pipeline not using the mcheck can
% set it to 0
mcheckflag = 1;

% convert data & btensor to double
data = double(data);
btensors = double(btensors);

% decide pipeline & do fit
switch pipeline
    
    case 0 % SDPdc
        fprintf('Selected step: SDPdc \n')
        fprintf('Fitting...\n')
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        
    case 1 % SDPdc & NLLSdc
        fprintf('Selected steps: SDPdc & NLLSdc \n')
        fprintf('Fitting...\n')
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        varargout{1} = model;
        model = qtip_pipe_NLLSdc(model,data,btensors,mask,ind,parallel);
        
    case 2 % SDPdc & NLLSdc & m-check & SDPdcm
        fprintf('Selected steps: SDPdc & NLLSdc & m-check & SDPdcm \n')
        fprintf('Fitting...\n')
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        varargout{1} = model;
        model = qtip_pipe_NLLSdc(model,data,btensors,mask,ind,parallel);
        varargout{2} = model;
        model = qtip_pipe_SDPdcm(model,data,btensors,mask,nvox,mcheckflag,ind,parallel,cvxsolver);
        
    case 3 % SDPdc & NLLSdc & SDPdcm
        fprintf('Selected steps: SDPdc & NLLSdc & SDPdcm \n')
        fprintf('Fitting...\n')
        mcheckflag = 0;
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        varargout{1} = model;
        model = qtip_pipe_NLLSdc(model,data,btensors,mask,ind,parallel);
        varargout{2} = model;
        model = qtip_pipe_SDPdcm(model,data,btensors,mask,nvox,mcheckflag,ind,parallel,cvxsolver);
        
    case 4 % SDPdc & m-check & SDPdcm
        fprintf('Selected steps: SDPdc & SDPdcm (with m-check) \n')
        fprintf('Fitting...\n')
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        varargout{1} = model;
        model = qtip_pipe_SDPdcm(model,data,btensors,mask,nvox,mcheckflag,ind,parallel,cvxsolver);
        
    case 5 % SDPdc & SDPdcm
        fprintf('Selected steps: SDPdc & SDPdcm \n')
        fprintf('Fitting...\n')
        mcheckflag = 0;
        model = qtip_pipe_SDPdc(data,btensors,mask,nvox,ind,parallel,cvxsolver);
        varargout{1} = model;
        model = qtip_pipe_SDPdcm(model,data,btensors,mask,nvox,mcheckflag,ind,parallel,cvxsolver);
        
end

% compute invariants
invariants = compute_invariants(model);
fprintf('...Done!\n')
end
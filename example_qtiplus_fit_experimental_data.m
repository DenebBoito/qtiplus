%% Example to reproduce partial results from Herberthosn et al. 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script guides through the replication of part of the results       %
% presented in Herberthson et al. 2021. To run this script it is necessary% 
% to download the data published in Szczepankiewic et al. 2019.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Download, load, and prepare the data
% Download, unzip, and add to the Matlab path the data from  
% https://github.com/filip-szczepankiewicz/Szczepankiewicz_DIB_2019/tree/master/DATA/brain/MD-dMRI
% The data we are going to use here are the "BRAIN_FWF_MERGED_mc.nii.gz" 
% (or "BRAIN_FWF_MERGED_mc.nii" depending on the unzipping procedure)
% which is contained in the MD-dMRI.zip file
% The b-tensors are saved in the "BRAIN_FWF_MERGED_mc_xps.mat"

% load the data and the experimental parameters (xps)
% data = niftiread('BRAIN_FWF_MERGED_mc.nii.gz');
data = niftiread('BRAIN_FWF_MERGED_mc.nii');
load('BRAIN_FWF_MERGED_mc_xps.mat')

% extract the b-tensors from the xps structure
btens = xps.bt;

% first we order the data according to b-value and tensor-encoding shape
% for each b-shell, the data will be in the order LTE, PTE, STE.
load('indx_order')
data = data(:,:,:,indx_order);
btens = btens(indx_order,:);

% load the indexes used to create the various datasets in Herberthson et.
% al 2021

% p217
load('indx_p217')

% p81
load('indx_p81')

% p56
load('indx_p56')

% p39
load('indx_p39')


%% Compute the model parameters for all the datasets
% here we compute some of the fit to then reproduce part of the results
% in Herberthson et al. 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these computation will take some time!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a mask
mask = simple_mask(data);

% for the sake of this example, we will only fit one slice of the data.
% we pick slice 13.
mask(:,:,1:12) = 0;
mask(:,:,14:end) = 0;


% fit the various protocols with SDP(dcm)
% use the 'ind' keyword to specify which volumes to use for the fit
% use the 'pipeline' keyword to specify which steps of the QTI+ framemork
% to perform
% use the 'mask' keyword to pass the mask to the function                                                                                                                                              
[model_p217, invariants_p217] = qtiplus_fit(data, btens,...
                                                     'ind', indx_p217,...
                                                     'mask', mask, ...
                                                     'pipeline', 2);
                                                
[model_p81, invariants_p81] = qtiplus_fit(data, btens,...
                                                   'ind', indx_p81,...
                                                   'mask', mask, ...
                                                   'pipeline', 2);
                                                 
[model_p56, invariants_p56, model_SDPdc_p56, model_NLLSdc_p56] = ...
                          qtiplus_fit(data, btens, 'ind', indx_p56,...
                                                   'mask', mask, ...
                                                   'pipeline', 2);
                                                
[model_p39, invariants_p39] = qtiplus_fit(data, btens, ...
                                                   'ind', indx_p39,...
                                                   'mask', mask, ...
                                                   'pipeline', 2);                                                
%% recreate some of the figures

% figure 6, last three rows
invariants_p56_SDPdc = compute_invariants(model_SDPdc_p56);
invariants_p56_NLLSdc = compute_invariants(model_NLLSdc_p56);
invs_fig6{1} = invariants_p56_SDPdc;
invs_fig6{2} = invariants_p56_NLLSdc;
invs_fig6{3} = invariants_p56;
plot_qti_invariants(invs_fig6);

% figure 7b, first 4 rows
invs_fig7b{1} = invariants_p217;
invs_fig7b{2} = invariants_p81;
invs_fig7b{3} = invariants_p56;
invs_fig7b{4} = invariants_p39;
plot_qti_invariants(invs_fig7b);

% figure 8a, second row
plot_qti_invariant(invs_fig7b, "FA", [0 1])

% figure 9a, second row
plot_qti_invariant(invs_fig7b, "uFA", [0 1])

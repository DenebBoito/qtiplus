% this script shows how to use the qtiplus software

% the example dataset is based on signals rising from a non-central
% Wishart distribution whose parameters are set as explained in the
% QTI+ paper by Herberthson et al. Here we used the same parameters as in
% the paper, i.e. one distribution originates from an isotropic mean
% diffusion tensor, while the second distribution is based on an
% anisotropic mean diffusion tensor. Rician distributed signals are then
% obtained by adding Gaussian noise to the real and imaginary parts of the
% noiseless data. 
%
% The dataset has matrix size [64 x 64 x 6 x 217], where the last dimension
% indicates the number of diffusion volumes. Each of the 6 slices contains,
% respectively:
% 1 - Noiseless signals for the isotropic mean diffusion tensor
% 2 - Signal from 1) with SNR = 30
% 3 - Signal from 1) with SNR = 15
% 4 - Noiseless signals for the anisotropic mean diffusion tensor
% 5 - Signal from 4) with SNR = 30
% 6 - Signal from 4) with SNR = 15

clear 

% load the example dataset and the b-tensors
% the b-tensors in this example are part of the protocol from
% Szczepankiewicz et al. Data in Brief 2019
% 
load('bten');
data = niftiread('example_dataset.nii');

% create a mask to limit the computation to portion of the dataset
% by selecting a [(2*range)+1 x (2*range)+1 ] window
range = 5;
mask = zeros(64,64,6); 
mask(end/2-range:end/2+range, end/2-range:end/2+range,:) = 1;

% the first 13 volumes are "acquired" without diffusion weighting, so
% we exclude them from the measurements to be used during the fitting
% ind == 1 -> keep
% ind == 0 -> discard
ind = ones(217,1);
ind(1:13) = 0;

% perform the fit
% we pass the indices of the measurements we want to keep during the fit
% and the mask as key-value pairs
[model,invariants] = qtiplus_fit(data, bten, ...
                                 'ind', ind, ...
                                 'mask', mask);


%% look at the uFA maps
% the maps for the first and fourth slice can be considered as ground truth
% for what the QTI model can return for the signals gerenated from a non
% central wishart distribution with the adopted parameters.
% here we look at a central window of the micro-FA for the two cases,
% isotropic and anisotropic mean diffusion tensor
uFA = invariants.uFA(end/2-range:end/2+range,end/2-range:end/2+range,:);
gt_iso = uFA(1,1,1);
gt_ani = uFA(1,1,4);
x_values = 0:0.01:1;


figure()
subplot(2,4,1)
imagesc(uFA(:,:,1), [0 1]), colormap gray, title('GT'),...
    ylabel('Isotropic $\widehat{\mathbf{D}}_{ij}$', 'Interpreter','latex')
ax = gca;
ax.YTick = [];
ax.YTickLabel = [];
ax.XTick = [];
ax.XTickLabel = [];
ax.FontSize = 20;


subplot(2,4,2)
imagesc(uFA(:,:,2), [0 1]), colormap gray, title('SNR = 30',...
    'FontSize', 20), axis off

subplot(2,4,3)
imagesc(uFA(:,:,3), [0 1]), colormap gray, title('SNR = 15',...
    'FontSize', 20), axis off


subplot(2,4,4)
tmp = uFA(:,:,2);
pd1 = fitdist(tmp(:),'Kernel');
y1 = pdf(pd1,x_values);
tmp = uFA(:,:,3);
pd2 = fitdist(tmp(:),'Kernel');
y2 = pdf(pd2,x_values);
plot(x_values, y1, 'Color', [1 0.6 0], 'LineWidth', 2),
hold on,
plot(x_values, y2, 'Color', [0 0.6 0.7], 'LineWidth', 2)
ax = gca;
ax.YTick = [];
ax.YTickLabel = [];
ax.FontSize = 20;
yl = max(ax.YLim);
line([gt_iso gt_iso],[0 yl], 'Color', [0.3 0.3 0.3], 'LineWidth', 2,...
    'Linestyle', '--')
title('\muFA value distribution')


subplot(2,4,5)
imagesc(uFA(:,:,4), [0 1]), colormap gray,...
    ylabel('Anisotropic $\widehat{\mathbf{D}}_{ij}$', 'Interpreter','latex'), 
ax = gca;
ax.YTick = [];
ax.YTickLabel = [];
ax.XTick = [];
ax.XTickLabel = [];
ax.FontSize = 20;
colorbar('southoutside')


subplot(2,4,6)
imagesc(uFA(:,:,5), [0 1]), colormap gray, axis off

subplot(2,4,7)
imagesc(uFA(:,:,6), [0 1]), colormap gray, axis off

subplot(2,4,8)
tmp = uFA(:,:,5);
pd1 = fitdist(tmp(:),'Kernel');
y1 = pdf(pd1,x_values);
tmp = uFA(:,:,6);
pd2 = fitdist(tmp(:),'Kernel');
y2 = pdf(pd2,x_values);
plot(x_values, y1, 'Color', [1 0.6 0], 'LineWidth', 2),
hold on,
plot(x_values, y2, 'Color', [0 0.6 0.7], 'LineWidth', 2)
ax = gca;
ax.YTick = [];
ax.YTickLabel = [];
ax.FontSize = 20;
yl = max(ax.YLim);
line([gt_ani gt_ani],[0 yl], 'Color', [0.3 0.3 0.3], 'LineWidth', 2,...
    'Linestyle', '--')
legend('SNR = 30', 'SNR = 15', 'GT', 'Location', 'northwest', 'box', 'off')
 
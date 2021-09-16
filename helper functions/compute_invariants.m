function invariants = compute_invariants(model)
% function compute_invariants()
%
% returns the invariants derived from the estimated D and C tensors

% reshape model for convenience
siz      = size(model);
mreshape = reshape(model, prod(siz(1:3)), 28);

% extract D and C & change units
D = mreshape(:,2:7) * 1e9;
C = mreshape(:,8:28)* 1e18;

% define matrices for contracting tensors
[Eiso,Ebulk,Eshear] = get_contracting_matrices();
d_1x21              = convert_1x6_to_1x21(D);
m_1x21              = C + d_1x21; % second moment tensor

%% S0
invariants.s0 = squeeze(model(:,:,:,1)); 

%% FA, MD, RD, AD, ax_dir
[invariants.FA, invariants.AD, invariants.RD, invariants.ax_dir, invariants.FA_col] = compute_lamda_params(D,siz);
invariants.MD = 1/3 * sum(model(:,:,:,2:4),4) * 1e9; % change units to um2/ms


%% Variances and Normalized Variances
% Nomenclature as in Westin et al. Neuroimage 2016

invariants.V_MD    = reshape(C * Ebulk', siz(1), siz(2), siz(3));
invariants.V_iso   = reshape(C * Eiso', siz(1), siz(2), siz(3));
invariants.V_shear = invariants.V_iso - invariants.V_MD;

invariants.C_MD = invariants.V_MD ./ reshape(max(m_1x21 * Ebulk', eps), siz(1), siz(2), siz(3));
invariants.C_mu = (3/2) * reshape(m_1x21 * Eshear' ./ max(m_1x21 * Eiso', eps), siz(1), siz(2), siz(3));
invariants.C_M  = (3/2) * reshape(d_1x21 * Eshear' ./ max(d_1x21 * Eiso', eps), siz(1), siz(2), siz(3));
invariants.C_c  = invariants.C_M ./ max(invariants.C_mu, eps);

invariants.uFA  = real(sqrt(invariants.C_mu)); % this shouldn't happen, but better be prepared

invariants.OP2 = reshape(d_1x21 * Eshear' ./ max(m_1x21 * Eiso', eps), siz(1), siz(2), siz(3)); % order parameter
invariants.OP  = sqrt(invariants.OP2);

%% Kurtosis
invariants.K_bulk  = 3 * invariants.V_MD ./ reshape(max(d_1x21 * Ebulk', eps), siz(1), siz(2), siz(3));
invariants.K_shear = (6/5) *  reshape(C * Eshear' ./ max(d_1x21 * Ebulk', eps), siz(1), siz(2), siz(3));
invariants.K_mu    = (6/5) *  reshape(m_1x21 * Eshear' ./ max(d_1x21 * Ebulk', eps), siz(1), siz(2), siz(3));
invariants.MK      = invariants.K_bulk + invariants.K_shear;
invariants.MKt     = invariants.K_bulk + invariants.K_mu;                   % 

end

function [Eiso,Ebulk,Eshear,Etsym] = get_contracting_matrices()
% returns matrices for computing invariants

Eiso   = 1/3 * eye(6);

Ebulk  = 1/9 * blkdiag(ones(3,3), zeros(3,3));

Eshear = Eiso - Ebulk;

Etsym  = Ebulk + 2/5 * Eshear;

% convert to 1x21 for calculations
Eiso = convert_6x6_to_1x21(Eiso);
Ebulk = convert_6x6_to_1x21(Ebulk);
Eshear = convert_6x6_to_1x21(Eshear);
Etsym = convert_6x6_to_1x21(Etsym);
end
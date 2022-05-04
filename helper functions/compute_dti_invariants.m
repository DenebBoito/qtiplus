function invariants = compute_dti_invariants(model)
siz = size(model);
mreshape = reshape(model, prod(siz(1:3)), 7);
D = mreshape(:,2:7) * 1e9;
invariants.s0 = squeeze(model(:,:,:,1));
[invariants.FA, invariants.AD, invariants.RD, invariants.ax_dir, invariants.FA_col] = compute_lamda_params(D,siz);
invariants.MD = clip_invariant_min_max(1/3 * sum(model(:,:,:,2:4),4) * 1e9, [0 4]); % change units to um2/ms
end
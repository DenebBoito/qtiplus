function out = cholesky_factorization(in)
% function cholesky_factorization()
try
    out = chol(in,'upper');
catch
    [L,Di] = mchol(in);
    try
        out = chol(L*Di*L','upper');
    catch
        out = (L*sqrtm(Di))';
    end
end
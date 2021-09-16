function [L,D]= mchol(G)
% function mchol(G)
%  
% returns the modified cholesky representation of a matrix G
% 
% the implementation follows the steps of the algorithm at page 111 of
%    "Practical Optimization" 1982
% by Philip E. Gill, Walter Murray, and Margareth H. Wright 
% ...

% MC1 [Compute the bound on the elements of the factors.]
n    = max(size(G));
v    = max([1, sqrt(n^2-1)]);
indx = eye(n) > 0;
gamma   = max(G(indx)); % maximum element of the diagonal of G
Epsilon = max(G((1-indx)>0)); % maximum off diagonal element of G
beta2   = max([gamma, Epsilon/v, eps]);

% MC2 [Initialize.]
C = diag(G(indx)); % set c_{ii} to g_{ii}

L = zeros(n);
D = zeros(n);
theta = zeros(1,n);

for j = 1:n % column index
    
    % set indexes to be used for calculations
    s = 1 : j-1;
    i = j+1 : n;

    % MC4 [Compute the j-th row of L. Find the maximum modulus of lijdj.]
    % set l_{js} = c_{js} / d_{s}
    if (j > 1)
        L(j,s) = C(j,s) ./ diag( D(s,s) )';
    end
    
    % c_{ij} = g_{ij} - sum_{s=1}^{j-1} l_{js} c_{is}
    if (j >= 2)
        if (j < n )
            C(i,j) = G(i,j) - ( L(j,s) * C(i,s)' )';
        end
    else
        C(i,j)=G(i,j);
    end
    
    % set theta_{j} = max_{i} |c_{ij}|
    if (j == n)
        theta(j) = 0;
    else
        theta(j) = max( abs( C(i,j) ) );
    end

    % MC5 [Compute the j-th diagonal element of D.]
    D(j,j) = max( [ eps, abs(C(j,j)), theta(j)^2/beta2 ] );
    
    % MC6 [Update the prospective diagonal elements and the column index.]
    % skip diagonal modification of E as we do not need it
    C(i,i) = C(i,i) - C(i,j).^2 ./ D(j,j);
    
end

% set diagonal element of L to 1
L(indx) = 1;

function out = qtipm_mcheck(model,D0)
% function qtipm_mcheck()
%
% checks if condition (m-) is satisfied using SDP
%
% assumes model is [28,nvox]

nvox = size(model,2);

% get d and c and change units
d = arrayfun(@(i)convert_1x6_to_3x3(model(2:7,i).*1e9), 1:nvox, 'UniformOutput', false);
c = arrayfun(@(i)convert_1x21_to_3x3x3x3(model(8:28,i).*1e18), 1:nvox,'UniformOutput', false);

% get D0^2 * dijdkl
dij = eye(3);
dij_1x6 = tm_3x3_to_1x6(dij);
DijDkl = D0^2 * tm_1x21_to_3x3x3x3(tm_1x6_to_1x21(dij_1x6));

cvx_begin sdp quiet

variable l(9,nvox) % independent lambdas

minimize 0    % feasibility problem: no objective function

subject to
for i = 1:nvox
    [[-c{i}(1, 1, 1, 1) - d{i}(1, 1)^2 + DijDkl(1,1,1,1), -c{i}(1, 1, 1, 2) - d{i}(1, 1) * d{i}(1, 2) + DijDkl(1,1,1,2), -c{i}(1, 1, 1, 3) - d{i}(1, 1) * d{i}(1, 3) + DijDkl(1,1,1,3),  -c{i}(1, 1, 1, 2) - d{i}(1, 1) * d{i}(1, 2) + DijDkl(1,1,1,2), l(1,i), l(2,i), -c{i}(1, 1, 1, 3) - d{i}(1, 1) * d{i}(1, 3) + DijDkl(1,1,1,3), l(3,i), l(4,i)];
    [-c{i}(1, 1, 1, 2) - d{i}(1, 1) * d{i}(1, 2) + DijDkl(1,1,1,2), -c{i}(1, 1, 2, 2) - d{i}(1, 1) * d{i}(2, 2) + DijDkl(1,1,2,2), -c{i}(1, 1, 2, 3) - d{i}(1, 1) * d{i}(2, 3) + DijDkl(1,1,2,3), -l(1,i) - 2 * c{i}(1, 2, 1, 2) - 2 * d{i}(1, 2)^2 + DijDkl(1,2,1,2), -c{i}(1, 2, 2, 2) - d{i}(1, 2) * d{i}(2, 2) + DijDkl(1,2,2,2), l(5,i), -l(3,i) - 2 * c{i}(1, 2, 1, 3) - 2 * d{i}(1, 2) * d{i}(1, 3) + DijDkl(1,2,1,3), -c{i}(1, 3, 2, 2) - d{i}(1, 3) * d{i}(2, 2) + DijDkl(1,3,2,2), l(6,i)];
    [-c{i}(1, 1, 1, 3) - d{i}(1, 1) * d{i}(1, 3) + DijDkl(1,1,1,3), -c{i}(1, 1, 2, 3) - d{i}(1, 1) * d{i}(2, 3) + DijDkl(1,1,2,3), -c{i}(1, 1, 3, 3) - d{i}(1, 1) * d{i}(3, 3) + DijDkl(1,1,3,3), -l(2,i) - 2 * c{i}(1, 2, 1, 3) - 2 * d{i}(1, 2) * d{i}(1, 3) + DijDkl(1,2,1,3), -l(5,i) - 2 * c{i}(1, 2, 2, 3) - 2 * d{i}(1, 2) * d{i}(2, 3) + DijDkl(1,2,2,3), -c{i}(1, 2, 3, 3) - d{i}(1, 2) * d{i}(3, 3) + DijDkl(1,2,3,3), -l(4,i) - 2 * c{i}(1, 3, 1, 3) - 2 * d{i}(1, 3)^2 + DijDkl(1,3,1,3), -l(6,i) - 2 * c{i}(1, 3, 2, 3) - 2 * d{i}(1, 3) * d{i}(2, 3) + DijDkl(1,3,2,3), -c{i}(1, 3, 3, 3) - d{i}(1, 3) * d{i}(3, 3) + DijDkl(1,3,3,3)];
    [-c{i}(1, 1, 1, 2) - d{i}(1, 1) * d{i}(1, 2) + DijDkl(1,1,1,2), -l(1,i) - 2 * c{i}(1, 2, 1, 2) - 2 * d{i}(1, 2)^2 + DijDkl(1,2,1,2), -l(2,i) - 2 * c{i}( 1, 2, 1, 3) - 2 * d{i}(1, 2) * d{i}(1, 3) + DijDkl(1,2,1,3), -c{i}(1, 1, 2, 2) - d{i}(1, 1) * d{i}(2, 2) + DijDkl(1,1,2,2), -c{i}(1, 2, 2, 2) - d{i}(1, 2) * d{i}(2, 2) + DijDkl(1,2,2,2), -c{i}(1, 3, 2, 2) - d{i}(1, 3) * d{i}(2, 2) + DijDkl(1,3,2,2), -c{i}(1, 1, 2, 3) - d{i}(1, 1) * d{i}(2, 3) + DijDkl(1,1,2,3), +l(7,i), +l(8,i)];
    [+l(1,i),-c{i}(1, 2, 2, 2) - d{i}(1, 2) * d{i}(2, 2) + DijDkl(1,2,2,2), -l(5,i) - 2 * c{i}(1, 2, 2, 3) - 2 * d{i}(1, 2) * d{i}(2, 3) + DijDkl(1,2,2,3), -c{i}(1, 2, 2, 2) - d{i}(1, 2) * d{i}(2, 2) + DijDkl(1,2,2,2), -c{i}(2, 2, 2, 2) - d{i}(2, 2)^2 + DijDkl(2,2,2,2), -c{i}(2, 2, 2, 3) - d{i}(2, 2) * d{i}(2, 3) + DijDkl(2,2,2,3), -l(7,i) - 2 * c{i}(1, 2, 2, 3) - 2 * d{i}(1, 2) * d{i}(2, 3) + DijDkl(1,2,2,3), -c{i}(2, 2, 2, 3) - d{i}(2, 2) * d{i}(2, 3) + DijDkl(2,2,2,3), +l(9,i)];
    [+l(2,i), +l(5,i),-c{i}(1, 2, 3, 3) - d{i}(1, 2) * d{i}(3, 3) + DijDkl(1,2,3,3), -c{i}(1, 3, 2, 2) - d{i}(1, 3) * d{i}(2, 2) + DijDkl(1,3,2,2),  -c{i}(2, 2, 2, 3) - d{i}(2, 2) * d{i}(2, 3) + DijDkl(2,2,2,3), -c{i}(2, 2, 3, 3) - d{i}(2, 2) * d{i}(3, 3) + DijDkl(2,2,3,3), -l(8,i) - 2 * c{i}(1, 3, 2, 3) - 2 * d{i}(1, 3) * d{i}( 2, 3) + DijDkl(1,3,2,3), -l(9,i) - 2 * c{i}(2, 3, 2, 3) - 2 * d{i}(2, 3)^2 + DijDkl(2,3,2,3), -c{i}(2, 3, 3, 3) - d{i}(2, 3) * d{i}(3, 3) + DijDkl(2,3,3,3)];
    [-c{i}(1, 1, 1, 3) - d{i}(1, 1) * d{i}(1, 3) + DijDkl(1,1,1,3), -l(3,i) - 2 * c{i}(1, 2, 1, 3) - 2 * d{i}(1, 2) * d{i}(1, 3) + DijDkl(1,2,1,3), -l(4,i) - 2 * c{i}(1, 3, 1, 3) - 2 * d{i}(1, 3)^2 + DijDkl(1,3,1,3), -c{i}(1, 1, 2, 3) - d{i}(1, 1) * d{i}(2, 3) + DijDkl(1,1,2,3), -l(7,i) - 2 * c{i}(1, 2, 2, 3) - 2 * d{i}(1, 2) * d{i}(2, 3) + DijDkl(1,2,2,3), -l(8,i) - 2 * c{i}(1, 3, 2, 3) - 2 * d{i}(1, 3) * d{i}(2, 3) + DijDkl(1,3,2,3), -c{i}(1, 1, 3, 3) - d{i}(1, 1) * d{i}(3, 3) + DijDkl(1,1,3,3), -c{i}(1, 2, 3, 3) - d{i}(1, 2) * d{i}(3, 3) + DijDkl(1,2,3,3), -c{i}(1, 3, 3, 3) - d{i}(1, 3) * d{i}(3, 3) + DijDkl(1,3,3,3)];
    [+l(3,i), -c{i}(1, 3, 2, 2) - d{i}(1, 3) * d{i}(2, 2) + DijDkl(1,3,2,2), -l(6,i) - 2 * c{i}(1, 3, 2, 3) - 2 * d{i}(1, 3) * d{i}(2, 3) + DijDkl(1,3,2,3), +l(7,i),  -c{i}(2, 2, 2, 3) - d{i}(2, 2) * d{i}(2, 3) + DijDkl(2,2,2,3), -l(9,i) - 2 * c{i}(2, 3, 2, 3) - 2 * d{i}(2, 3)^2 + DijDkl(2,3,2,3), -c{i}(1, 2, 3, 3) - d{i}(1, 2) * d{i}(3, 3) + DijDkl(1,2,3,3), -c{i}(2, 2, 3, 3) - d{i}(2, 2) * d{i}(3, 3) + DijDkl(2,2,3,3), -c{i}(2, 3, 3, 3) - d{i}(2, 3) * d{i}(3, 3) + DijDkl(2,3,3,3)];
    [+l(4,i), +l(6,i), -c{i}(1, 3, 3, 3) - d{i}(1, 3) * d{i}(3, 3) + DijDkl(1,3,3,3), +l(8,i), +l(9,i), -c{i}(2, 3, 3, 3) - d{i}(2, 3) * d{i}(3, 3) + DijDkl(2,3,3,3), -c{i}(1, 3, 3, 3) - d{i}(1, 3) * d{i}(3, 3) + DijDkl(1,3,3,3), -c{i}(2, 3, 3, 3) - d{i}(2, 3) * d{i}(3, 3) + DijDkl(2,3,3,3), -c{i}(3, 3, 3, 3) - d{i}(3, 3)^2 + DijDkl(3,3,3,3)]] >=0 ;
end
cvx_end
 
 % feasibility problem in cvx: if solver finds a feasible solution, cvx_optval will be 0
% otherwise cvx_optval will be +Inf
if cvx_optval == 0
    out = 0;
else
    out = 1;
end 
end
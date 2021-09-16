function out = convert_3x3_to_1x6(in)
% function convert_3x3_to_1x6
%
% converts a tensor from a 3x3 to a 1x6 representation
% picks up the upper triangular elements to be used in the Cholesky
% formulation

xx = in(1,1);
yy = in(2,2);
zz = in(3,3);
xy = in(1,2);
xz = in(1,3);
yz = in(2,3);

f = sqrt(2);
out = [xx yy zz f * [xy xz yz]];

end
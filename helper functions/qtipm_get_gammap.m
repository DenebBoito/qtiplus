function gamma_p = qtipm_get_gammap(c,l,D0)
%
% returns the 9x9 form of the matrix Gamma plus from the C 3x3x3x3 tensor
% C is expected to have already the correct units
%

gij = eye(3);
gij_1x6 = convert_3x3_to_1x6(gij);
scale = 0.25 * (D0^2);
d02gijgkl = scale * convert_1x21_to_3x3x3x3(convert_1x6_to_1x21(gij_1x6));


% this is d02gijgkl + Cijkl 
 gamma_p = [
             [+c(1, 1, 1, 1) + d02gijgkl(1,1,1,1), +c(1, 1, 1, 2) + d02gijgkl(1,1,1,2), +c(1, 1, 1, 3) + d02gijgkl(1,1,1,3),...
             +c(1, 1, 1, 2) + d02gijgkl(1,1,1,2), l(1), l(2), +c(1, 1, 1, 3) + d02gijgkl(1,1,1,3), l(3), l(4)];
     
 
             [+c(1, 1, 1, 2) + d02gijgkl(1,1,1,2), +c(1, 1, 2, 2) + d02gijgkl(1,1,2,2), +c(1, 1, 2, 3) + d02gijgkl(1,1,2,3),...
             -l(1) + 2 * c(1, 2, 1, 2) + d02gijgkl(1,2,1,2), +c(1, 2, 2, 2) + d02gijgkl(1,2,2,2), l(5), ...
             -l(3) + 2 * c(1, 2, 1, 3) + d02gijgkl(1,2,1,3), +c(1, 3, 2, 2) + d02gijgkl(1,3,2,2), l(6)];
     
             [+c(1, 1, 1, 3) + d02gijgkl(1,1,1,3), +c(1, 1, 2, 3) + d02gijgkl(1,1,2,3), +c(1, 1, 3, 3) + d02gijgkl(1,1,3,3),...
             -l(2) + 2 * c(1, 2, 1, 3) + d02gijgkl(1,2,1,3), -l(5) + 2 * c(1, 2, 2, 3) + d02gijgkl(1,2,2,3),...
             +c(1, 2, 3, 3) + d02gijgkl(1,2,3,3), -l(4) + 2 * c(1, 3, 1, 3) + d02gijgkl(1,3,1,3),...
             -l(6) + 2 * c(1, 3, 2, 3) + d02gijgkl(1,3,2,3), +c(1, 3, 3, 3) + d02gijgkl(1,3,3,3)];
     
            [+c(1, 1, 1, 2) + d02gijgkl(1,1,1,2), -l(1) + 2 * c(1, 2, 1, 2) + d02gijgkl(1,2,1,2),...
            -l(2) + 2 * c( 1, 2, 1, 3) + d02gijgkl(1,2,1,3), +c(1, 1, 2, 2) + d02gijgkl(1,1,2,2), ...
            +c(1, 2, 2, 2) + d02gijgkl(1,2,2,2), +c(1, 3, 2, 2) + d02gijgkl(1,3,2,2), +c(1, 1, 2, 3) + d02gijgkl(1,1,2,3), +l(7), +l(8)];
     
            [+l(1),+c(1, 2, 2, 2) + d02gijgkl(1,2,2,2), -l(5) + 2 * c(1, 2, 2, 3) + d02gijgkl(1,2,2,3), +c(1, 2, 2, 2) + d02gijgkl(1,2,2,2),...
            +c(2, 2, 2, 2) + d02gijgkl(2,2,2,2), +c(2, 2, 2, 3) + d02gijgkl(2,2,2,3), ...
            -l(7) + 2 * c(1, 2, 2, 3) + d02gijgkl(1,2,2,3), +c(2, 2, 2, 3) + d02gijgkl(2,2,2,3), +l(9)];
     
            [+l(2), +l(5),+c(1, 2, 3, 3) + d02gijgkl(1,2,3,3), +c(1, 3, 2, 2) + d02gijgkl(1,3,2,2),...
            +c(2, 2, 2, 3) + d02gijgkl(2,2,2,3), +c(2, 2, 3, 3) + d02gijgkl(2,2,3,3),...
            -l(8) + 2 * c(1, 3, 2, 3) + d02gijgkl(1,3,2,3), -l(9) + 2 * c(2, 3, 2, 3) + d02gijgkl(2,3,2,3), +c(2, 3, 3, 3) + d02gijgkl(2,3,3,3)];
     
            [+c(1, 1, 1, 3) + d02gijgkl(1,1,1,3), -l(3) + 2 * c(1, 2, 1, 3) + d02gijgkl(1,2,1,3),...
            -l(4) + 2 * c(1, 3, 1, 3) + d02gijgkl(1,3,1,3), +c(1, 1, 2, 3) + d02gijgkl(1,1,2,3), ...
            -l(7) + 2 * c(1, 2, 2, 3) + d02gijgkl(1,2,2,3), -l(8) + 2 * c(1, 3, 2, 3) + d02gijgkl(1,3,2,3),...
            +c(1, 1, 3, 3) + d02gijgkl(1,1,3,3), +c(1, 2, 3, 3) + d02gijgkl(1,2,3,3), +c(1, 3, 3, 3) + d02gijgkl(1,3,3,3)];
     
            [+l(3), +c(1, 3, 2, 2) + d02gijgkl(1,3,2,2), -l(6) + 2 * c(1, 3, 2, 3) + d02gijgkl(1,3,2,3), +l(7),...
            +c(2, 2, 2, 3) + d02gijgkl(2,2,2,3), -l(9) + 2 * c(2, 3, 2, 3) + d02gijgkl(2,3,2,3), ...
            +c(1, 2, 3, 3) + d02gijgkl(1,2,3,3), +c(2, 2, 3, 3) + d02gijgkl(2,2,3,3), +c(2, 3, 3, 3) + d02gijgkl(2,3,3,3)];
     
            [+l(4), +l(6), +c(1, 3, 3, 3) + d02gijgkl(1,3,3,3), +l(8), +l(9), +c(2, 3, 3, 3) + d02gijgkl(2,3,3,3),...
            +c(1, 3, 3, 3) + d02gijgkl(1,3,3,3), +c(2, 3, 3, 3) + d02gijgkl(2,3,3,3), +c(3, 3, 3, 3) + d02gijgkl(3,3,3,3)]
            
            ];
end
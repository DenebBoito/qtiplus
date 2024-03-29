function m_tilde = get_M_9x9_from_d_c_SL(d,c,l,D0)
%
% returns the 9x9 form of the matrix M from the 3x3 and 3x3x3x3 tensors
%

dij = eye(3);
dij_1x6 = convert_3x3_to_1x6(dij);
DijDkl = D0^2 * convert_1x21_to_3x3x3x3(convert_1x6_to_1x21(dij_1x6));

% this is D0^2IijIkl - Mijkl 
% switching to Magnus convetion for the lambdas
  m_tilde = [...
   
    [-c(1, 1, 1, 1) - d(1, 1)^2 + DijDkl(1,1,1,1),... 
     -c(1, 1, 1, 2) - d(1, 1) * d(1, 2) + DijDkl(1,1,1,2),... 
     -c(1, 1, 1, 3) - d(1, 1) * d(1, 3) + DijDkl(1,1,1,3),... 
     -c(1, 1, 1, 2) - d(1, 1) * d(1, 2) + DijDkl(1,1,1,2),... 
     l(1),...
     l(2),... 
     -c(1, 1, 1, 3) - d(1, 1) * d(1, 3) + DijDkl(1,1,1,3),...
     l(3),...
     l(4)];...
     
     
     [-c(1, 1, 1, 2) - d(1, 1) * d(1, 2) + DijDkl(1,1,1,2),...
     -c(1, 1, 2, 2) - d(1, 1) * d(2, 2) + DijDkl(1,1,2,2),...
     -c(1, 1, 2, 3) - d(1, 1) * d(2, 3) + DijDkl(1,1,2,3),...
     -l(1) - 2 * c(1, 2, 1, 2) - 2 * d(1, 2)^2 + 2 * DijDkl(1,2,1,2),...
     -c(1, 2, 2, 2) - d(1, 2) * d(2, 2) + DijDkl(1,2,2,2),...
     l(5),...
     -l(3) - 2 * c(1, 2, 1, 3) - 2 * d(1, 2) * d(1, 3) + 2 * DijDkl(1,2,1,3),...
     -c(1, 3, 2, 2) - d(1, 3) * d(2, 2) + DijDkl(1,3,2,2), l(6)];...
     
     
     [-c(1, 1, 1, 3) - d(1, 1) * d(1, 3) + DijDkl(1,1,1,3),...
     -c(1, 1, 2, 3) - d(1, 1) * d(2, 3) + DijDkl(1,1,2,3),...
     -c(1, 1, 3, 3) - d(1, 1) * d(3, 3) + DijDkl(1,1,3,3),...
     -l(2) - 2 * c(1, 2, 1, 3) - 2 * d(1, 2) * d(1, 3) + 2 * DijDkl(1,2,1,3),...
     -l(5) - 2 * c(1, 2, 2, 3) - 2 * d(1, 2) * d(2, 3) + 2 * DijDkl(1,2,2,3),...
     -c(1, 2, 3, 3) - d(1, 2) * d(3, 3) + DijDkl(1,2,3,3),...
     -l(4) - 2 * c(1, 3, 1, 3) - 2 * d(1, 3)^2 + 2 * DijDkl(1,3,1,3),...
     -l(6) - 2 * c(1, 3, 2, 3) - 2 * d(1, 3) * d(2, 3) + 2 * DijDkl(1,3,2,3),...
     -c(1, 3, 3, 3) - d(1, 3) * d(3, 3) + DijDkl(1,3,3,3)];...
     
     
     [-c(1, 1, 1, 2) - d(1, 1) * d(1, 2) + DijDkl(1,1,1,2),...
     -l(1) - 2 * c(1, 2, 1, 2) - 2 * d(1, 2)^2 + 2 * DijDkl(1,2,1,2),...
     -l(2) - 2 * c( 1, 2, 1, 3) - 2 * d(1, 2) * d(1, 3) + 2 * DijDkl(1,2,1,3),...
     -c(1, 1, 2, 2) - d(1, 1) * d(2, 2) + DijDkl(1,1,2,2),...
     -c(1, 2, 2, 2) - d(1, 2) * d(2, 2) + DijDkl(1,2,2,2),...
     -c(1, 3, 2, 2) - d(1, 3) * d(2, 2) + DijDkl(1,3,2,2),...
     -c(1, 1, 2, 3) - d(1, 1) * d(2, 3) + DijDkl(1,1,2,3),...
     +l(7),...
     +l(8)];...
     
     
     [+l(1),...
     -c(1, 2, 2, 2) - d(1, 2) * d(2, 2) + DijDkl(1,2,2,2),...
     -l(5) - 2 * c(1, 2, 2, 3) - 2 * d(1, 2) * d(2, 3) + 2 * DijDkl(1,2,2,3),...
     -c(1, 2, 2, 2) - d(1, 2) * d(2, 2) + DijDkl(1,2,2,2),...
     -c(2, 2, 2, 2) - d(2, 2)^2 + DijDkl(2,2,2,2),...
     -c(2, 2, 2, 3) - d(2, 2) * d(2, 3) + DijDkl(2,2,2,3),...
     -l(7) - 2 * c(1, 2, 2, 3) - 2 * d(1, 2) * d(2, 3) + 2 * DijDkl(1,2,2,3),...
     -c(2, 2, 2, 3) - d(2, 2) * d(2, 3) + DijDkl(2,2,2,3),...
     +l(9)];...
     
     
     [+l(2),...
     +l(5),...
     -c(1, 2, 3, 3) - d(1, 2) * d(3, 3) + DijDkl(1,2,3,3),...
     -c(1, 3, 2, 2) - d(1, 3) * d(2, 2) + DijDkl(1,3,2,2),...
     -c(2, 2, 2, 3) - d(2, 2) * d(2, 3) + DijDkl(2,2,2,3),...
     -c(2, 2, 3, 3) - d(2, 2) * d(3, 3) + DijDkl(2,2,3,3),...
     -l(8) - 2 * c(1, 3, 2, 3) - 2 * d(1, 3) * d( 2, 3) + 2 * DijDkl(1,3,2,3),...
     -l(9) - 2 * c(2, 3, 2, 3) - 2 * d(2, 3)^2 + 2 * DijDkl(2,3,2,3),...
     -c(2, 3, 3, 3) - d(2, 3) * d(3, 3) + DijDkl(2,3,3,3)];...
     
     
     [-c(1, 1, 1, 3) - d(1, 1) * d(1, 3) + DijDkl(1,1,1,3),...
     -l(3) - 2 * c(1, 2, 1, 3) - 2 * d(1, 2) * d(1, 3) + 2 * DijDkl(1,2,1,3),...
     -l(4) - 2 * c(1, 3, 1, 3) - 2 * d(1, 3)^2 + 2 * DijDkl(1,3,1,3),...
     -c(1, 1, 2, 3) - d(1, 1) * d(2, 3) + DijDkl(1,1,2,3),...
     -l(7) - 2 * c(1, 2, 2, 3) - 2 * d(1, 2) * d(2, 3) + 2 * DijDkl(1,2,2,3),...
     -l(8) - 2 * c(1, 3, 2, 3) - 2 * d(1, 3) * d(2, 3) + 2 * DijDkl(1,3,2,3),...
     -c(1, 1, 3, 3) - d(1, 1) * d(3, 3) + DijDkl(1,1,3,3),...
     -c(1, 2, 3, 3) - d(1, 2) * d(3, 3) + DijDkl(1,2,3,3),...
     -c(1, 3, 3, 3) - d(1, 3) * d(3, 3) + DijDkl(1,3,3,3)];...
     
     
     [+l(3),...
     -c(1, 3, 2, 2) - d(1, 3) * d(2, 2) + DijDkl(1,3,2,2),...
     -l(6) - 2 * c(1, 3, 2, 3) - 2 * d(1, 3) * d(2, 3) + 2 * DijDkl(1,3,2,3),...
     +l(7),...
     -c(2, 2, 2, 3) - d(2, 2) * d(2, 3) + DijDkl(2,2,2,3),...
     -l(9) - 2 * c(2, 3, 2, 3) - 2 * d(2, 3)^2 + 2 * DijDkl(2,3,2,3),...
     -c(1, 2, 3, 3) - d(1, 2) * d(3, 3) + DijDkl(1,2,3,3),...
     -c(2, 2, 3, 3) - d(2, 2) * d(3, 3) + DijDkl(2,2,3,3),...
     -c(2, 3, 3, 3) - d(2, 3) * d(3, 3) + DijDkl(2,3,3,3)];...
     
     [+l(4),...
     +l(6),...
     -c(1, 3, 3, 3) - d(1, 3) * d(3, 3) + DijDkl(1,3,3,3),...
     +l(8),...
     +l(9),...
     -c(2, 3, 3, 3) - d(2, 3) * d(3, 3) + DijDkl(2,3,3,3),...
     -c(1, 3, 3, 3) - d(1, 3) * d(3, 3) + DijDkl(1,3,3,3),...
     -c(2, 3, 3, 3) - d(2, 3) * d(3, 3) + DijDkl(2,3,3,3),...
     -c(3, 3, 3, 3) - d(3, 3)^2 + DijDkl(3,3,3,3)]...
     
     ];
 
 
end
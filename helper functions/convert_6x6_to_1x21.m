function out = convert_6x6_to_1x21(in)
% function convert_6x6_to_1x21
%
% converts a tensor from a 6x6 to a 1x21 representation
% picks up the upper triangular elements to be used in the Cholesky
% formulation
%
% in_6x6 = 
%           [       xx*xx       xx*yy         xx*zz         sqrt(2)*xx*xy sqrt(2)*xx*xz sqrt(2)*xx*yz
%                   xx*yy       yy*yy         yy*zz         sqrt(2)*yy*yx sqrt(2)*yy*xz sqrt(2)*yy*yz
%                   xx*zz       yy*zz         zz*zz         sqrt(2)*zz*xy sqrt(2)*zz*zx sqrt(2)*zz*zy
%              sqrt(2)*xx*xy sqrt(2)*yy*yx sqrt(2)*zz*xy        2*xy*xy      2*xy*xz       2*xy*yz
%              sqrt(2)*xx*xz sqrt(2)*yy*xz sqrt(2)*zz*zx        2*xy*xz      2*xz*xz       2*xz*yz
%              sqrt(2)*xx*yz sqrt(2)*yy*yz sqrt(2)*zz*zy        2*xy*yz      2*xz*yz       2*yz*yz
%           ]
% out_1x21 = [
%            [xx*xx yy*yy zz*zz] 
%            [xx*yy xx*zz yy*zz]                                               * sqrt(2)
%            [xx*yz yy*xz zz*xy xx*xy xx*xz yy*xy yy*yz zz*xz zz*yz] * sqrt(2) * sqrt(2)   
%            [xy*xy xz*xz yz*yz]                                     *    2  
%            [xy*xz xy*yz xz*yz]                                     *    2    * sqrt(2) 
%            ]
%

xxxx = in(1,1);
yyyy = in(2,2);
zzzz = in(3,3);

xxyy = in(1,2);
xxzz = in(1,3);
yyzz = in(2,3);
xxyz = in(1,6);
yyxz = in(2,5);
zzxy = in(3,4); 
xxxy = in(1,4);
xxxz = in(1,5);
yyxy = in(2,4);
yyyz = in(2,6);
zzxz = in(3,5);
zzyz = in(3,6);
xyxy = in(4,4);
xzxz = in(5,5);
yzyz = in(6,6);
xyxz = in(4,5);
xyyz = in(4,6);
xzyz = in(5,6);

f = sqrt(2);

out = [ ...
       [xxxx yyyy zzzz] ...
       [xxyy xxzz yyzz] * f ...
       [xxyz yyxz zzxy xxxy xxxz yyxy yyyz zzxz zzyz] * f ...
       [xyxy xzxz yzyz] ...
       [xyxz xyyz xzyz] * f ...
      ];

end
function out = convert_1x21_to_6x6(in)
% function convert_1x21_to_6x6()
%
% converts a tensor from a 1x21 to a 6x6 representation
%
% in_1x21 = [
%            [xx*xx yy*yy zz*zz] 
%            [xx*yy xx*zz yy*zz]                                               * sqrt(2)
%            [xx*yz yy*xz zz*xy xx*xy xx*xz yy*xy yy*yz zz*xz zz*yz] * sqrt(2) * sqrt(2)   
%            [xy*xy xz*xz yz*yz]                                     *    2  
%            [xy*xz xy*yz xz*yz]                                     *    2    * sqrt(2) 
%            ]
% out_6x6 = 
%           [       xx*xx       xx*yy         xx*zz         sqrt(2)*xx*xy sqrt(2)*xx*xz sqrt(2)*xx*yz
%                   xx*yy       xx*xx         yy*zz         sqrt(2)*yy*yx sqrt(2)*yy*xz sqrt(2)*yy*yz
%                   xx*zz       yy*zz         zz*zz         sqrt(2)*zz*xy sqrt(2)*zz*zx sqrt(2)*zz*zy
%              sqrt(2)*xx*xy sqrt(2)*yy*yx sqrt(2)*zz*xy        2*xy*xy      2*xy*xz       2*xy*yz
%              sqrt(2)*xx*xz sqrt(2)*yy*xz sqrt(2)*zz*zx        2*xy*xz      2*xz*xz       2*xz*yz
%              sqrt(2)*xx*yz sqrt(2)*yy*yz sqrt(2)*zz*zy        2*xy*yz      2*xz*yz       2*yz*yz
%           ]
%



if size(in,1) == 21
    in = in';
end

xxxx = in(1); 
yyyy = in(2);
zzzz = in(3);
xxyy = in(4); 
xxzz = in(5);
yyzz = in(6);
xxyz = in(7); 
yyxz = in(8);
zzxy = in(9);
xxxy = in(10);
xxxz = in(11);
yyxy = in(12);
yyyz = in(13);
zzxz = in(14);
zzyz = in(15);
xyxy = in(16);
xzxz = in(17);
yzyz = in(18);
xyxz = in(19);
xyyz = in(20);
xzyz = in(21);

f = 1 / sqrt(2);

out = [ ...
                 xxxx f * xxyy f * xxzz     f * [xxxy xxxz xxyz]   ;...
             f * xxyy     yyyy f * yyzz     f * [yyxy yyxz yyyz]   ;...
             f * xxzz f * yyzz     zzzz     f * [zzxy zzxz zzyz]   ;...
                 
            f * [xxxy yyxy zzxy]                 xyxy  f * xyxz  f * xyyz;...
            f * [xxxz yyxz zzxz]             f * xyxz      xzxz  f * xzyz;...
            f * [xxyz yyyz zzyz]             f * xyyz  f * xzyz      yzyz ...
     ];
end
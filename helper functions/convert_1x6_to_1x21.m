function out = convert_1x6_to_1x21(in)
% function convert_1x6_to_1x21()
%
% converts a tensor from a 1x6 to a 1x21 representation
%
% in_1x6   = [xx, yy, zz, sqrt(2)xy, sqrt(2)xz, sqrt(2)yz]
%
% tmp_6x6 = 
%           [       xx*xx       xx*yy         xx*zz         sqrt(2)*xx*xy sqrt(2)*xx*xz sqrt(2)*xx*yz
%                   xx*yy       xx*xx         yy*zz         sqrt(2)*yy*yx sqrt(2)*yy*xz sqrt(2)*yy*yz
%                   xx*zz       yy*zz         zz*zz         sqrt(2)*zz*xy sqrt(2)*zz*zx sqrt(2)*zz*zy
%              sqrt(2)*xx*xy sqrt(2)*xx*xz sqrt(2)*xx*yz        2*xy*xy      2*xy*xz       2*xy*yz
%              sqrt(2)*yy*yx sqrt(2)*yy*xz sqrt(2)*yy*yz        2*xy*xz      2*xz*xz       2*xz*yz
%              sqrt(2)*zz*xy sqrt(2)*zz*zx sqrt(2)*zz*zy        2*xy*yz      2*xz*yz       2*yz*yz
%           ]
%
% eqn (48) from Westin 16, different order
% out_1x21 = [
%            [xx*xx yy*yy zz*zz] 
%            [xx*yy xx*zz yy*zz]                                               * sqrt(2)
%            [xx*yz yy*xz zz*xy xx*xy xx*xz yy*xy yy*yz zz*xz zz*yz] * sqrt(2) * sqrt(2)   
%            [xy*xy xz*xz yz*yz]                                     * sqrt(4)  
%            [xy*xz xy*yz xz*yz]                                     * sqrt(4) * sqrt(2) 
%            ]
%             


if size(in,2) ~= 6
    in = in';
end


xx = in(:,1);
yy = in(:,2);
zz = in(:,3);
xy = in(:,4);
xz = in(:,5);
yz = in(:,6);

f = sqrt(2);

out = [ ...
        [xx.*xx yy.*yy zz.*zz]... 
        [xx.*yy xx.*zz yy.*zz]                                           * f ...
        [xx.*yz yy.*xz zz.*yz xx.*xy xx.*xz yy.*xy yy.*yz zz.*xz zz.*yz] * f ...   
        [xy.*xy xz.*xz yz.*yz]...                                       
        [xy.*xz xy.*yz xz.*yz]                                           * f ... 
       ];  
function bt = get_btensors_from_bvec_bval(bvec, bval, bdelta)
% function get_btensor_from_bvec_bval(bvec, bval, bdelta)
%
% returns the b-tensors (in Voigt format) associated to the bvec, bval, 
% and bdelta. If bdelta is not provided, it is assumed to be 1 
% (LTE or Stejskal-Tanner experiments only)
%
% the implementation follows that in Dipy

nmeas = length(bval);
if nargin < 3
    bdelta = ones(nmeas,1);
end

% define the tensors' shape
linear = [1 0 0; 0 0 0; 0 0 0];
planar = [0 0 0; 0 1 0; 0 0 1] / 2 ;
spherical = eye(3) / 3;

bt = zeros(nmeas,6);
for i = 1 : nmeas
    if bdelta(i) == 0
        
        tmp = spherical * bval(i) * 1e6;
        
    else if bdelta(i) == 1
            R = vec2rotmat(bvec(:,i));
            tmp = R * linear * R' * bval(i) * 1e6;
        else
            R = vec2rotmat(bvec(:,i));
            tmp = R * planar * R' * bval(i) * 1e6;
        end
    end
    
    bt(i,:) = convert_3x3_to_1x6(tmp);
    
end
end

function out = vec2rotmat(bvec)
% returns a rotation matrix that aligns [1 0 0] (aka the symmetry axis of 
% the defined shapes) to bvec
%
% check this for algebra
% https://math.stackexchange.com/questions/180418/
% calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

if size(bvec,1) == 3
    bvec = bvec';
end
if norm(bvec) == 0
    out = eye(3);
    return
end

w = cross([1 0 0], bvec);
wn = norm(w, 'fro');
w = w / wn;

vp = (bvec - (dot([1 0 0], bvec)) * [1 0 0]);
vp = vp / norm(vp, 'fro');

P = vertcat([1 0 0], vp, w);
Pt = P';
cosa = dot([1 0 0], bvec);
sina = sqrt(1 - cosa^2);
R = [cosa -sina 0;...
    sina   cosa 0;...
    0      0  1];
out = Pt * (R * P);

end

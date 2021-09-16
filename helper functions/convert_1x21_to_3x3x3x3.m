function out = convert_1x21_to_3x3x3x3(in)
% function out = convert_1x21_to_3x3x3x3_cvx()
%

if size(in,1) == 21
    in = in';
end

% here we remove all scaling factor used for the voigt format
xxxx = in(1); 
yyyy = in(2);
zzzz = in(3);

f = sqrt(2);
xxyy = in(4) / f; 
xxzz = in(5) / f;
yyzz = in(6) / f;

f = sqrt(4);
xxyz = in(7) / f; 
yyxz = in(8) / f;
zzxy = in(9) / f;
xxxy = in(10) / f;
xxxz = in(11) / f;
yyxy = in(12) / f;
yyyz = in(13) / f;
zzxz = in(14) / f;
zzyz = in(15) / f;
xyxy = in(16) / f;
xzxz = in(17) / f;
yzyz = in(18) / f;

f = sqrt(8);
xyxz = in(19) / f;
xyyz = in(20) / f;
xzyz = in(21) / f;

%% compose the 3x3x3x3 tensor using its symmetries
% xxxx, yyyy, zzzz
out(1,1,1,1) = xxxx;
out(2,2,2,2) = yyyy;
out(3,3,3,3) = zzzz;

% xxyy, xxzz, yyzz
out(1,1,2,2) = xxyy;
out(2,2,1,1) = out(1,1,2,2);
out(1,1,3,3) = xxzz;
out(3,3,1,1) = out(1,1,3,3);
out(2,2,3,3) = yyzz;
out(3,3,2,2) = out(2,2,3,3);

% xxyz
out(1,1,2,3) = xxyz;
out(1,1,3,2) = out(1,1,2,3);
out(2,3,1,1) = out(1,1,2,3);
out(3,2,1,1) = out(1,1,2,3);

% yyxz
out(2,2,1,3) = yyxz;
out(2,2,3,1) = out(2,2,1,3);
out(1,3,2,2) = out(2,2,1,3);
out(3,1,2,2) = out(2,2,1,3);

% zzxy
out(3,3,1,2) = zzxy;
out(3,3,2,1) = out(3,3,1,2);
out(1,2,3,3) = out(3,3,1,2);
out(2,1,3,3) = out(3,3,1,2);

% xxxy
out(1,1,1,2) = xxxy;
out(1,1,2,1) = out(1,1,1,2);
out(1,2,1,1) = out(1,1,1,2);
out(2,1,1,1) = out(1,1,1,2);

% xxxz
out(1,1,1,3) = xxxz;
out(1,1,3,1) = out(1,1,1,3);
out(1,3,1,1) = out(1,1,1,3);
out(3,1,1,1) = out(1,1,1,3);

% yyxy
out(2,2,1,2) = yyxy;
out(2,2,2,1) = out(2,2,1,2);
out(1,2,2,2) = out(2,2,1,2);
out(2,1,2,2) = out(2,2,1,2);

% yyyz
out(2,2,2,3) = yyyz;
out(2,2,3,2) = out(2,2,2,3);
out(2,3,2,2) = out(2,2,2,3);
out(3,2,2,2) = out(2,2,2,3);

% zzxz
out(3,3,1,3) = zzxz;
out(3,3,3,1) = out(3,3,1,3);
out(1,3,3,3) = out(3,3,1,3);
out(3,1,3,3) = out(3,3,1,3);

% zzyz
out(3,3,2,3) = zzyz;
out(3,3,3,2) = out(3,3,2,3);
out(2,3,3,3) = out(3,3,2,3);
out(3,2,3,3) = out(3,3,2,3);

% xyxy
out(1,2,1,2) = xyxy;
out(1,2,2,1) = out(1,2,1,2);
out(2,1,1,2) = out(1,2,1,2);
out(2,1,2,1) = out(1,2,1,2);


%xzxz
out(1,3,1,3) = xzxz;
out(1,3,3,1) = out(1,3,1,3);
out(3,1,1,3) = out(1,3,1,3);
out(3,1,3,1) = out(1,3,1,3);

% yzyz
out(2,3,2,3) = yzyz;
out(2,3,3,2) = out(2,3,2,3);
out(3,2,2,3) = out(2,3,2,3);
out(3,2,3,2) = out(2,3,2,3);

% xyxz
out(1,2,1,3) = xyxz;
out(1,2,3,1) = out(1,2,1,3);
out(2,1,1,3) = out(1,2,1,3);
out(2,1,3,1) = out(1,2,1,3);
out(1,3,1,2) = out(1,2,1,3);
out(1,3,2,1) = out(1,2,1,3);
out(3,1,1,2) = out(1,2,1,3);
out(3,1,2,1) = out(1,2,1,3);

% xyyz
out(1,2,2,3) = xyyz;
out(1,2,3,2) = out(1,2,2,3);
out(2,1,2,3) = out(1,2,2,3);
out(2,1,3,2) = out(1,2,2,3);
out(2,3,1,2) = out(1,2,2,3);
out(2,3,2,1) = out(1,2,2,3);
out(3,2,1,2) = out(1,2,2,3);
out(3,2,2,1) = out(1,2,2,3);

% xzyz
out(1,3,2,3) = xzyz;
out(1,3,3,2) = out(1,3,2,3);
out(3,1,2,3) = out(1,3,2,3);
out(3,1,3,2) = out(1,3,2,3);
out(2,3,1,3) = out(1,3,2,3);
out(2,3,3,1) = out(1,3,2,3);
out(3,2,1,3) = out(1,3,2,3);
out(3,2,3,1) = out(1,3,2,3);

end
function out = convert_1x6_to_3x3(in)
if size(in,1) == 6
    in = in';
end

f = 1/sqrt(2);

out = [  in(1)    f * in(4)  f * in(5); ...
       f * in(4)    in(2)    f * in(6);...
       f * in(5)  f * in(6)    in(3)  ];

end
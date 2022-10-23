function dot = inner_product(s1,s2)
%In C, the user must specify the length of the sequences.
[r1,c1] = size(s1);
[r2,c2] = size(s2);
el1 = r1*c1;
el2 = r2*c2;
if el1 ~= el2
    error('Input sequences must be of the same length.')
end
dot = 0; %Initialize buffer.
for ii = 1:el1 %Append until it is done.
    dot = dot + s1(ii)*s2(ii);
end
end
function val = horner(coeff,x)
[rc,cc] = size(coeff);
if rc ~= 1 && cc ~= 1 %
    error('Must input a column/row array for the coefficients!');
end
n_coeff = length(coeff);
val = ones(size(x))*coeff(n_coeff); %Vector initialized to last coefficient.
if n_coeff == 1 %Input is a scalar.
    return;
end
val = val.*x;
for ii = (n_coeff - 1):-1:2
    val = (val + coeff(ii)).*x;
end
val = val + coeff(1);
end
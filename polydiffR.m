%syms x
%f = (1 + x + x^2)/(2 + 6*x + x^2 + 5*x^3);
%df = diff(f,x);
%ddf = diff(df,x);
%dddf = diff(ddf,x);
%G = [1,1,1];
%H = [2,6,1,5];
%d = 3;
%val = 2;
%dF =  PolyderR(G,H,d,val);
function dF = polydiffR(G,H,d,val)
%Order of numerator polynomial
[rG,cG] = size(G); %In C, the user must provide the size.
if rG ~= 1 && cG ~= 1
    error('Input numerator coefficients either as a row or column array.');
end
if rG > cG 
    pG = rG - 1;
else
    pG = cG - 1;
end
%Order of denominator polynomial
[rH,cH] = size(H); %In C, the user must provide the size.
if rH ~= 1 && cH ~= 1
    error('Input numerator coefficients either as a row or column array.');
end
if rH > cH
    pH = rH - 1; %Order of numerator polynomial.
else
    pH = cH - 1;
end

%Memory allocation and initialization.
dHb = zeros(1 + pH,1); %Polynomial coefficient buffers for denominator.
for ii = 1:(1 + pH)
    dHb(ii) = H(ii);
end
dGb = zeros(1 + pG,1); %Polynomial coefficient buffers for numerator.
for ii = 1:(1 + pG)
    dGb(ii) = G(ii);
end
dH = zeros(1+d,1); %Stores intermediate values of h^(k).
dF = zeros(1+d,1); %Stores output values of F^(k).
% Computation
dH(1) = horner(H,val);
if dH(1) == 0
    error('Input denominator polynomial evaluates to zero!');
end
dF(1) = horner(G,val)/dH(1);

for ii = 2:(1 + d) %First derivative occurs with ii = 2.
    dGb = polydiff(dGb,1,true); %Compute next numerator derivative inplace.
    dHb = polydiff(dHb,1,true); %Compute next denominator derivative inplace.
    dG = horner(dGb,val); %Computes g^(ii-1)
    dH(ii) = horner(dHb,val); %Computes h^(ii-1)
    
    dF(ii) = dG;
    for jj = 2:ii
        dF(ii) = dF(ii) -  nCr(ii-1,jj-1)*dH(jj)*dF(ii-jj+1);
    end
    dF(ii) = dF(ii)/dH(1);
end
end
function c = nCr(n,r)
if n <= 0
    error('n must be positive!');
end
if r < 0
    error('r cannot be negative!');
end
if r > n
    error('r must be less than n!')
end
c = 1;
if r == n || r == 0
    return;
end
%Inspect the two factors of the denominator of the nCr formula. Determine
%which is the bigger.
if n-r > r
    small = r;
    big = n-r;
else
    small = n-r;
    big = r;
end
%Use the "big" factor to compute the numerator using the least
%multiplications.
for ii = big + 1:n
    c = c*ii;
end
%Divide
for ii = 1:small
    c = c/ii;
end
%  1 2 3 4 5 6 7 8 9
% (1 2 3 4) (1 2 3 4 5)
end
function dim = mustbeflat(array)
[r,c] = size(array);
if r ~= 1 && c ~= 1 || r == c
    error('Input array must be flat (row or column)');
elseif r > c
    dim = r;
elseif c > r
    dim = c;
end
end
function [array,r,c] = stand(array,side)
if nargin < 2
    side = 'tall';
end
[r,c] = size(array);
if strcmp(side,'tall') == 1
    %If less rows than columns, array is standing on its short side.
    if r < c %Transpose to tall side.
        array = array';
        buffer = r;
        r = c; %"c" is the new # of rows.
        c = buffer;
    end
elseif strcmp(side,'short') == 1
    %If less columns than rows, array is standing on its long side.
    if c < r %Transpose to short side.
        array = array';
    end
    buffer = c;
    c = r; %"r" is the new # of columns.
    r = buffer;
else
    error('Invalid side input!');
end
end
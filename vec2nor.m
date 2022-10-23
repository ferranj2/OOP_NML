function [un,uv] = vec2nor(vector)
mag = sqrt(inner_product(vector,vector));
uv = vector/mag; %Produce unit vector.
un = [...
    uv(1),uv(2);... %Dot product between the two is zero.
    uv(2),-uv(1)]... %Cross product returns i = j = 0, k = +1.
    \[0;1];
end
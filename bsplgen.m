function [B,k_range,k_idx] = bsplgen(kn,p)
if p < 1 
    error('This routine is for order 1 or above.')
end
[rk,ck] = size(kn);
if rk ~= 1 && ck ~= 1
    error('Must input a row or column array for the knots.');
end
if rk == 1 %Knots input as row vector.
    knots = ck;
end
if ck == 1 %Knots input as column vector.
    knots = rk;
end
if knots < p+2
    error(['Need "p+2" =',num2str(p+2),'knots. Only',num2str(knots),'input']);
end
if knots > p+2
    warning(['Only "p+2" =',num2str(p+2),'knots needed. Detected',num2str(knots-p-2),' in excess. Ignoring.']);
end
for ii = 1:knots-1 %Check if the knot sequence is non-decreasing.
    if kn(ii+1) - kn(ii) < 0
        error('Knot sequence must be non-decreasing.')
    end
end
%Allocate buffer arrays.
B = cell(p+1,1); %Need pointers to buffer arrays.
for ii = 1:p+1
    B{ii} = zeros(p+2-ii,p+2-ii);
    B{ii}(1,1) = 1; %This one is th B^0 spline.
end
%Generate the splines.
for kk = 1:p %Counter "kk" denotes the intermediate order of B-spline being made.
    idx = 1:kk;
    for jj = 1:p-kk+1 %Compute the "k^th" order B-splines into the respective buffers.
        den = kn(jj+kk) - kn(jj); %Denominator
        if den ~= 0 %"Bootstrap" component.
            B{jj}(idx,idx) = -B{jj}(idx,idx)*kn(jj)/den;
            B{jj}(idx,idx+1) = B{jj}(idx,idx+1) - B{jj}(idx,idx)/kn(jj);
        end
        den = kn(jj+kk+1) - kn(jj+1); %Denominator
        if den ~= 0 %"Next" spline component.
            B{jj}(idx+1,idx) = B{jj}(idx+1,idx) + B{jj+1}(idx,idx)*kn(jj+1+kk)/den;
            B{jj}(idx+1,idx+1) = B{jj}(idx+1,idx+1) - B{jj+1}(idx,idx)/den;
        end
    end
end
if nargout > 1
    k_range = zeros(knots-1,2);
    for kk = 1:knots-1
        k_range(kk,1) = kn(kk);
        k_range(kk,2) = kn(kk+1);
    end
end
if nargout > 2
    k_idx = zeros(knots-1,2);
    for kk = 1:knots-1
        k_idx(kk,1) = kk;
        k_idx(kk,2) = kk + 1;
    end
end
end
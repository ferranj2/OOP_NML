function B = bsplregen(B,kn)
buffers = length(B);
%[rB1,cB1] = size(B{1});
%p = cB1 - 1;
p = buffers -1 ;
%These copies the reusable intermediate splines to one buffer ahead.
for ii = 1:(buffers - 1) %Only up the B^{1}_{0} buffer is copied.
    for jj = 1:(p + 1 - ii)
        for kk = 1:(p + 1 - ii)
            B{ii}(jj,kk) = B{ii + 1}(jj,kk);
        end
    end
    %Zero excess column
    for jj = 1:p+2-ii
        B{ii}(jj,p+2-ii) = 0;
    end
    %Zero excess row.
    for jj = 1:p+1-ii %One common zero already taken care off.
        B{ii}(p+2-ii,jj) = 0;
    end
end

%Use the Cox-DeBoor formula to generate new splines.
for pp = 1:p %Order of the intermediate spline.
    idx = 1:pp; %CHEATING
    jj = p + 1 - pp; %Target Buffer.
    den = kn(jj + pp) - kn(jj); %Denominator
    if den ~= 0 %"Bootstrap" component.
        B{jj}(idx,idx) = -B{jj}(idx,idx)*kn(jj)/den;
        B{jj}(idx,idx+1) = B{jj}(idx,idx+1) - B{jj}(idx,idx)/kn(jj);
    end
    den = kn(jj+pp+1) - kn(jj+1); %Denominator
    if den ~= 0 %"Next" spline component.
        B{jj}(idx+1,idx) = B{jj}(idx+1,idx) + B{jj+1}(idx,idx)*kn(jj+1+pp)/den;
        B{jj}(idx+1,idx+1) = B{jj}(idx+1,idx+1) - B{jj+1}(idx,idx)/den;
    end
end
end
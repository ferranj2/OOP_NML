function O = polydiff(I,der,inplace)
if der == 0
    warning('der = 0 was input. Nothing done.');
    O = I;
    return;
end
if der < 1
    error('Derivative must be at least 1.');
end
if nargin < 3
    inplace = false;
end
switch inplace
    case true %Reuse buffer and shift elements accordingly.
        if der > length(I) - 1
            for ii = 1:length(I)
                I(ii) = 0;
            end
            O = I;
            return;
        end
        for ii = 1:der
            for jj = 1:(length(I) - ii) %This pushes latter coefficients to the left.
                I(jj) = I(jj+1)*jj; 
            end
            for jj = (length(I) - ii + 1):length(I)
                I(jj) = 0; %As derivatives are taken, latter entries get zeroed.
            end
        end
        O = I;
    case false %Create a new custom-sized Buffer.
        if der > length(I) - 1
            O = 0;
            return;
        end
        O = zeros(length(I)-der,1);
        for ii = 1:length(I)-der
            O(ii) = I(ii+der); %Copy the coefficients respective coefficients.
        end
        mult = der; %Initial multiplier.
        for ii = 1:der
            for jj  = 1:length(O)
                O(jj) = O(jj)*(mult + jj - 1);
            end
            mult = mult -1; %As powers are lowered, so is the multiplier.
        end
    otherwise
        error('Invalid input for the "inplace" flag.');
end
end
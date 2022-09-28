%% FUNCTION: Finite Difference Scheme Deducer
% * INPUT(*P*): Number of points around the "ith" point where
% derivative is to be approximated. The scheme will formulate an expression
% in terms of "i + points."
% * INPUT(*D*): If the number of points input is odd, this
% intrinsically corresponds to either a forward or backward difference.
% The "direction" of the finite difference must, therefore, be specified
% when inputting an odd number of points. This input is optional when the
% input points are even.
% * INPUT(*S*): This allows users to derive finite difference schemes with
% a "shift". Purely forward or purely backward differences make no
% reference to points before or after the "ith" respectively. This inputs
% shifts "S" points forward "+S" or backward "-S" to derive mostly forward
% mostly backward schemes that use "S" points before or after the "ith"
% one. This input is only valid when "D" is either "+1" or "-1".
% * INPUT(*dx*): Uneven spacings the coefficients that are usually returned
% by the routine are implicitly in terms of a constant spacing between the
% nodes. If the user needs schemes that account for uneven spacing, they
% may input them. The output coefficients will, this time have no implict
% spacing.
% * OUTPUT(*coeffs*): Coefficients of the nodes.
% * OUTPUT(*dfdx*): Order od the derivative.
% * OUTPUT(*order*): Order of the truncation error (meaningful only in even spacing).
%%
% $f(x) = f(a) + \frac{df}{}$
%%
% $f(x + dx) = f(x) + \frac{df}{dx} dx + \frac{d^{2}f}{dx^{2}}\frac{dx^{2}}{2!} + ... + \frac{d^{n}f}{dx^{n}}\frac{dx^{n}}{n!}$
function [coeffs,dfdx,order] = FDC(varargin)
if mod(nargin,2) == 1 %Incomplete input pair.
    error('Uneven number of inputs! Must supply complete pairs.')
elseif nargin > 8 %Absolutely too many inputs.
    error('Excessive number of input pairs! This function deals with only five(5).')
elseif nargin < 6 %Absolutely too few inputs
    error('This functions needs atleast 3 inputs!')
end
P = 2; %Need at least two points for any scheme.
D = 1; %Forward difference is the default.
S = 0; %No stencil shifting by default.
dx = []; %Even spacing by default.
for ii = 2:2:nargin %For each argument pair.
    switch varargin{ii-1}
        case 'P'
            if isnumeric(varargin{ii}) == 0
                error('INVALID TYPE. Stencil points must be numeric.');
            end
            if varargin{ii} < 2
                error('INVALID INPUT. Need at least 2 stencil points!');
            end
            P = varargin{ii};
        case 'D'
            if isnumeric(varargin{ii}) == 0
                error('INVALID TYPE. Direction of the differentiation must be numeric!');
            end
            if varargin{ii} ~= 1 && varargin{ii} ~= -1 && varargin{ii} ~= 0
                error('INVALID INPUT. Direction must be either forward "+1", backward "-1", or central "0"!');
            end
            D = varargin{ii};
        case 'S'
            if isnumeric(varargin{ii}) == 0
                error('INVALID TYPE. Stencil shift must be numeric!');
            end
            S = varargin{ii};
        case 'dx'
            if isnumeric(varargin{ii}) == 0
                error('INVALID TYPE. Stencil shift must be numeric!');
            end
            dx = varargin{ii};
        otherwise
            error('INVALID INPUT: Cannot recognize KEY name.')
    end
end
% Error checking
if D == 0 && mod(P,2) == 1
    error('ERROR: Central difference schemes require even stencil points!')
end
if D == 0 && abs(S) > P/2
    error('ERROR: Stencil shift abandons "ith" point.')
end
if D == 1 && S > 0 || S < -P
    error('ERROR: Stencil shift abandons "ith" point.')
end
if D == -1 && S < 0 || S > P
    error('ERROR: Stencil shift abandons "ith" point.')
end
% Computing the finite difference schema.
switch D
    case +1 %Forward difference.
        ith = 1;
        kk = 1;
    case 0 %Central difference.
        ith = P/2+1;
        kk = -P/2;% f_(i+k) = f_i + df/dx*(k*dx) + HOT.
    case -1 %backward difference.
        kk = -P;
        ith = P+1;
    otherwise
        error('Invalid input for D. Must be one of "+1", "0", or "-1" (numeric input, not string).')
end
if isempty(dx) == 0
    dx = dx(:); %Flatten the array of uneven spacings.
    if length(dx) > P
        error('ERROR: More spacings that stencil points input!')
    end
else
    if D == 0
        dx = [-P/2:+1:-1,1:+1:P/2];
    elseif D == -1
        dx = -P:+1:-1;
    elseif D == +1
        dx = +1:+1:+P;
    end
end

%kk is the factor of "dx" as in:
% f_(i+2) = f_i + 2*dx*(df/dx) + 2*(d^2f/dx^2)/2!
% kk above is equal to 2.
M = zeros(P,P);%Allocate memory for coefficient matrix.
coeffs = zeros(P,P+1);
dfdx = zeros(P,1); %Output indices of the derivative that corresponds to the row of the output.
order = zeros(P,1); %Output the orders if the derivatives.
kk = kk + S; %Apply shift (if any).
for ii = 1:P
    dfdx(ii) = ii; %Enumerate the derivative.
    order(ii) = P - ii + 1; %Compute the order.
    for jj = 1:P
        M(ii,jj) = power(dx(ii),jj + S)/factorial(jj + S); %Optimize this, no need to call factorial.
    end
    kk = kk + 1;
    if kk == 0%Skip k corresponding to f_i.
        kk = 1;
    end
end
Minv = M^-1;
coeffs(:,ith) = -sum(Minv,2); %Sum row-wise
coeffs(:,[1:ith-1,ith+1:P+1]) = Minv;
end

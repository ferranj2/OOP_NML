%% polynomial.m
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% A data structure used to represent univariate polynomials. This is meant
% to act as a building block for a greater data structure used to represent
% splines. This early version supports elementary algebra and calculus
% operations such as simple root factoring or differentiation. Custom
% creation routines for special interpolation polynomials such as the
% Bernstein and Lagrange polynomial constituents are available.
%% Formulae
%% Special Polynomials
% Lagrange Polynomials
%%
% Bernstein Polynomials
%% Root-finding formulae and algorithms 
% Quadratic formula
%%
% $x_{1,2} = \frac{-b \pm \sqrt{b^{2} - 4ac}}{2a}$
%% Class definition
classdef polynomial < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        N %Number of points to evaluate the curve at.
        x1%Lower limit of the visualization range.
        x2%Upper limit of the visualization range.
        
        color %Color of polynomial's graphics as will be rendered in an axes object.
        name %Name of the polynomial's graphics as they appear on axes legend.
        canvas %Axes object on which to draw the polynomial.
        sketches %Structure containing handles to the object's graphics.
    end%properties (public)
    %DEFINING DATA variables
    properties (SetAccess = protected)
        degree %Highest power of the polynomial's indeterminate.
        coefficients %Constant multipliers of the polynomial.
    end%properties
    %METRIC variables
    properties (SetAccess = protected)
        maxima
        minima
        roots
    end% properties
    %FLAG and STATE variables.
    properties (Hidden = true)
        roots_computed %Whether the polynomial had its roots computed.
        valid %Whether this polynomial is valid.
        
        %Graphics related flags.
        canvas_set
        graphics_initialized
    end%properties
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = polynomial
            %Related to the object.
            this.degree = [];
            this.coefficients = [];
            
            %State variables
            this.roots_computed = false;
            this.valid = false;
            
            %Graphics related
            this.canvas_set = false;
            this.graphics_initialized = false;
            this.sketches = struct(...
                'Curve',[],...
                'Roots',[],...
                'Inflection',[],...
                'Extrema',[]);
        end%function
        
        %Custom Creation routines.
        function poly = CreateFromRoots(n_roots,roots)
            %Create from roots
            poly = polynomial; %Call constructor.
            poly.degree = n_roots; %Degree is defined by
            poly.roots = zeros(1,poly.degree);
            for ii = 1:poly.degree
                poly.roots(ii) = roots(ii);
            end%ii
            poly.coefficients = zeros(1,poly.degree + 1);
            poly.coefficients = polynomial.ExpandRoots(n_roots,roots,poly.coefficients);
            
            
        end%function
        function poly = CreateFromCoefficients(n_coeff,coeff)
            %Create a polynomial from input coefficients.
            if n_coeff < 1
                error('Must input atleast one coefficient.')
            end%if
            poly.degree = n_coeff - 1;
            poly.coefficients = coeff;
        end%function
        function poly = CreateFromDifferentiation(original)
            %Create a polynomial by taking the derivative of an already
            %existing polynomial.
            poly = polynomial;
            poly.degree = original.degree - 1;
            poly.coefficients = polynomial.Differentiate(...
                original.degree + 1,...
                original.coefficients);
        end%function
        function poly = CreateFromIntegration(original,int_const)
            %Create a polynomial by integrating an already existing polynomial.
            %The user should input an integration constant, else, a default of
            %zero is assumed.
            poly = polynomial;
            poly.degree = original.degree + 1;
            if nargin == 1
                int_const = 0;
            end%if
            poly.coefficients = zeros(1,poly.degree + 1);
            poly.coefficients(1) = int_const;
            poly.coefficients = polynomial.Integrate(...
                poly.degree,... %#of coefficients in the original coefficient array.
                original.coefficients,...
                int_const);%Integration constant.
        end%function
        function poly = CreateFromAddition(polynomial_1,polynomial_2)
            %Create a polynomial by adding two already existing ones.
            poly = polynomial;
            poly.degree = polynomial.greater(...
                polynomial_1.degree,...
                polynomial_2.degree);
            poly.coefficients = polynomial.Addition(...
                polynomial_1.degree + 1,... %#Coefficients in polynomial_1.
                polynomial_1.coefficients,...
                polynomial_2.degree + 1,... %#Coefficients in polynomial_2.
                polynomial_2.coefficients);
            poly.valid = true;
        end%function
        
    end%methods (Static)
    %High-level instance MODIFICATION and QUERY routines.
    methods
        function Print(this)
            fprintf('%f\t',this.coefficients(1));
            if this.degree >= 1
                fprintf('+%fx\t',this.coefficients(2));
            end%if
            for ii = 3:(this.degree + 1)
                fprintf('+%fx^%i\t',this.coefficients(ii),ii-1);
            end%ii
            fprintf('\n');
        end%function
        
        %Evaluate polynomial at 
        function Y_out = Evaluate(this,n_x,X_in)
            Y_out = zeros(n_x,1);%Allocate output buffer.
            n_coeff = this.degree + 1; %So that the computation does not happen at each loop iteration.
            for ii = 1:n_x
                Y_out(ii) = Horner(...
                    n_coeff,...
                    this.coefficients,...
                    X_in(ii));
            end%ii
        end%function
        
        %Displace this polynomial vertically by some amount.
        function VerticalShift(this,shift)
            this.coefficients(1) = this.coefficients(1) + shift;
        end%function
        
        %Displace this polynomial horizontally by some amount
        function HorizontalShift(this,shift)
            
        end%shift
        
        %Differentiate inplace
        function DifferentiateInplace(this)
            if this.degree == 0
                this.coefficients(1) = 0;
                return;
            end%if
            this.coefficients = polynomial.Differentiate(...
                this.degree + 1,... %#Coefficients.
                this.coefficients);
            this.degree = this.degree - 1;
        end%if
        %Integrate inplace
        function IntegrateInplace(this,int_const)
            if nargin == 1
                int_const = 0;
            end%if
            this.coefficients = polynomial.Integrate(...
                this.degree + 1,...
                this.coefficients,...
                int_const);            
            this.degree = this.degree + 1;
        end%function
        %Draw polynomial
        function DrawRange(this,ax,x1,x2,N)
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch) 
                this.sketch = polynomial.DrawPolynomialRange(...
                    ax,...
                    this.degree + 1,...
                    this.coefficients,N,x1,x2);
            else
                this.sketch.XData = linspace(x1,x2,N);
                this.sketch.YData = polynomial.EvaluateRange(...
                    this.degree + 1,...
                    this.coefficients,...
                    N,...
                    x1...
                    ,x2);
            end%if
        end%function
    end% methods (Ordinary)
    %Graphical setups.
    methods
        %Whether the polygon is to be drawn.
        function Show(this)
            %Make sure this object has a definition that is valid.
            if ~this.valid
                warning('The polygon object does not have a valid definition!');
                return;
            end%if
            %Check if the object has a designated canvas to be drawn on.
            if ~this.canvas_set %nope?
                warning('No canvas set prior to drawing query. New canvas created.');
                this.canvas = custom_axis; %Then make one.
                this.canvas_set = true;
            end%if
            %If graphics do not exist, create them.
            if ~this.graphics_initialized 
                InitializeGraphics(this);
            end%if
        end%function
        
        %Set the axes object onto which to draw the graphics.
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.graphics_initialized
                this.sketches.Curve.Parent = ax;
                this.sketches.Roots.Parent = ax;
                this.sketches.Extrema.Parent = ax;
                this.sketches.Inflection.Parent = ax;
            end%if
        end%function
        function SetName(this,string)
            this.name = string;
            if graphics.initialized
                this.sketches.Curve.DisplayName = string;
                this.sketches.Roots.DisplayName = [string,'Roots'];
            end%if
        end%function
        function SetColor(this,RGB)
            this.color = RGB;
            if graphics.initialized
                this.sketches.Curve.Color = RGB;
                this.sketches.Roots.MarkerFaceColor = RGB;
            end%if
        end%function
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',0,... %Do not initalize with empty ("[]" ) because...
                'YData',0,... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','-',...
                'Color',this.color,...
                'DisplayName',this.name);
            
            %Show the roots.
            this.sketches.Roots = line(...
                'Parent',this.canvas,...
                'XData',NaN(1,this.degree),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(1,this.degree),... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','none',...
                'MarkerFaceColor',this.color,...
                'DisplayName',[this.name,'Roots']);             
            this.graphics_initialized  = true;
        end%function.           
        function TerminateGraphics(this)
            %Release system resources used to render graphics.
            delete(this.sketches.Curve);
            delete(this.sketches.Roots);
            
            %Update the flag.
            this.graphics_initialized  = false;
        end%function
    end%methods (Graphics)
    %Graphical demonstrations
    methods
    end %methods (Demonstrations)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Root finding formulae and algorithms
        function [xr1,xr2] = QuadraticFormula(a,b,c)
            xr1 = -0.5*b/a;
            xr2 = xr1 - sqrt(xr1*xr1 - c/a);
            xr1 = 2*xr1 - xr2;
        end%function
        function [xr1,xr2,xr3] = CubicFormula(a0,a1,a2,a3)
            p = a1 - a2*a2/3;
            q = a1*a2/3 - a0 -2*a2*a2*a2/27;
        end%function
        
        %Custom operations involving roots of the polynomial.
        function coeff = ExpandRoots(n_roots,roots,coeff)
            if isempty(coeff) %Is this bad practice in C and in general?
                coeff = zeros(1,n_roots + 1);%Allocate output memory for polynomial coefficients.
            end%if
            
            %Initialize the coefficient buffer.
            coeff(1) = -roots(1);
            coeff(2) = 1;
            for ii = 2:n_roots
                iip1 = ii + 1;
                %Shift coefficients to the right. This is akin to
                %multiplying by "x."
                for jj = iip1:-1:2
                    coeff(jj) = coeff(jj-1);
                end%jj
                coeff(1) = 0; %Need to zero the first coefficient.
                
                %Now, multiply by the "negative root" and add to the lower
                %powers of the indeterminates.
                for jj = 1:ii
                    coeff(jj) = coeff(jj) - roots(ii)*coeff(jj + 1);
                end%jj
            end%ii
            
        end%function
        function coeff = InPlaceRootDeletion(degree,coeff,root)
            %A specialized division routine for removing roots of
            %polynomials (if known). If "coeff" contains the coefficients
            %of polynomial "P", and "root" is a know root of "P", then this
            %routine carries out the division of P/(x-root).
            
            %WARNING: This routine cannot remove roots at x = 0.
            if root ~= 0
                for ii = 1:degree
                    coeff(ii) = coeff(ii)/-root;
                    coeff(ii+1) = coeff(ii+1) - coeff(ii);
                end%ii
                coeff(degree + 1) = 0;
            else
                for ii = 1:degree
                    coeff(ii) = coeff(ii+1);
                end%ii
                coeff(degree+1) = 0;
            end%if
        end%function
        function coeff = InPlaceRootAddition(degree,coeff,root)
            %WARNING: This routine is not aware of the memory allocated to
            %"coeff." Make sure that this is being called on arrays with
            %enough memory to store the coefficient resulting from the
            %newly added root.
            
            %WARNING: This routine cannot operate on arrays of with ALL
            %zeros. A zeroth degree polynomial must have a nonzero
            %constant.
            
            %Test:
            %[  6,-5,-2, 1] = (x-3)(x-1)(x+2)
            %[-30,31, 5,-7,1] = (x-3)(x-1)(x+2)(x-5)
            
            %P(x) = Q(x)*(x - root) 
            %Multiply by "x".
            for ii = (degree + 2):-1:2
                coeff(ii) = coeff(ii - 1); 
            end%ii
            coeff(1) = 0;
            %if root ~= 0
                %Multiply by "-root"
                for ii = 1:(degree + 1)
                    coeff(ii) = coeff(ii) - coeff(ii+1)*root;
                end%ii
            %end%if
        end%function
        
        %Polynomial Long Division.
        function [res,rem] = OutOfPlaceDivision(degree1,coeff1,degree2,coeff2)
            %"res" is an output buffer that contains the coefficients of
            %the dividend polynomial.
            %"rem" is an output buffer that contains the coefficients of
            %the remainder of the polynomial.
        
            
            %[res,rem] = polynomial.OutOfPlaceDivision(4,[-4,4,0,5,6],2,[-1,1,2])
            %[res,rem] = polynomial.OutOfPlaceDivision(4,[28,-46,27,-10,1],1,[-7,1])
            
            %Divide polynomial 1 by polynomial 2
            res = zeros(1,degree1 + 1); %Allocate memory for the output buffer.
            
            %Make a copy of the coefficients being divided.
            rem = zeros(1,degree1 + 1);
            for ii = 1:(degree1 + 1)
                rem(ii) = coeff1(ii);
            end%ii
            
            [res,rem] = polynomial.InPlaceDivision(degree1,rem,degree2,coeff2,res);
            %{
            %%%%%%%%%%%%%%%% 
            %This becomes inplace polynomial division.
            for ii = 1:(degree1 - degree2 + 1)
                idx1 = degree1 + 2 - ii;
                idxr = idx1 - degree2;
                res(idxr) = rem(idx1)/coeff2(degree2 + 1);
                for jj = 1:(degree2 + 1)
                    idxr1 = idx1 + 1 - jj;
                    rem(idxr1) = rem(idxr1) - res(idxr)*coeff2(degree2 + 2 - jj);
                end%jj
            end%ii
            %%%%%%%%%%%%%%%
            %}
        end%function
        function [res,coeff1] = InPlaceDivision(degree1,coeff1,degree2,coeff2,res)            
            %Divide the polynomial coefficients encoded inside "coeff1"
            %inplace by the polynomial coefficients encoded inside
            %"coeff2".
            
            %"res" needs to be a sufficiently large buffer to store the
            %result of the division. At the conclusion of this routine,
            %"coeff1" will contain the coefficients of the remainder of the
            %division.
            for ii = 1:(degree1 - degree2 + 1)
                idx1 = degree1 + 2 - ii;
                idxr = idx1 - degree2;
                res(idxr) = coeff1(idx1)/coeff2(degree2 + 1);
                for jj = 1:(degree2 + 1)
                    idxr1 = idx1 + 1 - jj;
                    coeff1(idxr1) = coeff1(idxr1) - res(idxr)*coeff2(degree2 + 2 - jj);
                end%jj
            end%ii
        end%function
        
        %Lagrange Polynomials
        function coeff = LagrangeCoefficients(points,X,Y)
            %points = #Of interpolation points.
            %X = Values of the indeterminate.
            %Y = Values of the determinate.
            coeff = zeros(1,points); %Allocate memory for the final output.
            basis = zeros(1,points); %Allocate memory for an intermediate output.
            for ii = 1:points
                basis = polynomial.LagrangeBasis(points,X,ii,basis);
                for jj = 1:points
                    coeff(jj) = coeff(jj) + basis(jj)*Y(ii);
                end%jj
            end%ii
            clear basis;
        end%function
        function coeff = LagrangeConstituent(points,X,number,coeff)
            %points = #Data points to fit.
            %X = values of the indeterminates.
            %Y = values of the determinates.
            %coeff = buffer for the coefficients of the Lagrange constituent.
            if isempty(coeff)
                coeff = zeros(1,points);
            end%if
            R_buffer = 1;
            roots = zeros(points,1);
            switch number
                case 1
                    for ii = 2:points
                        R_buffer = R_buffer*(X(number) - X(ii));
                        roots(ii - 1) = X(ii);
                    end%ii
                case points
                    for ii = 1:(points - 1)
                        R_buffer = R_buffer*(X(number) - X(ii));
                        roots(ii) = X(ii);
                    end%ii
                otherwise
                    for ii = 1:(number - 1)
                        roots(ii) = X(ii);
                        R_buffer = R_buffer*(X(number) - X(ii));
                    end%ii
                    for ii = (number + 1):points
                        roots(ii-1) = X(ii);
                        R_buffer = R_buffer*(X(number) - X(ii));
                    end%ii
            end%switch
            coeff = polynomial.ExpandRoots(points - 1,roots,coeff);
            for ii = 1:points
                coeff(ii) = coeff(ii)/R_buffer;
            end%ii
            clear roots;            
        end%function
        function coeff = LagrangeRebuildConstituent(points,X,old_number,new_number,coeff)
            %Obtain Lagrange polynomial constituents from an already
            %existing constituent and a number sequence. Computations done
            %inplace.
            if old_number == new_number
                return;
            end%if
            
            for kk = 1:2
                number = old_numer*(kk == 1) + new_number*(kk == 2);
                switch number
                end
            end%kk
            
            
            %First, rescale the coefficients by multiplying by 
            for ii = 1:points
                coeff(ii) = coeff(ii)*(1)
            end%ii
            
        end%function
        

        function weight = LagrangeBarycentricWeight(points,X,number)
            %Generates the weights for the 1st Barycentric form of the
            %Lagrange basis.
            weight = 1;
            switch number
                case 1
                    for ii = 2:points
                        weight = weight*(X(number) - X(ii));
                    end%ii
                case points
                    for ii = 1:(points - 1)
                        weight = weight*(X(number) - X(ii));
                        roots(ii) = X(ii);
                    end%ii
                otherwise
                    for ii = 1:(number - 1)
                        roots(ii) = X(ii);
                        weight = weight*(X(number) - X(ii));
                    end%ii
                    for ii = (number + 1):points
                        roots(ii-1) = X(ii);
                        weight = weight*(X(number) - X(ii));
                    end%ii
            end%switch
            
        end%function
        function coeff = LagrangeFiniteDifference(points,sequence,derivative)
            coeff = zeros(1,points); %Buffer to store finite difference coefficients.
            C = zeros(1,derivative); %Buffer to store integration constants.
            
            %Generate the first barycentric form of the Lagrange basis
            %polynomials.
            barycenter = zeros(1,points+1);
            barycenter(1) = 1;
            degree = 0;
            for ii = 1:points
                barycenter = polynomial.InPlaceRootAddition(degree,barycenter,sequence(ii));
                degree = degree + 1;
            end%ii
            %barycenter = polynomial.ExpandRoots(points,sequence,barycenter);
            for ii = 1:points
                degree = points;
                weight = polynomial.LagrangeBarycentricWeight(points,sequence,ii);
                %Turn the barycenter into an unweighted Lagrange basis
                %constituent.
                barycenter = polynomial.InPlaceRootDeletion(...
                    degree,...
                    barycenter,...
                    sequence(ii));
                degree = degree - 1;
                
                %THIS NEEDS TO BE UPGRADED. INSTEAD OF DIFFERENTIATING AND
                %INTEGRATING BACK, ONE SHOULD JUST USE A MODIFIED VERSION
                %OF HORNER'S RULE.
                %%%%
                %Take as many derivatives as needed.
                for jj = 1:derivative
                    C(jj) = barycenter(1); %Remember the coefficients lost to differentiation.
                    barycenter = polynomial.InPlaceDifferentiation(degree,barycenter);
                    degree = degree - 1; %Differentiation lowers the degree of the polynomial.
                end%jj
                coeff(ii) = barycenter(1)/weight;
                for jj = derivative:-1:1 %Undo the derivatives by integrating
                    barycenter = polynomial.InPlaceIntegration(degree,barycenter,C(jj));
                    degree = degree + 1;
                end%jj
                %%%%
                
                
                barycenter = polynomial.InPlaceRootAddition(degree,barycenter,sequence(ii));
            end%ii
        end%function
        
        %Bernstein Polynomials
        function coeff = BernsteinConstituent
        end%function
        
        %Elementary calculus operations.
        function coeff = OutOfPlaceDifferentiation(degree,coeff_in)
            coeff = zeros(1,degree + 1);
            for ii = 1:(degree + 1)
                coeff(ii) = coeff_in(ii);
            end%ii
            coeff = InPlaceDifferentiation(degree,coeff);
        end%function
        function coeff = InPlaceDifferentiation(degree,coeff)
            for ii =  2:(degree + 1)
                iim1 = ii - 1;
                coeff(iim1) = coeff(ii)*(iim1);
            end%ii
            coeff(ii) = 0; %Must zero the leading coefficient.
        end%function
        
        function coeff = OutOfPlaceIntegration(degree,coeff_in,C)
            coeff = zeros(1,degree + 2);
            for ii = 1:(degree+2)
                coeff(ii) = coeff_in(ii);
            end%ii
            coeff = polynomial.InPlaceIntegration(degree,coeff,C);
        end%function
        function coeff = InPlaceIntegration(degree,coeff,C)
            for ii = (degree + 2):-1:2
                iim1 = ii - 1;
                coeff(ii) = coeff(iim1)/iim1;
            end%ii
            coeff(1) = C; %Include the integration constant.
        end%function
    

        %O(n) horizontal shift reflected on a polynomial's coefficients.
        function coeff_out = TaylorShift(n_coeff,coeff_in,shift)
            coeff_out = zeros(1,n_coeff); %Allocate output buffer.
            
            %First compute g(x) = P(shift*x)
            pb = 1; %Power buffer.
            for ii = 2:n_coeff
                pb = pb*shift;
                coeff_out(ii) = pb*coeff_in(ii);
            end%ii
            
            %Now, compute f(x) = g(x+1) 
            %WARNING: NEED A suitable replacement for the "binexp"
            %function.
            for ii = 2:n_coeff
                
                coeff = binexp(ii-1);
                pb = coeff_out(ii); %Repurpose the buffer to remember the coefficient.
                for jj = 1:(ii - 1)
                    coeff_out(jj) = coeff_out(jj) + pb*coeff(jj);
                end%jj
            end%ii
            
            %Finally compute q(x) = f(x/x0)
            pb = 1; %Reset the power buffer.
            for ii = 2:n_coeff
                pb = pb/shift;
                coeff_out(ii) = pb*coeff_out(ii);
            end%ii
        end%function
        %Horner's method for O(n) polynomial evaluation & variants.
        function Y_out = Horner(n_coeff,coeff,x)
            Y_out = coeff(n_coeff)*x;
            for ii = (n_coeff - 1):-1:2
                Y_out = (Y_out + coeff(ii)).*x;
            end%ii
            Y_out = Y_out + coeff(1);
        end%function
        %Generate points evenly on the polynomial over a specified range.
        function Y_out = EvaluateRange(n_coeff,coeff,N,x1,x2)
            Y_out = zeros(N,1); %Allocate output buffer.
            dx = (x2 - x1)/(N - 1); %Uniform spacing.
            x = x1;
            for ii = 1:N
                Y_out(ii) = polynomial.Horner(n_coeff,coeff,x);
                x = x + dx;
            end%ii
        end%function
        %Generate points on the polynomial using a custom list of values.
        function Y_out = EvaluateCustom(n_coeff,coeff,n_x,x_in)
            Y_out = zeros(n_x,1); %Allocate output buffer.
            for ii = 1:n_x
                Y_out(ii) = Horner(n_coeff,coeff,x_in(ii));
            end%ii
        end%function
        %Add two polynomials
        function coeff_out = Addition(n_coeff_1,coeff_1,n_coeff_2,coeff_2)
            if n_coeff_1 >= n_coeff_2
                coeff_out = zeros(1,n_coeff_1);
                for ii = 1:n_coeff_2
                    coeff_out(ii) = coeff_1(ii) + coeff_2(ii);
                end%ii
                for ii = (ii + 1):n_coeff_1
                    coeff_out(ii) = coeff_1(ii);
                end%ii
            elseif n_coeff_1 < n_coeff_2
                coeff_out = zeros(1,n_coeff_2);
                for ii = 1:n_coeff_1
                    coeff_out(ii) = coeff_1(ii) + coeff_2(ii);
                end%if
                for ii = (ii + 1):n_coeff_2
                    coeff_out(ii) = coeff_2(ii);
                end%ii
            end%if     
        end%function
        %Create line plot of a polynomial evaluated at evenly spaced points
        %in a range
        function graphics = DrawPolynomialRange(ax,n_coeff,coeff,N,x1,x2)
            graphics = custom_line(...
                'Parent',ax,...
                'DisplayName',['Degree-',num2str(n_coeff - 1),'Polynomial'],...
                'X',linspace(x1,x2,N),...
                'Y',polynomial.EvaluateRange(n_coeff,coeff,N,x1,x2));
        end%function
        %Create line plot of a polynomial evaluated at custom values.
        function graphics = DrawPolynomialCustom(n_coeff,coeff,n_x,x_in)
            graphics = custom_line(...
                'Parent',ax,...
                'DisplayName',['Degree-',num2str(n_coeff - 1),'Polynomial'],...
                'X',x_in,...
                'Y',EvaluateCustom(n_coeff,coeff,n_x,x_in));
        end%function
    end%methods (Static)  
    %Elementary low-level functions that are not unique in application to
    %the class.
    methods (Static)
        %Find the greater of two numbers.
        function maxAB = greater(A,B)
            maxAB = A*(A >= B) + B*(B > A);
        end%function
    end%methods (Static)
end%classdef
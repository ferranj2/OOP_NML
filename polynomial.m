classdef polynomial < handle
    properties (SetAccess = public)
        canvas
        sketches
        
        sketch %Line graph of this polynomial %DEPRECATE
    end%properties (public)
    properties (SetAccess = protected)
        degree
        coefficients
    end%properties
    properties (Hidden = true)
        valid %Whether this polynomial is valid.
        
        %Graphics related flags.
        canvas_set
        graphics_initialized
    end%properties
    %High-level functions that MODIFY specific instances of the class
    %object.
    methods
        %Constructor
        function this = polynomial
            
            %Graphics related
            this.canvas_set = false;
            this.graphics_initialized = false;
            this.sketches = struct(...
                'Curve',[],...
                'Roots',[],...
                'Inflection',[],...
                'Extrema',[]);
            
            %Related to the object.
            this.degree = [];
            this.coefficients = [];
            this.sketch = line.empty;
            this.valid = [];
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
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',0,... %Do not initalize with empty ("[]" ) because...
                'YData',0,... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','-',...
                'DisplayName',[num2str(this.degree),'-degree polynomial']);
                            
   
            
            this.graphics_initialized  = true;
        end%function.           
        
    end%methods (Graphics)
    
    
    
    
    %High-level functions that CREATE instances of the polynomial objects
    %from low-level code.Error checking involved.
    methods (Static)
        %Create a polynomial from input coefficients.
        function poly = CreateFromCoefficients(n_coeff,coeff)
            if n_coeff < 1
                error('Must input atleast one coefficient.')
            end%if
            poly.degree = n_coeff - 1;
            poly.coefficients = coeff;
        end%function
        %Create a polynomial by taking the derivative of an already
        %existing polynomial.
        function poly = CreateFromDifferentiation(original)
            poly = polynomial;
            poly.degree = original.degree - 1;
            poly.coefficients = polynomial.Differentiate(...
                original.degree + 1,...
                original.coefficients);
        end%function
        %Create a polynomial by integrating an already existing polynomial.
        %The user should input an integration constant, else, a default of
        %zero is assumed.
        function poly = CreateFromIntegration(original,int_const)
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
        %Create a polynomial by adding two already existing ones.
        function poly = CreateFromAddition(polynomial_1,polynomial_2)
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
        %Create from roots
        function poly = CreateFromRoots(n_roots,roots)
            poly = polynomial;
            poly.degree = n_roots;
        end%function
    end%methods (Static)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Quadratic formula for the roots of degree-2 polynomials.
        function [xr1,xr2] = QuadraticFormula(a,b,c)
            xr1 = -0.5*b/a;
            xr2 = xr1 - sqrt(xr1*xr1 - c/a);
            xr1 = 2*xr1 - xr2;
        end%function
        %Cubic formula for the roots of degree-3 polynomials.
        function [xr1,xr2,xr3] = CubicFormula(a,b,c,d)
        end%function
        %Expand polynomial coefficients from prescribed roots.
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
        
        %Compute the coefficients of a Lagrange Interpolating polynomial.
        %THIS IS BUGGED! NEEDS FIXIGN!!
        function coeff = LagrangeCoefficients(points,X,Y)
            coeff = zeros(1,points); %Allocate memory for the output.
            
            %The Lagrange constituents can be expanded as polynomial roots.
            roots = zeros(1,points-1);
            n_roots = points - 1;
            
            
            %Generate and append the first Lagrange constituent to the
            %running sum of coefficients.
            coeff_roots = zeros(1,points);

            %Construct the last Lagrange constituent and append.
            for ii = 1:n_roots%Initialize the "roots" buffer.
                roots(ii) = X(ii); %All roots except last. (n_roots + 1 = points)
            end%ii
            coeff_roots = polynomial.ExpandRoots(n_roots,roots,coeff_roots);
            R_buffer = 1;
            for ii = 1:n_roots
                R_buffer = R_buffer*(X(points) - X(ii));
            end%ii
            for ii = 1:points %Appending
                coeff(ii) = coeff(ii) + Y(points)*coeff_roots(ii)/R_buffer;
            end%ii

            %Construct the first Lagrange constuent.
            roots(1) = X(points);
            coeff_roots = polynomial.ExpandRoots(n_roots,roots,coeff_roots);
            R_buffer = 1;
            for ii = 2:points
                R_buffer = R_buffer*(X(1) - X(ii));
            end%ii
            for ii = 1:points %Appending
                coeff(ii) = coeff(ii) + Y(1)*coeff_roots(ii)/R_buffer;
            end%ii            
            roots(1) = X(1); %Undo the padding
            
            %Construct the Lagrange constituents in the middle if
            %necessary. Append just like before.
            if points > 2 
                %Assemble the remaining Lagrange constituents.
                for ii = (points - 1):-1:2
                    iip1 = ii + 1;
                    %Skip the point that would result division by zero.
                    %Reset the buffer and reassemble it.
                    R_buffer = 1;
                    for jj = points:-1:iip1
                        R_buffer = R_buffer*(X(ii) - X(jj));
                    end%ii
                    for jj = (ii - 1):-1:2
                        R_buffer = R_buffer*(X(ii) - X(jj));
                    end%jj
                    %Append to the running coefficient buffer.
                    coeff_roots = polynomial.ExpandRoots(n_roots,roots,coeff_roots);
                    for jj = 1:points
                        coeff(jj) = coeff(jj) + Y(ii)*coeff_roots(jj)/R_buffer;
                    end%jj
                    roots(ii) = X(iip1); %Permute the "roots" buffer.
                end%ii
            end%if
            
            
            clear coeff_roots;
        end%function
        %Simple O(n) Differentiation of coefficients
        function coeff_out = Differentiate(n_coeff,coeff_in)
            coeff_out = zeros(1,n_coeff - 1);%Allocate output buffer.
            for ii = 2:n_coeff
                coeff_out(ii-1) = coeff_in(ii)*(ii-1);
            end%ii
        end%function
        %Simple O(n) Integration of coefficients
        function coeff_out = Integrate(n_coeff,coeff_in,int_const)
            coeff_out = zeros(1,n_coeff + 1);%Allocate output buffer.
            coeff_out(1) = int_const;
            for ii = 1:n_coeff
                coeff_out(ii+1) = coeff_in(ii)/ii;
            end%ii
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
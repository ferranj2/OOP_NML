%% ellipse.m
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% A data structure used to generate and manipulates ellipses in the 2D
% plane. For use in geometry and CAD applications. Custom creation routines
% like elliptical fillet, from two arbitrary vectors and center (affine
% circle transform), and from valid quadric curve coefficients. Can
% generate elliptical arcs and perform point inclusion tests.
%% Formulae
% Area
%%
% $A = \pi a b$
%%
% Eccentricity
%%
% $e = \sqrt{1-\left(\frac{b}{a}\right)^{2}}$
%%
% "Linear" eccentricity (distance from center to a focus.
%%
% $c = \sqrt{a^{2} - b^{2}}$
%%
% Quadrics from parameters.
%%
% $\begin{array}{c}
% A = a^{2}\sin^{2}{(\theta)} + b^{2}\cos^{2}{(\theta)}\\
% B = 2(b^{2} - a^{2})\sin{(\theta)}\cos{(\theta)}\\
% C = a^{2}\cos^{2}{(\theta)} + b^[2}\sin^{2}{(\theta)}\\
% D = -2A x_{c} - B y_{c}\\
% E = -B x_{c} - 2C y_{c}\\
% F = A x_{c}^{2} + B x_{c} y_{x} + C y_{c}^2 - a^{2}b^{2}
% \end{array}$
%%
% Measures from parameters.
%%
% Effective polar radius (measured counterclockwise from major radius)
%%
% $r(\theta) = \frac{b}{\sqrt{1 - \left(e\cos{\left(\theta\right)}\right)^{2}}}$
%%
% Affine parametrization
%%
% $E(t) = x_{0} + d_{1} \cos{\left(t\right)} + d_{2} \sin{\left(t\right)}$
%%
% Vertex anomaly
%%
% $\cot{\left(2 t_{0}\right)} = \frac{\|d_{1}\|^{2} - \|d_{2}\|^{2}}{2 d_{1} \cdot d_{2}}$
%% Class definition
classdef ellipse < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the ellipse's graphics as they appear on axes legend.
        color %Color of ellipse's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the graphics.
        sketches %Structured array that contains all the plotting handles.
        
        sketch %DEPRECATE
        sketch_center %DEPRECATE
        sketch_a %DEPRECATE
        sketch_b %DEPRECATE
    end%properties
    %remove
    properties (SetAccess = protected)
        N %Number of points to plot (when calling plotting routines).
    end
    %remove
    %DEFINING DATA variables
    properties (SetAccess = protected)
        center %Center. %DEPRECATE
        a %major radius.
        b %minus radius.
        fa%Unit vector towards major radius from center. %DEPRECATE
        fb%Unit vector towards minor radius from center. %DEPRECATE
        
        xc %X-coordinate of the center.
        yc %Y-coordinate of the center.
        eax %X-cosine of the major radius.
        eay %Y-cosine of the major radius.
        ebx %X-cosine of the minor radius.
        eby %Y-cosine of the minor radius.
    end
    %METRIC variables
    properties (SetAccess = protected)
        A %(Quadric) coefficient of x^2;
        B %(Quadric) coefficient of x*y; 
        C %(Quadric) coefficient of y^2;
        D %(Quadric) coefficient of x;
        E %(Quadric) coefficient of y;
        F %(Quadric) constant coefficient;
        dim %Dimension.
        c %"linear" eccentricity.
        e %eccentricity.
        l %semi-latus rectum.
        area %area.
        d %direction cosines. DEPRECATE?
        n %Unit normal.
        kappa_a %Curvature at the major radius.
        kappa_b %Curvature at the minor radius.
        V1 %First Vertex.
        V2
        V3
        V4
        F1 %First focus
        F2 %Second focus
        generated %Flag to see if the ellipse is generated.
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        canvas_set
        graphics_initialized
    end%properties
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = ellipse(varargin)
            this.graphics_initialized = false;
            this.canvas_set = false;
            this.sketches = struct(...
                'Curve',[],...
                'Orientation',[],...
                'Center',[],...
                'Focci_Labels',[],...
                'Focci',[],...
                'Vertex_Labels',[],...
                'Vertices',[]);
            
            this.N = [];
            this.dim = [];
            this.a = [];
            this.b = [];
            this.area = [];
            this.center = [];
            this.sketch = [];
            this.sketch_center = [];
            this.sketch_a = [];
            this.sketch_b = [];
            this.generated = false;
            %{
            if nargin == 0
                this.generated = false;
                return;
            end
            for ii = 2:nargin
                switch varargin{2*(ii-1)}
                    case 'a'
                        this.a = varargin{2*(ii-1) + 1};
                    case 'b'
                        this.b = varargin{2*(ii-1) + 1};
                    case 'dim'
                        this.dim = varargin{2*(ii-1) + 1};
                end
            end
            %}
        end%function
        
        %Custom creation routines
        function elli = CreateXYABTH(xc,yc,a,b,th)
            elli = ellipse;
            elli.xc = xc;
            elli.yc = yc;
            
            sign_a = (a > 0) - (a < 0);
            sign_b = (b > 0) - (b < 0);
            if sign_a == 0 || sign_b == 0
                error('Either "a" or "b" were input as zero!');
            end%if
            %Set the direction cosines.
            elli.eax = cos(th);
            elli.eay = sin(th);
            elli.ebx = -elli.eay;
            elli.eby = elli.eax;
            
            %Account for the signs of the inputs
            elli.eax = elli.eax*sign_a;
            elli.eay = elli.eay*sign_a;
            elli.ebx = elli.ebx*sign_b;
            elli.eby = elli.eby*sign_b;
            elli.a = a*sign_a;
            elli.b = b*sign_b;
            
        end%function
        function elli = CreateFromQuadric(A,B,C,D,E,F)
            ellip = ellipse;
            
        end%function
        function elli = CreateSteinerEllipse(x0,y0,x1,y1,x2,y2)
            %Create a Steiner ellipse from 2D coordinates of a triangle.
            elli = ellipse;
            
            %Ellipse shares the 
            C = zeros(1,2);
            C(1) = (x0 + x1 + x2)/3;
            C(2) = (y0 + y1 + y2)/3;
            
            %Create the "SC" side.
            f1 = zeros(1,2);
            f1(1) = x2 - C(1);
            f1(2) = y2 - C(2);
            
            %Create the "AB" side.
            f2 = zeros(1,2);
            f2(1) = (x1 - x0)/sqrt(3);
            f2(2) = (y1 - y0)/sqrt(3);
            elli.genfrom2dir1c(f1,f2,C);
        end%function
        function elli = CreateSteinerInEllipse(x0,y0,x1,y1,x2,y2)
            elli = ellipse;
        end%function
        
    end
    %High-level instance MODIFICATION and QUERY routines.
    methods 
        function Measure(this)
            %Basic measures.
            this.area = pi*this.a*this.b;
            this.c = sqrt(this.a*this.a - this.b*this.b);
            this.e = this.c/this.a;
            this.l = this.a*(1 - this.e*this.e);
            
            %Quadric coefficients
            this.A = (this.a*this.eay)^2 + (this.b*this.eax)^2;
            this.B = 2*(this.b*this.b - this.a*this.a)*this.eax*this.eay;
            this.C = (this.a*this.eax)^2 + (this.b*this.eay)^2;
            this.D = -2*this.A*this.xc - this.B*this.yc;
            this.E = -this.B*this.xc - 2*this.C*this.yc;
            this.F = this.A*this.xc^2 + this.B*this.xc*this.yc + this.C*this.yc^2 - (this.a*this.b)^2;
            
        end%function
        
        %Orientation
        function flip_a(this)
            this.fa(1) = this.fa(1)*-1;
            this.fa(2) = this.fa(2)*-1;
            if ~isempty(this.sketch_a)
                this.sketch_a.UData  = this.fa(1)*this.a;
                this.sketch_a.VData  = this.fa(2)*this.a;
            end%if
        end%function
        function flip_b(this)
            this.fb(1) = this.fb(1)*-1;
            this.fb(2) = this.fb(2)*-1;
            if ~isempty(this.sketch_b)
                this.sketch_b.UData  = this.fb(1)*this.b;
                this.sketch_b.VData  = this.fb(2)*this.b;
            end%if
        end%function
        
        %Generate from 2 directions and a center.
        function genfrom2dir1c(this,f1,f2,C)
            dim1 = mustbeflat(f1);
            dim2 = mustbeflat(f2);
            if dim1 ~= dim2
                error('Input directions are of different dimensions!');
            end
            %If no center is specified, assume the center is the origin.
            if nargin < 2
                dimC = dim1;
                this.center = zeros(1,dim1);
            else
                dimC = mustbeflat(C);
                this.center = C;
            end
            if dimC ~= dim1
                error('Dimension of center does not mathc dimension of directions!');
            end
            this.dim = dimC;
            
            %Solve for parameter t
            t0 = acot((inner_product(f1,f1) - inner_product(f2,f2))/(2*inner_product(f1,f2)))/2;
            
            %Allocate space for the vertices.
            this.V1 = zeros(1,this.dim); %Vertex at t = t0
            this.V2 = zeros(1,this.dim); %Vertex at t = t0 + pi/2 
            this.V3 = zeros(1,this.dim); %Vertex at t = t0 + pi/2 
            this.V4 = zeros(1,this.dim); %Vertex at t = t0 + pi/2 

            d1 = zeros(1,this.dim);
            d2 = zeros(1,this.dim);

            for ii = 1:this.dim
                this.V1(ii) = this.center(ii) + f1(ii)*cos(t0) + f2(ii)*sin(t0);
                this.V2(ii) = this.center(ii) + f1(ii)*cos(t0 + pi/2) + f2(ii)*sin(t0 + pi/2);
                this.V3(ii) = this.center(ii) + f1(ii)*cos(t0 + pi) + f2(ii)*sin(t0 +pi);
                this.V4(ii) = this.center(ii) + f1(ii)*cos(t0 - pi/2) + f2(ii)*sin(t0 - pi/2);
                %Compute vectors that point from one oppossing vertex to
                %the other.
                d1(ii) = this.V4(ii) - this.V2(ii);
                d2(ii) = this.V1(ii) - this.V3(ii);
            end
            
            %Make the vectors d1 and d2 unit.
            mag1 = sqrt(inner_product(d1,d1));
            mag2 = sqrt(inner_product(d2,d2));
            for ii = 1:this.dim
                d1(ii) = d1(ii)/mag1;
                d2(ii) = d2(ii)/mag2;
            end

            %Distinguish the major radius from the minor one.
            if mag1 > mag2
                this.a = mag1/2;
                this.b = mag2/2;
                this.fa = d1; %In C, do this with pointers.
                this.fb = d2; %In C, do this with pointers.
            elseif mag2 > mag1
                this.a = mag2/2;
                this.b = mag1/2;
                this.fb = d1; %In C, do this with pointers.
                this.fa = d2; %In C, do this with pointers.
            else %They are equal.
                this.a = mag1/2;
                this.b = this.a;
                this.fb = d2;
                this.fa = d1;
            end
            %Compute other useful ellipse properties.
            this.area = this.a*this.b*pi;
            this.e = sqrt(1-(this.b/this.a)^2);
            this.c = sqrt(this.a^2 - this.b^2);
            
            %Flag that the ellipse has been successfully created.
            this.generated = true;
        end
        
        %Generate from 2 directions and two points.
        function genfrom2dir2p(this,A,dA,B,dB)
            this.dim = mustbeflat(A);
            %Determine where the two points intersect.
            O = intersect2p2d(A,dA,B,dB); %This already error checks.
            
                      
            %Compute vectors from the intersection back to the points.
            OA = zeros(1,this.dim); %Vector from O to A.
            OB = zeros(1,this.dim); %Vector from O to B.
            %f1 = zeros(this.dim,1); %Vector from C to A.
            %f2 = zeros(this.dim,1); %Vector from C to B.
            for ii = 1:this.dim
                OA(ii) = A(ii) - O(ii);
                OB(ii) = B(ii) - O(ii);
                this.center(ii) = O(ii) + OA(ii) + OB(ii); % Default center.
                %Compute vectors from center to intersection points.
                %f1(ii) = A(ii) - this.center(ii);
                %f2(ii) = B(ii) - this.center(ii);
            end
            
            %Find reciprocal basis to p1 and p2.
            mat = [...
            OA(1),OA(2),    0,    0;...
                0,    0,OB(1),OB(2);...
            OB(1),OB(2),    0,    0;...
                0,    0,OA(1),OA(2)];
            RHS = [1,1,0,0]';
            q1q2 = mat\RHS;
            q1 = q1q2([1,2]);
            q2 = q1q2([3,4]);
           

            %Alpha for minimal eccentricity.
            alpha = 2*inner_product(OA,OB)/(inner_product(OA,OA) + inner_product(OB,OB));
            
            %Quadric coefficients.
            AA = (q1(1)^2 + 2*alpha*q1(1)*q2(1) + q2(1)^2);
            BB = 2*(q1(1)*q1(2) + alpha*(q1(1)*q2(2) + q2(1)*q1(2)) + q2(1)*q2(2));
            CC = (q1(2)^2 + 2*alpha*q1(2)*q2(2) + q2(2)^2);
            DD = -2*(q1(1) + q2(1));
            EE = -2*(q1(2) + q2(2));
            FF = 1;
            this.center = O' + (OA + OB)/(1+alpha);

            this.genfromquadric(AA,BB,CC,DD,EE,FF,O' + (OA + OB)/(1+alpha));
            %this.center = O + (OA + OB)/(1+alpha);
            %this.center = O;

            this.genfrom2dir1c(this.a*this.fa,this.b*this.fb,O' + (OA + OB)/(1+alpha));
            %this.e
            this.generated = true;
        end%function
        
        %Generate the ellipse from knowledge of its quadric form.
        function genfromquadric(this,A,B,C,D,E,F,cent)
            discriminant = B*B - 4*A*C;
            %determinant4 = discriminant*F + B*E*D
            this.center(1) = (2*C*D - B*E)/discriminant; %WHY DOES THIS NOT WORK?
            this.center(2) = (2*A*E - B*D)/discriminant; %WHY DOES THIS NOT WORK?
            this.center = cent; %WHY DO I NEED THIS?
            this.a = -sqrt(2*(A*E*E + C*D*D - B*D*E + discriminant*F)*(A + C + sqrt((A-C)^2 + B^2)))/discriminant;
            this.b = -sqrt(2*(A*E*E + C*D*D - B*D*E + discriminant*F)*(A + C - sqrt((A-C)^2 + B^2)))/discriminant;
            this.e = sqrt(1 - (this.b/this.a)^2);
            
            if B ~= 0
                theta = atan((C - A - sqrt((A-C)^2 + B^2))/B);
            end%if
            if B == 0
                if A < C
                    theta = 0;
                end
                if A > C
                    theta = pi/4;
                end
            end
            %{
            this.fa(1) = sin(theta);
            this.fa(2) = cos(theta);
            this.fb(1) = sin(theta + pi/2);
            this.fb(2) = cos(theta + pi/2);
            %}
            
            this.fa(1) = cos(theta);
            this.fa(2) = sin(theta);
            this.fb(1) = cos(theta + pi/2);
            this.fb(2) = sin(theta + pi/2);
            %}
            this.V1 = this.center + this.fa*this.a;
            this.V2 = this.center + this.fb*this.b;
            this.V3 = this.center - this.fa*this.a;
            this.V4 = this.center - this.fb*this.b;

            
            this.generated = true;
        end
        
        %"Imploder" type operation
        function xy_out = attract(this,xy_in)
            %Can only attract point to the ellipse if it is generated.
            if this.generated == false
                error('Ellipse must have been generated prior to this operation!');
            end%if
            %Input validation.
            [row,col] = size(xy_in);
            if this.dim ~= col
                error('Columns of input array do not match dimensionality of the ellipse.')
            end%if
            
            %Attract the points.
            xy_out = zeros(row,col); %Allocate memory for output buffer.
            dir = zeros(1,this.dim); %Allocate a buffer to compute directions to center.
            for ii = 1:row
                %Compute direction from center to exterior point.
                for jj = 1:this.dim
                    dir(jj) =  xy_in(ii,jj) - this.center(jj);
                end%jj
                %Make said direction a unit vector.
                mag = sqrt(inner_product(dir,dir));
                for jj = 1:this.dim
                    dir(jj) = dir(jj)/mag;
                end%jj
                %Find the direction cosine between direction and major radius.:
                cos_th = inner_product(dir,this.fa);
                %Find the distance from the center to the border along dir.
                row = polar_radius_C(this,cos_th);
                for jj = 1:this.dim
                    xy_out(ii,jj) = this.center(jj) + row*dir(jj);
                end%jj
            end%ii
        end%function
        
        %Generate an elliptical sector (angles measured CCW from major
        %radius)
                
        %Generate an elliptical sector subtended by lines connecting two 
        %distinct points each to the center of the ellipse.
        function xy_out = sector2P(this,p1,p2,N)
            %{
            %Define vectors from the points to the center of the ellipse.
            dir1 =  direction2p(this.dim,p1,this.center);
            dir2 =  direction2p(this.dim,p2,this.center);
            
            %Compute the magnitudes of the vectors.
            mag_1 = sqrt(inner_product(dir1,dir1));
            mag_2 = sqrt(inner_product(dir2,dir2));
            
            %Cosines of angles between major and minor radii with the
            %vectors defined above.
            cos_th1a = inner_product(this.fa,dir1)/mag_1; %dir1 vs. fa.
            cos_th2a = inner_product(this.fa,dir2)/mag_2; %dir2 vs. fa.
            %}
            
            %NOTE: IN THE NEAR FUTURE, "p1" and "p2" WILL BE REPLACED WITH
            %"x1", "y1", AND "x2", "y2" RESPECTIVELY.
            %*NEW
            dx1 = p1(1) - this.center(1);
            dy1 = p1(2) - this.center(2);
            dx2 = p2(1) - this.center(1);
            dy2 = p2(2) - this.center(2);
            mag_1 = sqrt(dx1*dx1 + dy1*dy1);
            mag_2 = sqrt(dx2*dx2 + dy2*dy2);
            cos_th1a = (dx1*this.fa(1) + dy1*this.fa(2))/mag_1;
            cos_th2a = (dx2*this.fa(1) + dy2*this.fa(2))/mag_2;
            %*NEW
            
            %Need to determine the hemispheres that the lines belong to.
            %Remember that the directions are defined from the points
            %towards the center.
            %if inner_product(this.fb,dir1) <= 0 
            if (this.fb(1)*dx1 + this.fb(2)*dy1) <= 0
                angle1 = acos(-cos_th1a);
            else 
                angle1 = -acos(-cos_th1a);
            end%if
            %if inner_product(this.fb,dir2) <= 0 
            if (this.fb(1)*dx2 + this.fb(2)*dy2) <= 0
                angle2 = acos(-cos_th2a);
            else 
                angle2 = -acos(-cos_th2a);
            end%if            
            xy_out = sector2A(this,N,angle1,angle2);
        end%function
        
        %Generate discrete points along the ellipse to produce an
        %elliptical sector.
        function xy_out = sector2A(this,N,E1,E2)
            if N <= 1
                error('Must input a positive integer greater than 1.')
            end%if
            xy_out = zeros(N,this.dim); %Allocate output memory.

            %The expectation is that eccentric angle2 is greater than eccentric angle1.
            %Swap operation
            if E1 > E2
                buffer = E1;
                E1 = E2;
                E2 = buffer;
                clear buffer;
            end%if

            d_theta = (E2 - E1)/(N-1);
            theta = E1; %Starting angle.      
            %*NEW
            crossAB = this.fa(1)*this.fb(2) - this.fb(1)*this.fa(2);
            %*NEW
            for ii = 1:N %for all points
                cos_thA = cos(theta);
                cos_thB = cos(pi/2 - theta); %Should I replace this with sine?
                %*OLD
                %{
                %I forget what this is....
                dir = [this.fb(2),-this.fa(2);-this.fb(1),this.fa(1)]*[cos_thA;cos_thB]...
                    /(this.fa(1)*this.fb(2) - this.fb(1)*this.fa(2));
                R = polar_radius_C(this,cos_thA);
                for jj = 1:this.dim
                    xy_out(ii,jj) = this.center(jj) + R*dir(jj);
                end%jj
                %}
                %*OLD
                %*NEW
                R = polar_radius_C(this,cos_thA);
                xy_out(ii,1) = this.center(1) + R*(+this.fb(2)*cos_thA - this.fa(2)*cos_thB)/crossAB;
                xy_out(ii,2) = this.center(2) + R*(-this.fb(1)*cos_thA + this.fa(1)*cos_thB)/crossAB;
                %*NEW
                theta = theta + d_theta; %March the angle.
            end%ii
        end%function
        
        %Determine if a point is inside ellipse.
        function inside = is_inside(this,p1)
            %{
            %Direction from point to center.
            vec =  direction2p(this.dim,p1,this.center);
            
            %Compute the cosine of the angle between the direction just
            %computed and orientation of the major radius.
            cos_th = inner_product(this.fa,vec)/inner_product(vec,vec);
            
            inside = true; %"Innocent until proven innocent."
            if inner_product(vec,vec) > polar_radius_C(this,cos_th)^2
                inside = false;
            end%if
            %}
            dx = p1(1) - this.center(1);
            dy = p1(2) - this.center(2);
            mag2 = dx*dx + dy*dy; 
            cos_th = (dx*this.fa(1) + dy*this.fa(2))/mag2;
            if mag2 > polar_radius_C(this,cos_th)^2
                inside = false;
            end%if            
        end%function
                
        %Polar radius measured from center.
        function R = polar_radius_C(this,cos_th)
           R = this.b/sqrt(1 - (this.e*cos_th)^2);
        end%function
        
        function Draw(this,ax,N)
            if ~this.generated
                error('Cannot draw an invalid ellipse!');
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if nargin < 3
                N = 100;
            end%if
            
            
            if isempty(this.sketch)
                
                
                data = zeros(N,2);
                t = 0;
                dt = 2*pi/(N-1);
                for ii = 1:N
                    cosine = cos(t);
                    sine = sin(t);
                    data(ii,1) = this.center(1) + this.a*this.fa(1)*cosine + this.b*this.fb(1)*sine ;
                    data(ii,2) = this.center(2) + this.a*this.fa(2)*cosine + this.b*this.fb(2)*sine ;
                    t = t+dt;
                end%ii
                this.sketch = custom_line(...
                    'Parent',ax,...
                    'DisplayName','Ellipse',...
                    'Color',[0,0,0]);
                this.sketch.XData = data(:,1);
                this.sketch.YData = data(:,2);
            end%if

        end%function
        
        %Draw the ellipse's major radius.
        function DrawMajorRadius(this,ax)
            if ~this.generated
                error('Cannot draw an invalid ellipse')
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch_a) 
                this.sketch_a = quiver(ax,...
                    this.center(1),... %XData
                    this.center(2),... %YData
                    this.fa(1)*this.a,... %UData
                    this.fa(2)*this.a,... %VData
                    'Color',[1,0,0],...
                    'DisplayName','Major Radius');
            end%if
            
        end%function
        
        %Draw the ellipse's major radius.
        function DrawMinorRadius(this,ax)
            if ~this.generated
                error('Cannot draw an invalid ellipse')
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch_b) 
                this.sketch_b = quiver(ax,...
                    this.center(1),... %XData
                    this.center(2),... %YData
                    this.fb(1)*this.b,... %UData
                    this.fb(2)*this.b,... %VData
                    'Color',[0,0,1],...
                    'DisplayName','Minor Radius');
            end%if
            
        end%function
        
        %Draw the ellipse's center
        function DrawCenter(this,ax)
            if ~this.generated
                error('Cannot draw an invalid ellipse')
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch_center) 
                this.sketch_center = scatter(ax,...
                    this.center(1),... %XData
                    this.center(2),... %YData
                    'Marker','o',... 
                    'MarkerFaceColor',[0,0,1],...
                    'MarkerEdgeColor',[0,0,0]);
            end%if
            
        end%function
        
        %Draw everything about the ellipse.
        function DrawAll(this,ax)
            if ~this.generated
                error('Cannot draw an invalid ellipse')
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            this.Draw(ax);
            this.DrawCenter(ax);
            this.DrawMajorRadius(ax);
            this.DrawMinorRadius(ax);
        end%function
        
        %Draw this ellipse.
        function p1 = plot(this,ax)
            if isempty(this.N) ==1
                this.N = 200;
            end%if
            if this.generated == true
                p1 = line;
                p1.Parent = ax;
                p1.XData = zeros(this.N,1);
                p1.YData = zeros(this.N,1);
                p1.DisplayName = 'Ellipse';
                dt = 2*pi/(this.N - 1);
                t = 0;
                for ii = 1:this.N
                   cosine = cos(t);
                   sine = sin(t);
                   p1.XData(ii) = this.center(1) + this.a*this.fa(1)*cosine + this.b*this.fb(1)*sine ;
                   p1.YData(ii) = this.center(2) + this.a*this.fa(2)*cosine + this.b*this.fb(2)*sine ;
                   t = t+dt;
                end%ii
                
                %Plot center
                s1 = scatter(ax,this.center(1),this.center(2),'k','filled');
                s1.DisplayName = 'Center';
                
                %Plot the major radius.
                qa = quiver(ax,...
                    this.center(1),...
                    this.center(2),...
                    this.a*this.fa(1),...
                    this.a*this.fa(2),...
                    'r');
                qa.DisplayName = 'Major radius';
                %Plot the minor radius.
                qb = quiver(ax,...
                    this.center(1),...
                    this.center(2),...
                    this.b*this.fb(1),...
                    this.b*this.fb(2),...
                    'b');
                qb.DisplayName = 'Minor radius';
                
                %Plot vertices
                %{
                sv = scatter(ax,...
                    [this.V1(1),this.V2(1),this.V3(1),this.V4(1)],...
                    [this.V1(2),this.V2(2),this.V3(2),this.V4(2)],...
                    'm',...
                    'filled');
                sv.DisplayName = 'Vertices';
                
                plot(ax,...
                    [this.V1(1),this.V3(1),this.center(1),this.V2(1),this.V4(1)],...
                    [this.V1(2),this.V3(2),this.center(2),this.V2(2),this.V4(2)],...
                    '--k');
                %}
            else
                warning('Ellipse was not generated. Cannot plot.')
                p1 = [];
            end%if
        end%function

    end%methods (Ordinary)
    %Graphical setups.
    methods
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',[this.XY(:,1);this.XY(1,1)],... %Do not initalize with empty ("[]" ) because...
                'YData',[this.XY(:,2);this.XY(1,2)],... %MATLAB won't allow ANY property access otherwise.
                'Linewidth',1,...
                'LineStyle','-',...
                'DisplayName',this.name);
         
          
            
            this.graphics_initialized  = true;
        end%function.           
        
        
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
    end%methods (Demonstrations)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Compute the quadric coefficients from the ellipse's measures.
        function [A,B,C,D,E,F] = Metrics2Quadrics(xc,yc,a,b,eax,eay)
            A = a*a*eay*eay + b*b*eax*eax;
            B = 2*(b*b - a*a)*eax*eay;
            C = a*a*eax*eax + b*b*eay*eay;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc*xc + B*xc*yc + C*yc*yc - a*a*b*b;
        end%function
        
        %Compute the ellipses' measure from quadric coefficients.
        function [xc,yc,a,b,eax,eay] = Quadrics2Metrics(A,B,C,D,E,F)
            disc = B*B - 4*A*C;
            if disc > 0
                %If discriminant greater than 0, this quadric belongs to a
                %conic other than an ellipse.
                xc = NaN;
                yc = NaN;
                a = NaN;
                b = NaN;
                eax = NaN;
                eay = NaN;
                return;
            end%if
            det4 = A*E*E + C*D*D - B*D*E + disc*F; %This is 4 times the determinant.
            
            xc = (2*C*D - B*E)/disc;
            yc = (2*A*E - B*D)/disc;
            a = -sqrt(2*(det4)*(A + C  + sqrt((A - C)^2 + B*B)))/disc;
            b = -sqrt(2*(det4)*(A + C  - sqrt((A - C)^2 + B*B)))/disc;
            
            th = atan((C- A - sqrt((A - C)^2 + B*B))/B)*(B ~= 0) + pi/4*(B == 0 & A > C);
            eax = cos(th);
            eay = sin(th);
            
        end%function
        
    end%methods
    
    
    
end%classdef
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
        refresh_rate
        
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
        dim %DEPRECATEs
        
        c %"linear" eccentricity.
        e %eccentricity.
        l %semi-latus rectum.
        area %area.
        d %direction cosines. DEPRECATE?
        n %Unit normal.
        kappa_a %Curvature at the major radius.
        kappa_b %Curvature at the minor radius.
        V1 %First Vertex. %DEPRECATE
        V2 %DEPRECATE
        V3 %DEPRECATE
        V4 %DEPRECATE
        F1 %First focus
        F2 %Second focus
        
        %Vertex Coordinates.
        Vx1
        Vy1
        Vx2
        Vy2
        Vx3
        Vy3
        Vx4
        Vy4
        
        %Focci Coordinates
        Fx1
        Fy1
        Fx2
        Fy2
        
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        canvas_set
        graphics_initialized
        
        refresh
        updated
        generated %Flag to see if the ellipse is generated.
        valid %Overall validity of the object.
    end%properties
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = ellipse(varargin)
            this.graphics_initialized = false; %DEPRECATE
            this.canvas_set = false;
            this.name = 'Ellipse';
            this.color = [0,0,0];
            this.sketches = struct(...
                'Curve',[],...
                'MajorRadius',[],...
                'MinorRadius',[],...
                'Center',[],...
                'FocciLabels',[],...
                'Focci',[],...
                'VertexLabels',[],...
                'Vertices',[]);
            this.generated = struct(...
                'Curve',false,...
                'MajorRadius',false,...
                'MinorRadius',false,...
                'Center',false,...
                'FocciLabels',false,...
                'Focci',false,...
                'VertexLabels',false,...
                'Vertices',false);
            this.refresh = struct(...
                'Curve',false,...
                'MajorRadius',false,...
                'MinorRadius',false,...
                'Center',false,...
                'FocciLabels',false,...
                'Focci',false,...
                'VertexLabels',false,...
                'Vertices',false);
            this.updated = struct(...
                'Curve',false,...
                'MajorRadius',false,...
                'MinorRadius',false,...
                'Center',false,...
                'FocciLabels',false,...
                'Focci',false,...
                'VertexLabels',false,...
                'Vertices',false);
            this.refresh_rate = 1;
            
            this.N = []; %DEPRECATE
            this.a = [];
            this.b = [];
            this.e = [];
            this.c = [];
            
            this.Vx1 = [];
            this.Vy1 = [];
            this.Vx2 = [];
            this.Vy2 = [];
            this.Vx3 = [];
            this.Vy3 = [];
            this.Vx4 = [];
            this.Vy4 = [];
            
            this.Fx1  =[];
            this.Fy1  =[];
            this.Fx2  =[];
            this.Fy2  =[];
            
            this.area = [];
            
            
            this.center = []; %DEPRECATE
            this.sketch = []; %DEPRECATE
            this.sketch_center = []; %DEPRECATE
            this.sketch_a = []; %DEPRECATE
            this.sketch_b = []; %DEPRECATE
            this.valid = false;

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
        function elli = CreateFromFocciAndMajor(Fx1,Fy1,Fx2,Fy2,a)
            %The orientation of the Major Radius is from focus 1 to focus 2.
            elli = ellipse;
            elli.a = a;
            elli.Fx1 = Fx1;
            elli.Fy1 = Fy1;
            elli.Fx2 = Fx2;
            elli.Fy2 = Fy2;
            [elli.xc,elli.yc,elli.eax,elli.eay,elli.ebx,elli.eby,elli.b,elli.c] = ellipse.FocciAndMajor2Metrics(Fx1,Fy1,Fx2,Fy2,a);
            elli.MeasureBasic;
            elli.MeasureVertices;
        end%function
        function elli = CreateFromFocciAndMinor(Fx1,Fy1,Fx2,Fy2,b)
            elli = ellipse;
            
        end%function
        
        function elli = CreateFromQuadric(A,B,C,D,E,F)
            elli = ellipse;
            
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
        
        
        function elli = CreateLeastEccentricFillet(x1,y1,dx1,dy1,x2,y2,dx2,dy2)
            elli = ellipse;
            %Determine where the two lines intersect.
            [Ox,Oy,~,~] = ellipse.Intersect2p2d(...
                x1,y1,...
                dx1,dy1,...
                x2,y2,...
                dx2,dy2);
                  
            %Compute vectors from the intersection back to the points.
            OAx = x1 - Ox;
            OAy = y1 - Oy;
            OBx = x2 - Ox;
            OBy = y2 - Oy;

            %Find reciprocal basis to p1 and p2.
            %{
            

    [  A32/(A11*A32 - A12*A31),                        0, -A12/(A11*A32 - A12*A31),                        0]
    [ -A31/(A11*A32 - A12*A31),                        0,  A11/(A11*A32 - A12*A31),                        0]
    [                        0,  A44/(A23*A44 - A24*A43),                        0, -A24/(A23*A44 - A24*A43)]
    [                        0, -A43/(A23*A44 - A24*A43),                        0,  A23/(A23*A44 - A24*A43)]          
            
    +A32/(A11*A32 - A12*A31)
    -A31/(A11*A32 - A12*A31)
    +A44/(A23*A44 - A24*A43)
    -A43/(A23*A44 - A24*A43)
            %}

          
            q1x = +OBy/(OAx*OBy - OAy*OBx);
            q1y = -OBx/(OAx*OBy - OAy*OBx);
            q2x = +OAy/(OBx*OAy - OBy*OAx);
            q2y = -OAx/(OBx*OAy - OBy*OAx);


            %Alpha for minimal eccentricity.
            alpha = 2*(OAx*OBx + OAy*OBy)/(OAx*OAx + OAy*OAy + OBx*OBx + OBy*OBy);

            %Quadric coefficients.
            AA = (q1x^2 + 2*alpha*q1x*q2x + q2x^2);
            BB = 2*(q1x*q1y + alpha*(q1x*q2y + q2x*q1y) + q2x*q2y);
            CC = (q1y^2 + 2*alpha*q1y*q2y + q2y^2);
            DD = -2*(q1x + q2x);
            EE = -2*(q1y + q2y);
            FF = 1;
        
            elli.xc = Ox + (OAx + OBx)/(1 + alpha);
            elli.yc = Oy + (OAy + OBy)/(1 + alpha);

            %elli.genfromquadric(AA,BB,CC,DD,EE,FF,O' + (OA + OB)/(1+alpha));
            elli.genfromquadric(AA,BB,CC,DD,EE,FF,[elli.xc;elli.yc]);

            %this.center = O + (OA + OB)/(1+alpha);
            %this.center = O;

            elli.genfrom2dir1c(elli.a*elli.fa,elli.b*elli.fb,[elli.xc;elli.yc]);
            
            elli.eax = elli.fa(1);
            elli.eay = elli.fa(2);
            elli.ebx = elli.fb(1);
            elli.eby = elli.fb(2);
            
            
            %NEW
            %{
            [elli.xc,elli.yc,elli.a,elli.b,elli.eax,elli.eay,elli.ebx,elli.eby] = ellipse.Quadrics2Metrics(AA,BB,CC,DD,EE,FF);
            elli.fa(1) = elli.eax;
            elli.fa(2) = elli.eay;
            elli.fb(1) = elli.ebx;
            elli.fb(2) = elli.eby;
            %}
        end%function
        
        function elli = CreateLEEF(x1,y1,x2,y2,x3,y3,x4,y4)
            %Fit an ellipse that goes through two points "x1,y1" and
            %"x2,y2" such that it is also tangent to directions
            %"x2-x1,y2-y1" and "x4-x3,y4-y3."
            elli = ellipse;
            %Determine where the two lines intersect.
            [Ox,Oy,~,~] = ellipse.Intersect2p2p(...
                x1,y1,...
                x2,y2,...
                x3,y3,...
                x4,y4);          
                  
            %Compute vectors from the intersection back to the points.
            OAx = x1 - Ox;
            OAy = y1 - Oy;
            OBx = x3 - Ox;
            OBy = y3 - Oy;
            
            %Find reciprocal basis to p1 and p2.
            q1x = +OBy/(OAx*OBy - OAy*OBx);
            q1y = -OBx/(OAx*OBy - OAy*OBx);
            q2x = +OAy/(OBx*OAy - OBy*OAx);
            q2y = -OAx/(OBx*OAy - OBy*OAx);

            %Alpha for minimal eccentricity.
            alpha = 2*(OAx*OBx + OAy*OBy)/(OAx*OAx + OAy*OAy + OBx*OBx + OBy*OBy);

            %{
            syms q1x q1y q2x q2y x y a
            eq = (q1x*x + q1y*y - 1)^2 + (q2x*x + q2y*y-1)^2 + 2*a*(q1x*x + q1y*y)*(q2x*x+q2y*y) - 1
            q1x^2*x^2 + 2*a*q1x*q2x*x^2 + 2*q1x*q1y*x*y + 2*a*q1x*q2y*x*y - 2*q1x*x + q2x^2*x^2 + 2*a*q2x*q1y*x*y + 2*q2x*q2y*x*y - 2*q2x*x + q1y^2*y^2 + 2*a*q1y*q2y*y^2 - 2*q1y*y + q2y^2*y^2 - 2*q2y*y + 1
            
            (q1x^2+ 2*a*q1x*q2x + q2x^2)*x^2
            (2*q1x*q1y + 2*a*q1x*q2y+ 2*a*q2x*q1y + 2*q2x*q2y)*x*y
            (- 2*q1x - 2*q2x)*x
            (q1y^2 + 2*a*q1y*q2y + q2y^2)*y^2 - 2*q1y*y - 2*q2y*y + 1
            %}
            
            
            
            %Quadric coefficients.
            A = (q1x^2 + 2*alpha*q1x*q2x + q2x^2);
            B = 2*(q1x*q1y + alpha*(q1x*q2y + q2x*q1y) + q2x*q2y);
            C = (q1y^2 + 2*alpha*q1y*q2y + q2y^2);
            D = -2*(q1x + q2x);
            E = -2*(q1y + q2y);            
            F = 1;
        
            [elli.xc,elli.yc,elli.a,elli.b,elli.eax,elli.eay,elli.ebx,elli.eby] = ellipse.Quadrics2Metrics(A,B,C,D,E,F);

            elli.xc = Ox + (OAx + OBx)/(1 + alpha);
            elli.yc = Oy + (OAy + OBy)/(1 + alpha);            
            elli.MeasureAll;
            %elli.genfrom2dir1c(elli.a*elli.fa,elli.b*elli.fb,[elli.xc;elli.yc]);
        end%function
        
        function elli = CreateLEEF2(x1,y1,x2,y2,x3,y3,x4,y4)
            %Fit an ellipse that goes through two points "x1,y1" and
            %"x2,y2" such that it is also tangent to directions
            %"x2-x1,y2-y1" and "x4-x3,y4-y3."
            elli = ellipse;
            %Determine where the two lines intersect.
            [Ox,Oy,~,~] = ellipse.Intersect2p2p(...
                x1,y1,...
                x2,y2,...
                x3,y3,...
                x4,y4);          
                  
            %Compute vectors from the intersection back to the points.
            OAx = x1 - Ox;
            OAy = y1 - Oy;
            OBx = x3 - Ox;
            OBy = y3 - Oy;

            %Alpha for minimal eccentricity.
            alpha = 2*(OAx*OBx + OAy*OBy)/(OAx*OAx + OAy*OAy + OBx*OBx + OBy*OBy);
            elli.xc = Ox + (OAx + OBx)/(1 + alpha);
            elli.yc = Oy + (OAy + OBy)/(1 + alpha);
            %fprintf('%f\n',elli.xc)
            %fprintf('%f\n',elli.yc)

            [elli.a,elli.eax,elli.eay,elli.b,elli.ebx,elli.eby,elli.Vx1,elli.Vy1,elli.Vx2,elli.Vy2] = ellipse.RadiiFromDirections(...
                elli.xc,...
                elli.yc,...
                x1 - elli.xc, y1 - elli.yc,...
                x3 - elli.xc, y3 - elli.yc);
            %x1 - elli.xc, y1 - elli.yc,...
            %x3 - elli.xc, y3 - elli.yc);
            elli.MeasureAll;
            %elli.genfrom2dir1c(elli.a*elli.fa,elli.b*elli.fb,[elli.xc;elli.yc]);
        end%function
        
    end%methods (Static)
    %High-level instance MODIFICATION and QUERY routines.
    methods 
        
        function Authenticate(this)
            %This is meant to be used during the creation and redefinition
            %routines to make sure that the inputs are workable. 
            this.valid = false; %"Guilty until proven innocent."
            
            %Major radius inspection.
            if isempty(this.a)
                warning('Ellipse does not have a major radius set!')
                return;
            end%if
            if this.a <= 0
                warning('Ellipse does not have a positive major radius!')
                return;
            end%if
            
            %Minor radius inspection.
            if isempty(this.b)
                warning('Ellipse does not have a minor radius set!')
                return;
            end%if
            if this.b <= 0
                warning('Ellipse does not have a positive minor radius set!')
                return;
            end%if
            
            %Center coordinates inspection
            if isempty(this.xc)
                warning('Ellipse does not have an x-coordinate set!')
                return;
            end%if
            if isempty(this.yc)
                warning('Ellipse does not have an y-coordinate set!')
                return;
            end%if
            
            %Orientation inspection
            if isempty(this.eax)
                warning('Ellipse does not have an x-cosine for the major radius!')
                return;
            end%if
            if isempty(this.eay)
                warning('Ellipse does not have an y-cosine for the major radius!')
                return;
            end%if
            if isempty(this.ebx)
                warning('Ellipse does not have an x-cosine for the minor radius!')
                return;
            end%if
            if isempty(this.eby)
                warning('Ellipse does not have an y-cosine for the minor radius!')
                return;
            end%if
            if abs(this.ebx*this.eax + this.eby*this.eay) > 1e-8
                warning('Orthogonality has appreciable error!')
                %return; 
            end%if
            %If it makes it through, its valid!
            this.valid = true;
        end%function
        
        %Metric property queries
        function MeasureAll(this)
            MeasureBasic(this)
            MeasureVertices(this);
            MeasureFocci(this);
        end%function
        function MeasureBasic(this)
            this.area = pi*this.a*this.b; 
            this.e = this.c/this.a; %eccentricity
            this.l = this.a*(1 - this.e*this.e); %Semi latus rectum
        end%function
        function MeasureVertices(this)
            %Vertex in direction of major radius.
            this.Vx1 = this.xc + this.a*this.eax;
            this.Vy1 = this.yc + this.a*this.eay;

            %Vertex in direction of minor radius.
            this.Vx2 = this.xc + this.b*this.ebx;
            this.Vy2 = this.yc + this.b*this.eby;
            
            %Vertex opposite to major radius.
            this.Vx3 = this.xc - this.a*this.eax;
            this.Vy3 = this.yc - this.a*this.eay;
            
            %Vertex in direction of minor radius.
            this.Vx4 = this.xc - this.b*this.ebx;
            this.Vy4 = this.yc - this.b*this.eby;
        end%function
        function MeasureFocci(this)
            this.c = sqrt(this.a*this.a - this.b*this.b); %"linear" eccentricity
            this.Fx1 = this.xc + this.c*this.eax;
            this.Fy1 = this.yc + this.c*this.eay;
            this.Fx2 = this.xc - this.c*this.eax;
            this.Fy2 = this.yc - this.c*this.eay;
        end%function
        function MeasureBounds(this)
            %
            %WIP for measuring the AABB of the ellips.
            %
        end%function
        
        %Redefinition routines
        function RedefineFromQuadrics(this,A,B,C,D,E,F)
            %Morph the ellipse to a new elliptical shape according to
            %Quadric coefficients.
            discriminant = B*B - 4*A*C;
            if discriminant >= 0
                error('Quadric discriminant greater than 0! (Quadric is not an ellipse)');
            end%if
            Delta4 = A*E*E + C*D*D - B*D*E + discriminant*F;
            this.a = -sqrt(2*Delta4*(A + C + sqrt((A-C)^2 + B^2)))/discriminant;
            this.b = -sqrt(2*Delta4*(A + C - sqrt((A-C)^2 + B^2)))/discriminant;
            if B ~= 0
                theta = atan((C - A - sqrt((A-C)^2 + B^2))/B);
            else
                if A < C
                    theta = 0;
                end%if
                if A > C
                    theta = pi/4;
                end%if
            end%if
            this.eax = cos(theta);
            this.eay = sin(theta);
            this.ebx = cos(theta + pi/2);
            this.eay = sin(theta + pi/2);
            
            this.MeasureAll; %Get new metric data.
            this.Refresh; %Refresh any graphics.
        end%function
        function RedefineFromFocciAndMajor(this,Fx1,Fy1,Fx2,Fy2,a)
            [this.xc,this.yc,this.eax,this.eay,this.ebx,this.eby,this.b,this.c] = FocciAndMajor2Metrics(Fx1,Fy1,Fx2,Fy2,a);
            this.a = a;
            this.Fx1 = Fx1;
            this.Fy1 = Fy1;
            this.Fx2 = Fx2;
            this.Fy2 = Fy2;
            this.MeasureBasic;
            this.MeasureVertices;
            this.Refresh;
        end%function
        
        %Orientation
        function ReverseMinorRadius(this)
            this.ebx = -this.ebx;
            this.eby = -this.eby;

            RefreshMinorRadius(this);
        end%function
        function ReverseMajorRadius(this)
            this.eax = -this.eax;
            this.eay = -this.eay;
            
            RefreshMajorRadius(this);
        end%function
        
        %Query routines
        function XY = CreatePolarSector(this,th1,th2,N)
            %Creates an elliptical sector using the center-based polar
            %form. The sector is subtended between a minimum angle "th1",
            %and maximum angle "th2", both measured from the oriented major
            %radius. The sector will consist of "N" evenly spaced (by polar
            %angle) points.
            if nargin < 4
                N = 30;
            end%if
            if N <= 1
                error('Must input a positive integer greater than 1.')
            end%if
            XY = zeros(N,2); %Allocate output memory.
            %Make sure the inputs are in ascending order.
            if th1 > th2
                [th1,th2] = ellipse.swap(th1,th2);
            end%if

            d_theta = (th2 - th1)/(N-1);
            theta = th1; %Starting angle.      
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
                XY(ii,1) = this.center(1) + R*(+this.fb(2)*cos_thA - this.fa(2)*cos_thB)/crossAB;
                XY(ii,2) = this.center(2) + R*(-this.fb(1)*cos_thA + this.fa(1)*cos_thB)/crossAB;
                theta = theta + d_theta; %March the angle.
            end%ii
        end%function
        function XY = CreateEccentricSector(this,N,E1,E2) 
        end%function
        function XY = CreateFocalSector(this,N,E1,E2) 
        end%function
        
        
        %DEPRECATE
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
        %DEPRECATE


        
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
            this.dim = 2;
            
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
            if this.valid == false
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
            dx1 = this.center(1) - p1(1);
            dy1 = this.center(2) - p1(2);
            dx2 = this.center(1) - p2(1);
            dy2 = this.center(2) - p2(2);
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
        %DEPRECATE
        function R = polar_radius_C(this,cos_th)
           R = this.b/sqrt(1 - (this.e*cos_th)^2);
        end%function
        %DEPRECATE
        %DEPRECATE
        function Draw(this,ax,N)
            if ~this.valid
          %     error('Cannot draw an invalid ellipse!');
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
        %DEPRECATE
        %DEPRECATE
        %Draw the ellipse's major radius.
        function DrawMajorRadius(this,ax)
            if ~this.valid
            %    error('Cannot draw an invalid ellipse')
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
        %DEPRECATE
        %DEPRECATE
        %Draw the ellipse's major radius.
        function DrawMinorRadius(this,ax)
            if ~this.valid
             %   error('Cannot draw an invalid ellipse')
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
        %DEPRECATE
        %DEPRECATE
        %Draw the ellipse's center
        function DrawCenter(this,ax)
            if ~this.valid
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
        %DEPRECATE
        %DEPRECATE
        %Draw everything about the ellipse.
        function DrawAll(this,ax)
            if ~this.valid
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
        %DEPRECATE
        
        %DEPRECATE
        %Draw this ellipse.
        function p1 = plot(this,ax)
            if isempty(this.N) ==1
                this.N = 200;
            end%if
            if this.valid == true
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
        %DEPRECATE
        
    end%methods (Ordinary)
    %Graphical setups.
    methods
        %Graphical setting functions
        function SetRefreshRate(this,rate)
            %Set a target refresh rate for the Refresh routine.
            if rate == 0
                rate = 1;
            end%if
            this.refresh_rate = rate*(-1*(rate < 0));
        end%function
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.generated.Curve
                this.sketches.Curve.Parent = ax;
            end%if
            if this.generated.Center
                this.sketches.Center.Parent = ax;
            end%if
            if this.generated.MajorRadius
                this.sketches.MajorRadius.Parent = ax;
            end%if
            if this.generated.MinorRadius
                this.sketches.MinorRadius.Parent = ax;
            end%if
            if this.generated.Focci
                this.sketches.Focci.Parent = ax;
            end%if
            if this.generated.Vertices
                this.sketches.Vertices.Parent = ax;
            end%if
        end%function
        function SetName(this,string)
            this.name = string;
        end%function
        function SetColor(this,RGB)
            this.color = RGB;
        end%function
        function SetCurvePoints(this,N)
            %The default is 100 points. Finer/Coarser ellipses can be
            %rendered by inputting higher/lower values for N.
            if this.generated.Curve
                delta = N - this.sketches.Curve.UserData;
                if delta > 0
                    %Allocate more space.
                    this.sketches.Curve.XData = [this.sketches.Curve.XData,zeros(1,delta)];
                    this.sketches.Curve.YData = [this.sketches.Curve.YData,zeros(1,delta)];
                elseif delta < 0
                    this.sketches.Curve.XData(end-delta:end) = [];
                    this.sketches.Curve.YData(end-delta:end) = [];
                end%if
                this.sketches.Curve.UserData = N;
                RefreshCurve(this);
            else
                warning('Curve has not been created.')
            end%if
        end%function
      
        
        %Create the Graphical objects.
        %DEPRECATE
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',zeros(1,100),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(1,100),... %MATLAB won't allow ANY property access otherwise.
                'Linewidth',1,...
                'LineStyle','-',...
                'DisplayName',this.name);
            
            this.sketches.Major_Radius = quiver(...
                this.canvas,...
                this.xc,...
                this.yc,... %Do not initalize with empty ("[]" ) because...
                this.eax,... %MATLAB won't allow ANY property access otherwise.
                this.eay,...
                'DisplayName',this.name);
            
            this.sketches.Minor_Radius = quiver(...
                this.canvas,...
                this.xc,...
                this.yc,... %Do not initalize with empty ("[]" ) because...
                this.ebx,... %MATLAB won't allow ANY property access otherwise.
                this.eby,...
                'DisplayName',this.name);
            
            this.graphics_initialized  = true;
            
            %{
            this.sketches = struct(...
                'Curve',[],...
                'Major_Radius',[],...
                'Minor_Radius',[],...
                'Center',[],...
                'FocciLabels',[],...
                'Focci',[],...
                'VertexLabels',[],...
                'Vertices',[]);
            %}
        end%function.           
        %DEPRECATE
        
        %Visibility toggling functions
        function Toggle(this,varargin)
            %Sets the MATLAB "Visible" flag on the various graphics to "on"
            %or "off" depending on the flag's value at the time of calling
            %it. If a graphics object is "on", this routine will set it to
            %"off", the converse is true. If the instance of the object has
            %been subject to transformations in the background and the
            %corresponding graphics are not up to date, this routine will
            %call the appropriate refresh routine before toggling
            %visibility to guarantee that what is shown is up to date.
            
            %NOTE: It is possible to avoid the if-else "ladder" by using
            %MATLAB's notation struct.(*string*). C-programming does not
            %allow this, so here I settle with the "ladder" until I find
            %something more reminiscent of what is allowed in C.
            %NOTE: GNU has created something called "gperf" which may help
            %with this issue when porting to C and OpenGL.            
            
            for ii = 1:(nargin - 1)
                if strcmpi(varargin{ii},'Curve') == 1
                    if ~this.generated.Curve
                        GenerateCurve(this);
                    end%if
                    if ~this.updated.Curve 
                        RefreshCurve(this);
                    end%if
                    VisibilityToggle(this,this.sketches.Curve);
                elseif strcmpi(varargin{ii},'Vertices') == 1
                    if ~this.generated.Vertices
                        GenerateVertices(this);
                    end%if
                    if ~this.updated.Vertices
                        RefreshVertices(this);
                    end%if
                    VisibilityToggle(this,this.sketches.Vertices);
                elseif strcmpi(varargin{ii},'Center') == 1
                    if ~this.generated.Center
                        GenerateCenter(this);
                    end%if
                    if ~this.updated.Center
                        RefreshCenter(this);
                    end%if
                    VisibilityToggle(this,this.sketches.Center);
                elseif strcmpi(varargin{ii},'Focci') == 1
                    if ~this.generated.Focci
                        GenerateFocci(this);
                    end%if
                    if ~this.updated.Focci
                        RefreshFocci(this);
                    end%if
                    VisibilityToggle(this,this.sketches.Focci);
                elseif strcmpi(varargin{ii},'MajorRadius') == 1
                    if ~this.generated.MajorRadius
                        GenerateMajorRadius(this);
                    end%if
                    if ~this.updated.MajorRadius
                        RefreshMajorRadius(this);
                    end%if
                    VisibilityToggle(this,this.sketches.MajorRadius);
                elseif strcmpi(varargin{ii},'MinorRadius') == 1
                    if ~this.generated.MinorRadius
                        GenerateMinorRadius(this);
                    end%if
                    if ~this.updated.MinorRadius
                        RefreshMinorRadius(this);
                    end%if
                    VisibilityToggle(this,this.sketches.MinorRadius);
                elseif strcmpi(varargin{ii},'VertexLabels') == 1
                    if ~this.generated.VertexLabels
                        GenerateVertexLabels(this);
                    end%if
                    if ~this.updated.VertexLabels
                        RefreshVertexLabels(this);
                    end%if
                    VisibilityToggle(this,this.sketches.VertexLabels);
                elseif strcmpi(varargin{ii},'FocciLabels') == 1
                    if ~this.generated.FocciLabels
                        GenerateFocciLabels(this);
                    end%if
                    if ~this.updated.FocciLabels
                        RefreshFocciLabels(this);
                    end%if
                    VisibilityToggle(this,this.sketches.FocciLabels);
                else
                    warning('Unrecognizable graphics option. Valid options are: "Curve","Normals","Centroid", "VertexLabels", and "Segment_Labels"');
                end%if
            end%ii
        end%function
        function VisibilityToggle(this,handle)
            if strcmp(handle.Visible,'on') == 1
                handle.Visible = 'off';
            else
                handle.Visible = 'on';
            end%if
        end%function
        
        %Graphics Refresh Routines
        function Refresh(this)
            tic
            %Set all state of date flags to false. These are all set back
            %to true by the refresh routines if the refresh flag for each
            %curve is set to true.
            this.updated.Curve = false;
            this.updated.Vertices = false;
            this.updated.Center = false;
            this.updated.Focci = false;
            this.updated.MajorRadius = false;
            this.updated.MinorRadius = false;
            this.updated.VertexLabels = false;
            this.updated.FocciLabels = false;

            if this.refresh.Curve 
                RefreshCurve(this);
            end%if
            if this.refresh.Vertices 
                RefreshVertices(this);
            end%if
            if this.refresh.Center 
                RefreshCenter(this);
            end%if
            if this.refresh.Focci 
                RefreshFocci(this);
            end%if
            if this.refresh.MajorRadius 
                RefreshMajorRadius(this);
            end%if
            if this.refresh.MinorRadius 
                RefreshMinorRadius(this);
            end%if
            if this.refresh.VertexLabels 
                RefreshVertexLabels(this);
            end%if
            if this.refresh.FocciLabels 
                RefreshFocciLabels(this);
            end%if
            
            time = 1/this.refresh_rate - toc;
            pause(0 + time*(time > 0));
            drawnow
        end%function.
        function RefreshCurve(this)
            fax = this.a*this.eax;
            fay = this.a*this.eay;
            fbx = this.b*this.ebx;
            fby = this.b*this.eby;
            dt = 2*pi/(this.sketches.Curve.UserData - 1);
            t = 0;
            for ii = 1:(this.sketches.Curve.UserData)
                cosine = cos(t);
                sine = sin(t);
                this.sketches.Curve.XData(ii) = this.xc + fax*cosine + fbx*sine ;
                this.sketches.Curve.YData(ii) = this.yc + fay*cosine + fby*sine ;
                t = t+dt;
            end%ii
            this.updated.Curve = true;
        end%function
        function RefreshVertices(this)
            this.sketches.Vertices.XData(1) = this.Vx1;
            this.sketches.Vertices.YData(1) = this.Vy1;

            this.sketches.Vertices.XData(2) = this.Vx2;
            this.sketches.Vertices.YData(2) = this.Vy2;

            this.sketches.Vertices.XData(3) = this.Vx3;
            this.sketches.Vertices.YData(3) = this.Vy3;
            
            this.sketches.Vertices.XData(4) = this.Vx4;
            this.sketches.Vertices.YData(4) = this.Vy4;
            
            this.updated.Vertices = true;
        end%function
        function RefreshCenter(this)
            this.sketches.Center.XData = this.xc;
            this.sketches.Center.YData = this.yc;
            this.updated.Center = false;
        end%function
        function RefreshFocci(this)
            this.sketches.Focci.XData(1) = this.Fx1;
            this.sketches.Focci.YData(1) = this.Fy1;
            this.sketches.Focci.XData(2) = this.Fx2;
            this.sketches.Focci.YData(2) = this.Fy2;
            this.updated.Focci = true;
        end%function
        function RefreshMajorRadius(this)
            this.sketches.MajorRadius.XData = this.xc;
            this.sketches.MajorRadius.YData = this.yc;
            this.sketches.MajorRadius.UData = this.eax*this.a;
            this.sketches.MajorRadius.VData = this.eay*this.a;

            this.updated.MajorRadius = true;
        end%function
        function RefreshMinorRadius(this)
            this.sketches.MinorRadius.XData = this.xc;
            this.sketches.MinorRadius.YData = this.yc;
            this.sketches.MinorRadius.UData = this.ebx*this.b;
            this.sketches.MinorRadius.VData = this.eby*this.b;
            this.updated.MinorRadius = true;
        end%function
        function RefreshVertexLabels(this)
            %Vertex in direction of major radius.
            this.sketches.VertexLabels(1).String = {'V$_{1}$';num2str(this.Vx1),num2str(this.Vy1)};
            this.sketches.VertexLabels(1).Position(1) = this.Vx1;
            this.sketches.VertexLabels(1).Position(2) = this.Vy1;
            
            %Vertex in direction of minor radius.
            this.sketches.VertexLabels(2).String = {'V$_{2}$';num2str(this.Vx2),num2str(this.Vy2)};
            this.sketches.VertexLabels(2).Position(1) = this.Vx2;
            this.sketches.VertexLabels(2).Position(2) = this.Vy2;
            
            %Vertex in direction opposite to major radius.
            this.sketches.VertexLabels(3).String = {'V$_{3}$';num2str(this.Vx3),num2str(this.Vy3)};
            this.sketches.VertexLabels(3).Position(1) = this.Vx3;
            this.sketches.VertexLabels(3).Position(2) = this.Vy3;
            
            %Vertex in direction opposite to minor radius.
            this.sketches.VertexLabels(4).String = {'V$_{4}$';num2str(this.Vx4),num2str(this.Vy4)};
            this.sketches.VertexLabels(4).Position(1) = this.Vx4;
            this.sketches.VertexLabels(4).Position(2) = this.Vy4;
            
            this.updated.VertexLabels = true;
        end%function
        function RefreshFocciLabels(this)
            %Vertex in direction of major radius.
            this.sketches.FocciLabels(1).String = {'F$_{1}$';num2str(this.Fx1),num2str(this.Fy1)};
            this.sketches.FocciLabels(1).Position(1) = this.Fx1;
            this.sketches.FocciLabels(1).Position(2) = this.Fy1;
            %Vertex in direction of major radius.
            this.sketches.FocciLabels(2).String = {'F$_{2}$';num2str(this.Fx2),num2str(this.Fy2)};
            this.sketches.FocciLabels(2).Position(1) = this.Fx2;
            this.sketches.FocciLabels(2).Position(2) = this.Fy2;
            
            this.updated.FocciLabels = true;
        end%function.
        
        %Graphics generation functions.
        function GenerateCurve(this)
            GenerateDefaultCanvas(this);
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',zeros(1,100),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(1,100),... %MATLAB won't allow ANY property access otherwise.
                'Linewidth',1,...
                'UserData',100,...
                'LineStyle','-',...
                'Visible','off',...
                'DisplayName',this.name);
            this.generated.Curve = true;
        end%function
        function GenerateCenter(this)
            GenerateDefaultCanvas(this);
            this.sketches.Center = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',0,... %Do not initalize with empty ("[]" ) because...
                'YData',0,... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','none',...
                'MarkerFaceColor',this.color,...
                'MarkerEdgeColor',[0,0,0],...
                'Marker','o',...
                'DisplayName',[this.name, 'Center'],...
                'Visible','off');
            this.generated.Center = true;
        end%function
        function GenerateVertices(this)
            GenerateDefaultCanvas(this);
            this.sketches.Vertices = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',[0,0,0,0],... %Do not initalize with empty ("[]" ) because...
                'YData',[0,0,0,0],... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','none',...
                'MarkerFaceColor',this.color,...
                'MarkerEdgeColor',[0,0,0],...
                'Marker','s',...
                'DisplayName',[this.name, 'Vertices'],...
                'Visible','off');
            this.generated.Curve = true;
        end%function
        function GenerateFocci(this)
            GenerateDefaultCanvas(this);
            this.sketches.Focci = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',[0,0],... %Do not initalize with empty ("[]" ) because...
                'YData',[0,0],... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','none',...
                'MarkerFaceColor',this.color,...
                'MarkerEdgeColor',[0,0,0],...
                'Marker','*',...
                'DisplayName',[this.name, 'Focci'],...
                'Visible','off');
            this.generated.Focci = true;
        end%function
        function GenerateMajorRadius(this)
            GenerateDefaultCanvas(this);
            this.sketches.MajorRadius = quiver(...
                this.canvas,...
                0,... %XData
                0,... %YData
                0,... %UData
                0,... %VData
                'Color',[1,0,0],...
                'DisplayName',[this.name,'Major Radius'],...
                'AutoScale','off',...
                'Visible','off');
            this.generated.MajorRadius = true;
        end%function
        function GenerateMinorRadius(this)
            GenerateDefaultCanvas(this);
            this.sketches.MinorRadius = quiver(...
                this.canvas,...
                0,... %XData
                0,... %YData
                0,... %UData
                0,... %VData
                'Color',[0,0,1],...
                'DisplayName',[this.name,'Minor Radius'],...
                'AutoScale','off',...
                'Visible','off');
            this.generated.MinorRadius = true;
        end%function
        function GenerateVertexLabels(this)
            GenerateDefaultCanvas(this);
            this.sketches.VertexLabels = text(...
                this.canvas,...
                [0,0,0,0],...
                [0,0,0,0],...
                '');
            for ii = 1:4
                this.sketches.VertexLabels(ii).Visible = 'off';
                this.sketches.VertexLabels(ii).Interpreter = 'latex';
            end%ii
            this.generated.VertexLabels = true;
        end%function
        function GenerateFocciLabels(this)
            GenerateDefaultCanvas(this);
            this.sketches.FocciLabels = text(...
                this.canvas,...
                [0,0],...
                [0,0],...
                '');
            this.sketches.FocciLabels(1).Visible = 'off';
            this.sketches.FocciLabels(1).Interpreter = 'latex';
            this.sketches.FocciLabels(2).Visible = 'off';
            this.sketches.FocciLabels(2).Interpreter = 'latex';
        end%function
        function GenerateDefaultCanvas(this)
            if ~this.canvas_set
                warning('Default canvas created!');
                this.canvas = custom_axis;
                this.canvas_set = true;
            end%if
        end%function
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
        
        function [ax,elli] = TestLEEFillet
            ellipse.CleanSlate;
            %[ax,elli] = ellipse.TestLEEFillet
            ax = custom_axis;
            axis(ax,'equal');
            x1 = 1;
            y1 = 1;
            
            x2 = -5;
            y2 = -2;
            
            x3 = 1;
            y3 = -1;
            
            x4 = 3;
            y4 = 3;
            scatter(ax,[x1,x3],[y1,y3],'Marker','s')
            quiver(ax,[x1,x3],[y1,y3],[(x2-x1),(x4-x3)],[(y2-y1),(y4-y3)]);
            
            elli = ellipse.CreateLEEF(x1,y1,x2,y2,x3,y3,x4,y4);
            %elli = ellipse.CreateLEEF2(x1,y1,x2,y2,x3,y3,x4,y4);
            elli.SetCanvas(ax)
            elli.Toggle('Curve','MajorRadius','MinorRadius','Focci','Vertices')            
            
        end%function
        
    end%methods (Demonstrations)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Operations involving the quadric coefficients..
        function [A,B,C,D,E,F] = Metrics2Quadrics(xc,yc,a,b,eax,eay)
            A = a*a*eay*eay + b*b*eax*eax;
            B = 2*(b*b - a*a)*eax*eay;
            C = a*a*eax*eax + b*b*eay*eay;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc*xc + B*xc*yc + C*yc*yc - a*a*b*b;
        end%function
        function [xc,yc,a,b,eax,eay,ebx,eby] = Quadrics2Metrics(A,B,C,D,E,F)
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
            a = -sqrt(2*det4*(A + C  + sqrt((A - C)^2 + B*B)))/disc;
            b = -sqrt(2*det4*(A + C  - sqrt((A - C)^2 + B*B)))/disc;
            
            th = atan((C- A - sqrt((A - C)^2 + B*B))/B)*(B ~= 0) + pi/4*(B == 0 & A > C);
            eax = cos(th);
            eay = sin(th);
            ebx = -eay;
            eby = +eax;
        end%function
        function [a,eax,eay,b,ebx,eby,Vx1,Vy1,Vx2,Vy2] = RadiiFromDirections(f0x,f0y,f1x,f1y,f2x,f2y)
            %f0x: X-coordinate of center.
            %f0y: Y-coordinate of center.
            %f1x: X-component of some vector from center to ellipse.
            %f1y: Y-component of some vector from center to ellipse.
            %f2x: X-component of some other vector...
            %f2y: Y-component of some other vector...
            
            %Intermediate quantities.
            mag1 = f1x*f1x + f1y*f1y; %Square magnitudes.
            mag2 = f2x*f2x + f2y*f2y;
            dotP = f1x*f2x + f1y*f2y; %Dot product
            
            %Find a Vertex
            t0 = 0.5*acot(0.5*(mag1 - mag2)/dotP);
            cos_t0 = cos(t0);
            sin_t0 = sin(t0);
            Vx1 = f0x + f1x*cos_t0 + f2x*sin_t0;
            Vy1 = f0y + f1y*cos_t0 + f2y*sin_t0;
            
            %Find another Vertex
            t0 = t0 + 0.5*pi;
            cos_t0 = cos(t0);
            sin_t0 = sin(t0);
            Vx2 = f0x + f1x*cos_t0 + f2x*sin_t0;
            Vy2 = f0y + f1y*cos_t0 + f2y*sin_t0;
            
            dxa = Vx1 - f0x;
            dya = Vy1 - f0y;
            dxb = Vx2 - f0x;
            dyb = Vy2 - f0y;
            a = sqrt(dxa*dxa + dya*dya);
            b = sqrt(dxb*dxb + dyb*dyb);
            eax = dxa/a;
            eay = dya/a;
            ebx = dxb/b;
            eby = dyb/b;
            
            %In case the first vertex actually pointed to the minor axis.
            if b > a
                [Vx1,Vx2] = ellipse.swap(Vx1,Vx2);
                [Vy1,Vy2] = ellipse.swap(Vy1,Vy2);
                [eax,ebx] = ellipse.swap(eax,ebx);
                [eay,eby] = ellipse.swap(eay,eby);
                [a,b] = ellipse.swap(a,b);
            end%if    
        end%function
        
        %Sector Constructions
        function XY = PolarSector2P(XY,N,x1,y1,x2,y2,xc,yc,e,eax,eay,b,ebx,eby)
            %Given the ellipse's defining data and two points, will trace
            %the points towards the center of the ellipse and infer
            %subtended angles. This information is used to produce an
            %elliptical arc using a polar form. The details of that are
            %covered in "PolarSector2P"'s description.
            %{
            dx1 = xc - x1; 
            dy1 = yc - y1;
            dx2 = xc - x2;
            dy2 = yc - y2;
            mag_1 = sqrt(dx1*dx1 + dy1*dy1);
            mag_2 = sqrt(dx2*dx2 + dy2*dy2);
            cos_th1a = (dx1*eax + dy1*eay)/mag_1;
            cos_th2a = (dx2*eax + dy2*eay)/mag_2;
            
            if (ebx*dx1 + eby*dy1) <= 0
                th1 = acos(-cos_th1a);
            else 
                th1 = -acos(-cos_th1a);
            end%if
            
            if (ebx*dx2 + eby*dy2) <= 0
                th2 = acos(-cos_th2a);
            else 
                th2 = -acos(-cos_th2a);
            end%if     
            %}
            th1 = PolarAngleOfPoint(xc,yc,eax,eay,ebx,eby,x1,y1);
            th2 = PolarAngleOfPoint(xc,yc,eax,eay,ebx,eby,x2,y2);

            XY = PolarSector2A(XY,N,th1,th2,xc,yc,e,eax,eay,b,ebx,eby);
        end%function
        function XY = PolarSector2A(XY,N,th1,th2,xc,yc,e,eax,eay,b,ebx,eby)
            %Given the ellipse's defining data and two angles measured from
            %the oriented major radius, uses a polar form to generate an
            %elliptical sector/arc. The user specifies how many points to
            %generate via the parameter "N". The storage buffer must be
            %provided and must be large enough to store the "N" 2D points.
            d_theta = (th2 - th1)/(N-1);
            theta = 0;
            crossAB = eax*eby + eay*ebx;
            for ii = 1:N
                cos_thA = cos(theta);
                cos_thB = cos(pi/2 - theta); %Should I replace this with sine?
                R = ellipse.PolarRadius(b,e,cos_thA);
                XY(ii,1) = xc + R*(+eby*cos_thA - eay*cos_thB)/crossAB;
                XY(ii,2) = yc + R*(-ebx*cos_thA + eax*cos_thB)/crossAB;
                theta = theta + d_theta; %March the Angle
            end%ii
        end%function
        
        %Uselful parametrizations
        function R = PolarRadius(b,e,cos_th)
            R = b/sqrt(1 - (e*cos_th)^2);
        end%function
        function th = PolarAngleOfPoint(xc,yc,eax,eay,ebx,eby,x1,y1)
            %Given an ellipse's center and the direction cosines of its
            %major and minor radii, returns the polar angle of said point
            %as measured from the direction of the major radius.
            dx = xc - x1; 
            dy = yc - y1;
            mag = sqrt(dx*dx + dy*dy);
            cos_tha = (dx1*eax + dy1*eay)/mag;
            if (ebx*dx + eby*dy) <= 0
                th = acos(-cos_tha);
            else
                th = -acos(-cos_tha);
            end%if
        end%function
        function R = FocalRadius(a,e,cos_th)
            R = a*(1 - e*e)/(1 - cos_th);
        end%function
        
        %Elementary definition routines.
        function [xc,yc,eax,eay,ebx,eby,b,c] = FocciAndMajor2Metrics(Fx1,Fy1,Fx2,Fy2,a)
            dx = Fx2 - Fx1;
            dy = Fy2 - Fy1;
            xc = Fx1 + dx*0.5;
            yc = Fy1 + dy*0.5;
            %xc = 0.5*(Fx1 + Fx2); %Center is the average of the Focci.
            %yc = 0.5*(Fy1 + Fy2);
            c = 0.5*sqrt((Fx2 - Fx1)^2 + (Fy1 - Fy2)^2); %C is half the distance.
            eax = 0.5*dx/c;
            eay = 0.5*dy/c;
            b = sqrt(a*a - c*c);
            ebx = +eay;
            eby = -eax;
        end%function
        
    end%methods
    
    %Elementary low-level functions that are not unique in application to
    %the class.
    methods (Static)
        %Basic Geometry features.
        function [x,y,t1,t2] = Intersect2p2p(x1,y1,x2,y2,x3,y3,x4,y4)
            %Intersect two lines defined by two segments in two dimensions.
            [x,y,t1,t2] = ellipse.Intersect2p2d(...
                x1,y1,... %First point
                x2 - x1,y2 - y1,... %Direction at first point.
                x3,y3,... %Second point
                x4 - x3,y4 - y3); %Direction at second point.
        end%function
        function [x,y,t1,t2] = Intersect2p2d(x1,y1,dx1,dy1,x2,y2,dx2,dy2)
            %Intersect two points and 2 directions in two dimensions. 
            %WARNING: THIS DOS NOT CHECK IF THE LINES ARE PARALLEL!
            t2 = (dx1*(y2 - y1) - dy1*(x2 - x1))/(dx2*dy1 - dy2*dx1);
            t1 = (x2 - x1 + t2*dx2)/dx1;
            x = x1 + t1*dx1;
            y = y1 + t1*dy1;
        end%function

        %"Quality of Life" routines.
        function CleanSlate
            clc;
            clear;
            close all;
        end%function
        function [A,B] = swap(A,B)
            buffer = B;
            B = A;
            A = buffer;
            clear buffer;
        end%function
        function sign = Signum(x)
            %Branchless sign function.
            sign = (x > 0) - (x < 0);
        end%function       
    end%methods (Static)
    
end%classdef
%% polygon.m
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% A data structure used to generate and manipulate polygons in the 2D
% plane. Custom polygon generation routines, affine transformations, and a
% few nonlinear transformations are offered. Self-intersection detection is
% offered to help treat complex polygons.
%% Formulae
% Perimeter
%%
% $P = \sum_{i = 1}^{n}{\sqrt{(x_{i+1} - x_{i})^{2} + (y_{i+1} - y _{i})^{2}}}$
%%
% Shoelace Formula (Area and centroid measure)
%%
% $A_{i} = x_{i}y_{i+1} - x_{i+1}y_{i}$
%%
% $A = \frac{1}{2}\sum_{i=0}^{n-1}{A_{i}}$
%% 
% $X_{c}=\frac{1}{6A} \sum_{i=0}^{n-1}{(x_{i}+x_{i+1})A_{i}}$
%% 
% $Y_{c}=\frac{1}{6A} \sum_{i=0}^{n-1}{(y_{i}+y_{i+1})A_{i}}$
%%
% Laplacian Smoothing
%%
% $\hat{P} = P + \lambda \nabla^{2}(P)$
%% Class definition
classdef polygon < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the polygon's graphics as they appear on axes legend.
        color %Color of polygon's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the polygon.
        sketches %Structure containing handles to the object's graphics.
    end%properties
    %DEFINING DATA variables
    properties (SetAccess = private)
        sides %Number of sides (and coordinates).
        XY %XY coordinates (1st column is for x-values, 2nd for y-values.)
    end%properties
    %METRIC variables
    properties (SetAccess = protected)
        aabb %Axis-Aligned bounding box of this polygon.
        orientation %Whether the normals point inward (true) or outward (false).
        nxy %XY components of the normal vectors at each segment.
        area %Area enclosed by the polygon.
        perimeter %Perimeter subtended by polygon.
        simple %Flag to denote whether the polygon self-intersects.
        convex %Flag to denote whether the polygon is convex.
        regular %Flag to denote whether the polygon is regular.
        open %Flag to denote whether the polygon should be thought of as an open curve.
        xc %X-coordinate of the centroid.
        yc %Y-coordinate of the centroid.
    end    
    %FLAG and STATE variables.
    properties (Hidden = true)
        %Graphics-related flags.
        supress_updates %Whether the graphics should be updated.
        graphics_initialized %Whether "line" and or "quiver" have been generated.
        normals_computed %Whether the inward/outward normals have been computed.
        canvas_set %Whether an axes object has been initialized.
        AABB_present
        
        updated
        valid %Whether this instance has a valid definition.
    end%properties (Hidden)    
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = polygon(varargin)
            %Defining and metric data.
            this.XY = [];
            this.simple = [];
            this.convex = [];
            this.open = [];
            this.area = [];
            this.perimeter = [];
            this.xc = [];
            this.yc = [];
            
            %Graphics related.
            this.sketches = struct(...
                'Curve',[],...
                'Normals',[],...
                'Centroid',[],...
                'Vertex_Labels',[],...
                'Segment_Labels',[]);
            this.name = 'Polygon';
            this.color = [0,0,0];
            
            %Flag and state variables.
            this.graphics_initialized = false;
            this.normals_computed = false;
            this.orientation = [];
            this.supress_updates = false;
            this.AABB_present = (exist('AABB','file') == 2);
        end%function
        
        %Custom creation routines.
        function poly = CreateRegularByLength(x0,y0,sides,length)
            if sides < 3
                error('Need atleast 3 sides for a polygon.');
            end%if
            poly = polygon;
            poly.sides = sides;
            poly.XY = zeros(sides,2);
            
            %Can this be replaced with a "polar loop"?
            theta = 0;
            d_theta = 2*pi/sides;
            R = length/(2*sin(d_theta)); %This is the formula for a circle's segment.
            for ii = 1:sides
                poly.XY(ii,1) = x0 + R*cos(theta);
                poly.XY(ii,2) = y0 + R*sin(theta);
                theta = theta + d_theta;
            end%ii
            
            %Measure properties manually.
            poly.perimeter = sides*length;
            poly.area = 0.25*sides*length*length/cot(pi/sides); %Formula from circle.
            poly.xc = x0;
            poly.yc = y0;
            
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.ComputeNormals;
            
            %Set Flags
            poly.open = false;
            poly.orientation = true; %This is oriented inwards.
            poly.regular = true;
            poly.convex = true;
            poly.simple = true;
            poly.valid = true;
        end%function
        function poly = CreateRegularByRadius(R,sides,x0,y0)
            if sides < 3
                error(['A polygon needs at least three (3) sides! Input was:',num2str(sides)]);
            end%id
            if R == 0
                error('Cannot input a radius of zero.');
            end%if
            poly = polygon;
            poly.sides = sides;
            
            poly.XY = zeros(poly.sides,2);
            
            theta = 0;
            d_theta = 2*pi/sides;
            for ii = 1:poly.sides
                poly.XY(ii,1) = x0 + R*cos(theta);
                poly.XY(ii,2) = y0 + R*sin(theta);
                theta = theta + d_theta;
            end%ii
            
            %Set Flags
            poly.open = false;
            poly.regular = true;
            poly.convex = true;
            poly.orientation = true; %This is oriented inwards.
            poly.simple = true;
            poly.valid = true;
        end%function
        function poly = CreateSimpleStar(x0,y0,pairs,ri,ro)
            if pairs < 2
                error('Need at least two pairs for a simple star!');
            end%if
            poly = polygon; %Call constructor. 
            poly.sides = pairs*2;
            poly.XY = zeros(poly.sides,2); %Allocate memory for the coordinates.
            
            theta = 0;
            delta_theta = 2*pi/pairs;
            
            %Generate the outer vertices.
            for ii = 1:2:(2*pairs - 1)
                poly.XY(ii,1) = x0 + ro*cos(theta);
                poly.XY(ii,2) = y0 + ro*sin(theta);
                theta = theta + delta_theta;
            end%ii
            
            %Generate the inner vertices.
            theta = pi/pairs;
            for ii = 2:2:(2*pairs)
                poly.XY(ii,1) = x0 + ri*cos(theta);
                poly.XY(ii,2) = y0 + ri*sin(theta);
                theta = theta + delta_theta;
            end%ii
            poly.Measure;
             
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.SetName([num2str(poly.sides),'-sided Simple Star']);
            poly.ComputeNormals;

            %Set Flags
            poly.open = false;
            poly.regular = false; %In general, the simple star is not regular.
            poly.simple = true; %Simple stars do not self-intersect.
            poly.convex = false; %By definition, this is a concave polygon.
            poly.orientation = true; %These are oriented inward.
            poly.valid = true;
        end%function
        function poly = CreateComplexStar(x0,y0,sides,length,q)
            if mod(sides,2) == 0
                sides = sides + 1;
                warning('This routine works only for odd-sides polygons.')
            end%if
            if 2*q > sides
                q = 2;
                warning('The "q" value cannot exceed half the number of sides. Defaulting to q = 2.');
            end%if
            poly = polygon;
            poly.sides = sides;
            poly.xc = x0;
            poly.yc = y0;
            d_theta = 2*pi/sides;
            R = length/(2*sin(d_theta)); %This is the formula for a circle's segment.
            theta = 0;
            for ii = 1:sides
                poly.XY(ii,1) = x0 + R*cos(theta);
                poly.XY(ii,2) = y0 + R*sin(theta);
                theta = theta + q*d_theta;
                theta = theta*(theta <= 2*pi) + (theta - 2*pi)*(theta > 2*pi);
            end%ii
            poly.perimeter = sides*sqrt((poly.XY(2,1) - poly.XY(ii,1))^2 + (poly.XY(2,2)-poly.XY(1,2))^2);
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if    
            poly.SetName([num2str(poly.sides),'-sided Complex Star']);
            poly.ComputeNormals;
            
            %Set Flags
            poly.open = false;
            poly.regular = true; %Complex stars do have equal length sides.
            poly.simple = false; %Complex stars self intersect by definition.
            poly.convex = false; %Complex stars deviate from the convexity definition.
            poly.orientation = true; %These are oriented inward.
            poly.valid = true;
        end%function
        function poly = CreateFromOffset(progenitor,offset)
            if ~progenitor.valid
                error('Input polygon is flagged as invalid.');
            end%if
            poly = polygon;
            poly.sides = progenitor.sides;
            
            poly.XY = polygon.Offset(...
                progenitor.sides,... %Sides of the original polygon.
                progenitor.XY,...%Coordinates of original polygon.
                offset); %Thickness
            
            %{
            poly.XY = zeros(poly.sides,2);
            %First and last point are offset manually
            %First point
            [poly.XY(1,1),poly.XY(1,2)] = polygon.OffsetPointWithNormals(...
                        progenitor.XY(poly.sides,1),  progenitor.XY(poly.sides,2),...
                        progenitor.XY(1,1),           progenitor.XY(1,2),...
                        progenitor.XY(2,1),           progenitor.XY(2,2),...
                        progenitor.nxy(poly.sides,1), progenitor.nxy(poly.sides,2),...
                        progenitor.nxy(1,1),          progenitor.nxy(1,2),...
                        offset);
            %Last point
            [poly.XY(poly.sides,1),poly.XY(poly.sides,2)] = polygon.OffsetPointWithNormals(...
                        progenitor.XY(poly.sides-1,1), progenitor.XY(poly.sides-1,2),...
                        progenitor.XY(poly.sides  ,1), progenitor.XY(poly.sides  ,2),...
                        progenitor.XY(1,1),            progenitor.XY(1,2),...
                        progenitor.nxy(poly.sides-1,1),progenitor.nxy(poly.sides-1,2),...
                        progenitor.nxy(poly.sides,1),  progenitor.nxy(poly.sides,2),...
                        offset);
                    
            %The 2nd to 2nd to last points generalize with this loop.
            for ii = 2:(poly.sides - 1)
                iim1 = ii - 1;
                iip1 = ii + 1;
                [poly.XY(ii,1),poly.XY(ii,2)] = polygon.OffsetPointWithNormals(...
                        progenitor.XY(iim1,1), progenitor.XY(iim1,2),...
                        progenitor.XY(ii,1),   progenitor.XY(ii,2),...
                        progenitor.XY(iip1,1), progenitor.XY(iip1,2),...
                        progenitor.nxy(iim1,1),progenitor.nxy(iim1,2),...
                        progenitor.nxy(ii,1),  progenitor.nxy(ii,2),...
                        offset);
            end %ii
            %}
            

            
            
            
            poly.Measure;
            poly.regular = progenitor.regular;
            poly.valid = true;
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.ComputeNormals;
        end%function
        function poly = CreateGranulatedCopy(progenitor,granules)
            %Given a parent polygon, creates another instance of a polygon
            %with evenly spaced "granules" along the sides of the polygon. 
            poly = polygon;
            poly.sides = progenitor.sides + (progenitor.sides - 1)*granules;
            poly.XY = zeros(poly.sides,2);
            if progenitor.regular %Regular polygon's can skip some computations.
                
            else %Irregular polygons are more expensive to granularize.
                for kk = 1:2
                    %Use vomit-inducing branchless conditionals to avoid
                    %having to repeat the code in the "ii" loop. (*vomits*)
                    kkis1 = kk == 1;
                    kkis2 = kk == 2;
                    start = 1*kkis1 + progenitor.sides*kkis2;
                    step = 1*kkis1 + (1 - progenitor.sides)*kkis2;
                    finish = (progenitor.sides - 1)*kkis1 + (2)*kkis2; 
                    for ii = start:step:finish
                        iip1 = ii + step;
                        dx = progenitor.XY(iip1,1) - progenitor.XY(ii,1);
                        dy = progenitor.XY(iip1,2) - progenitor.XY(ii,2);
                        L = sqrt(dx*dx + dy*dy); %Length of the progenitor segment.
                        s = 0;
                        ds = L/granules;
                        for jj = 1:(1 + granules)
                            poly.XY()
                        end%jj
                    end%ii
                end%kk
            end%if
        end%function
        function poly = CreateRandom(sides)
            poly = polygon;
            poly.sides = sides;
            poly.XY = rand(sides,2);
            poly.valid = true;
            poly.ComputeNormals;
        end%function
        %Create a polygon from a list of XY-coordinates.
        %THIS REALLY NEEDS TO BE ABLE TO REMOVE REPEATED COORDINATES. AT
        %BARE MINIMUM IT SHOULD CHECK THAT NO TWO CONSECUTIVE COORDINATES
        %ARE IDENTICAL.
        function poly = CreateFromList(sides,XY)
            poly = polygon;
            if sides < 3
                error('Need at least 3 sides for a polygon.')
            end%if
            
            %Check the input list for consecutive repeated coordinates.
            %{
            stops = 0;
            stop = []; %IN THE FUTURE, REPLACE THIS WITH THE DOUBLY-LINKED LIST ONCE THAT MATURES.
            for ii = 1:(sides - 1)
                iip1 = ii + 1;
                if XY(ii,1) == XY(iip1,1) && XY(ii,2) == XY(iip1,2)
                    stops = stops + 1;
                    stop = [stop,ii];
                end%if
            end%ii
            %}
            
            %Initially, take the inputs for granted.
            
            
            if XY(1,1) == XY(sides,1)
                poly.sides = sides - 1;
                poly.XY = XY(1:end-1,:);
            else
                poly.sides = sides;
                poly.XY = XY;
            end%if
            
            %}
            %{
            %Need to ensure that the entries in the list are unique within
            %a certain tolerance.
            [X_sorted,X_perm] = sort(XY(:,1));
             Y_sorted = XY(X_perm,2);
             for ii = 1:sides
                 iip1 = ii + 1;
                 if X_sorted(ii) == X_sorted(iip1) && Y_sorted(ii) == Y_sorted(iip1)
                     
                 end%if
             end%ii
             
             clear X_sorted;
             clear Y_sorted;
             clear X_perm;
            %}
            
            
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.Measure;
            poly.regular = false;
            poly.valid = true;
            poly.ComputeNormals;
        end%function
    end%methods
    %High-level instance MODIFICATION and QUERY routines.
    methods
        %Redefinition routines (modification)
        function RedefineFromList(this,sides,XY_in)
            %If redefining reset the validity flags.
            this.XY = XY_in;
            this.sides = sides;
            this.valid = true;
            this.ComputeNormals;
            if this.graphics_initialized
                this.sketches.Curve.XData = XY_in(:,1);
                this.sketches.Curve.YData = XY_in(:,2);
            end%if
        end%function
        function RedefineAsRegular(this)
            length = this.perimeter/this.sides;
            
        end%function
        function RedefineAsSimpleStar(this)
            
            
        end%function
        function RedefineAsComplexStar(this)
            if mod(this.sides,2) == 0
                warning('Support for Complex star generation is limited to odd-sides polygons.');
                return;
            end%if
            
            length = this.perimeter/this.sides;
        end%function
        
        %Affine Transformations (modification)
        function Displace(this,dx,dy)
            %[x'] = [x] + [dx]
            %[y'] = [y] + [dy]
            if dx == 0 && dy == 0
                return;
            end%if
            
            for ii = 1:this.sides
                this.XY(ii,1) = this.XY(ii,1) + dx;
                this.XY(ii,2) = this.XY(ii,2) + dy;
            end%ii
            this.xc = this.xc + dx;
            this.yc = this.yc + dy;
            if this.graphics_initialized
                UpdateByTranslation(this,dx,dy);
            end%if
            if this.AABB_present
                 this.aabb.Displace([dx,dy]);
            end%if
        end%function
        function Rotate(this,angle)
            %[x'] = [+cos(TH),+sin(TH)][x]
            %[y'] = [-sin(TH),+cos(TH)][y]
            if angle == 0 || mod(angle,2*pi) == 0
                return;
            end%if
            
            %Rotate the polygon about its centroid.
            sin_th = sin(angle);
            cos_th = cos(angle);
            
            %Initialize trackers to redefine the bounding box. %MOVE BACK
            %TO UPDATE RAW.
            max_x = this.XY(1,1);
            min_x = this.XY(1,1);
            max_y = this.XY(1,2);
            min_y = this.XY(1,2);
            
            for ii = 1:this.sides
                shift_x = this.XY(ii,1) - this.xc; %Need an intermediate copy of the points.
                shift_y = this.XY(ii,2) - this.yc;
                this.XY(ii,1) = +shift_x*cos_th + shift_y*sin_th + this.xc;
                this.XY(ii,2) = -shift_x*sin_th + shift_y*cos_th + this.yc;
                
                old_nx = this.nxy(ii,1); %Need an intermediate copy of the points.
                old_ny = this.nxy(ii,2);
                this.nxy(ii,1) = +cos_th*old_nx + sin_th*old_ny;
                this.nxy(ii,2) = -sin_th*old_nx + cos_th*old_ny;
                
                %Branchless update of AABB dimension trackers.
                max_x = this.XY(ii,1)*(this.XY(ii,1) >= max_x) + max_x*(this.XY(ii,1) < max_x);
                min_x = this.XY(ii,1)*(this.XY(ii,1) <= min_x) + min_x*(this.XY(ii,1) > min_x);
                max_y = this.XY(ii,2)*(this.XY(ii,2) >= max_y) + max_y*(this.XY(ii,2) < max_y);
                min_y = this.XY(ii,2)*(this.XY(ii,2) <= min_y) + min_y*(this.XY(ii,2) > min_y);
            end%ii
            if this.AABB_present
                this.aabb.RedefinePoints(2,[min_x,min_y],[max_x,max_y]);
            end%if
            if this.graphics_initialized
                UpdateRaw(this);
            end%if
        end%function
        function Shear(this,Sx,Sy)
            %[x'] = [ 1, Sy][x]
            %[y'] = [Sx,  1][y]
            if Sx == 0 && Sy == 0
                return;
            end%if
            
            for ii = 1:this.sides
                %Update the points on the curve.
                old_x = this.XY(ii,1) - this.xc;
                old_y = this.XY(ii,2) - this.yc;
                this.XY(ii,1) = old_x + old_y*Sy + this.xc;
                this.XY(ii,2) = old_y + Sx*old_x + this.yc;
            end%ii
            this.ComputeNormals;
            if this.graphics_initialized
                UpdateRaw(this);
            end%if
        end%function
        function Scale(this,Sx,Sy)
            %[x'] = [Sx,  0][x]
            %[y'] = [ 0, Sy][y]
            if Sx == 1 && Sy == 1
                return;
            end%if
            for ii = 1:this.sides
                %Update the points on the curve.
                this.XY(ii,1) = (this.XY(ii,1) - this.xc)*Sx + this.xc;
                this.XY(ii,2) = (this.XY(ii,2) - this.yc)*Sy + this.yc;
            end%ii
            this.ComputeNormals;
            if this.graphics_initialized
                UpdateRaw(this);
            end%if
        end%function

        %Custom Transformations (modification)
        function Smooth(this,lambda,N)
            if ~this.valid
                error('Polygon is flagged as invalid!');
            end%if
            if ~this.simple
                warning('This routine was not meant to be used with complex polygons!')
            end%if
            if ~this.open
                this.SmoothClose(lambda,N);
            else
                this.SmoothOpen(lambda,N);
            end%if
        end%function
        function Disperse(this,lambda,N)
        end%function
        
        
        %Move to low-level
        function SmoothClose(this,lambda,N)
            %WIP: THIS FUNCTION HAS A LOT OF COPIED CODE. THE LOW LEVEL
            %FUNCTIONS THAT POWER IT NEED TO BE WRITTEN.
            
            %Applies a single pass of discrete Laplacian smoothing to the
            %polygon's coordinates inplace. "N" is the number of stencil
            %points used in the scheme. It must be an odd integer equal to
            %three or greater because the polygon is assumed to b closed.
            
            
            %Need to make copies of the first and last "N/2" points in the
            %sequence of coordinates.
            if this.open
                warning('This routine does not support "open" polygons.')
                return;
            end%if
            %If # of stencil points is not specified, default to 3.
            if nargin < 3
                N = 3;
            end%if
            switch N
                case 3
                    FD = [+1,-2,+1];
                case 5
                    FD = [-1/12,+4/3,-2.5,+4/3,-1/12];
                case 7
                    FD = [+1/90,-3/20,+1.5,-49/18,+1.5,-3/20,+1/90];
                case 9
                    FD = [-1/560,+8/315,-0.2,+1.6,-205/72,+1.6,-0.2,+8/315,-1/560];
                otherwise
                   if exist('polynomial.m','file') == 2
                       %For schema involving more than 9 stencil points,
                       %one needs to built the Lagrange interpolation
                       %sequence to compute the coefficients.
                       sequence = zeros(1,N);
                       sequence(1) = -N;
                       for ii = 2:N
                           sequence(ii) = sequence(ii-1) + 1;
                       end%ii
                       FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
                   else
                       warning('Cannot implement finite difference scheme with more than 9 stencil points without the file "polynomial.m"');
                   end%if
            end%switch
            
            Nm1o2 = (N-1)*0.5; %Nm1o2 = "N minus 1 over 2"
            Np1o2 = Nm1o2 + 1; %This index corresponds to the middle of the central difference scheme.
            
            %Allocate memory for the "remember" buffers and initialize
            %them.
            ahead_x = zeros(1,Nm1o2); %Stores the first (N-1)/2...
            ahead_y = zeros(1,Nm1o2); %... stencil points.
            trail_x  = zeros(1,Nm1o2); %Stores the trailing (N-1)/2...
            trail_y  = zeros(1,Nm1o2); %... stencil points.
            for ii = 1:Nm1o2
                ahead_x(ii) = this.XY(ii,1);
                ahead_y(ii) = this.XY(ii,2);
                trail_x(ii) = this.XY(this.sides - Nm1o2 + ii,1);
                trail_y(ii) = this.XY(this.sides - Nm1o2 + ii,2);
            end%ii
            %Begin smoothing. In this first loop, only the trailing
            %coefficients are referenced.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP#1
            spanner = 1;
            sm1 = spanner - 1;
            for ii = 1:(this.sides - Nm1o2)
                x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
                y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
                
                %Contribution of the current point.
                x_diffusion = x_diffusion + this.XY(ii,1)*FD(Np1o2);
                y_diffusion = y_diffusion + this.XY(ii,2)*FD(Np1o2);
                
                %for
                %
                for jj = 1:Nm1o2
                    %Account for the permutation of the trailing points.
                    t_idx = jj + sm1; %Index into the trailing buffer.
                    t_idx = t_idx - Nm1o2*(Nm1o2 - t_idx < 0); %Push it back if it overflows.
                    %Contribution of trailing points and points ahead.
                    x_diffusion = x_diffusion + trail_x(t_idx)*FD(jj) + this.XY(ii + jj,1)*FD(Np1o2 + jj);
                    y_diffusion = y_diffusion + trail_y(t_idx)*FD(jj) + this.XY(ii + jj,2)*FD(Np1o2 + jj);
                end%jj
                
                %Update the permuted buffer before expending the original
                %point.
                trail_x(spanner) = this.XY(ii,1);
                trail_y(spanner) = this.XY(ii,2);
                spanner = (spanner + 1)*(spanner < Nm1o2) + 1*(spanner >= Nm1o2);
                sm1 = spanner - 1;

                %Apply the smoothing operation InPlace (spend the original
                %point).
                this.XY(ii,1) = this.XY(ii,1) + lambda*x_diffusion;
                this.XY(ii,2) = this.XY(ii,2) + lambda*y_diffusion;        
            end%ii
            
            %This second loop references both the trailing and ahead
            %coefficients.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP#2
            for ii = (this.sides - Nm1o2 + 1):(this.sides - 1)
            %for ii = (this.sides - Nm1o2 + 1):(this.sides)
                x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
                y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
                
                %Contribution of the current point
                x_diffusion = x_diffusion + this.XY(ii,1)*FD(Np1o2);
                y_diffusion = y_diffusion + this.XY(ii,2)*FD(Np1o2);
                          
                %Contribution of the trailing points (be wary of permutation).
                for jj = 1:Nm1o2
                    t_idx = jj + sm1; %Index into the trailing buffer.
                    t_idx = t_idx - Nm1o2*(Nm1o2 - t_idx < 0); %Push it back if it overflows.
                    x_diffusion = x_diffusion + trail_x(t_idx)*FD(jj);
                    y_diffusion = y_diffusion + trail_y(t_idx)*FD(jj);
                end%jj
                %Contribution of the points ahead (not in the "ahead"
                %buffer).
                for jj = 1:(this.sides - ii)
                    x_diffusion = x_diffusion + this.XY(ii+jj,1)*FD(Np1o2 + jj);
                    y_diffusion = y_diffusion + this.XY(ii+jj,2)*FD(Np1o2 + jj);
                end%jj
                %Contribution of the points ahead (in the "ahead" buffer).
                for jj = (this.sides - ii+1):Nm1o2
                    x_diffusion = x_diffusion + ahead_x(jj)*FD(Np1o2 + jj);
                    y_diffusion = y_diffusion + ahead_y(jj)*FD(Np1o2 + jj);
                end%jj
                
                
                %Update the permuted buffer before expending the original
                %point.
                trail_x(spanner) = this.XY(ii,1);
                trail_y(spanner) = this.XY(ii,2);
                spanner = (spanner + 1)*(spanner < Nm1o2) + 1*(spanner >= Nm1o2);
                sm1 = spanner - 1;
                
                this.XY(ii,1) = this.XY(ii,1) + lambda*x_diffusion;
                this.XY(ii,2) = this.XY(ii,2) + lambda*y_diffusion; 
            end%ii
            
            %This third (and final loop is just for the end point.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP#3
            x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
            y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
            
            %Contribution of the current point
            x_diffusion = x_diffusion + this.XY(this.sides,1)*FD(Np1o2);
            y_diffusion = y_diffusion + this.XY(this.sides,2)*FD(Np1o2);
            for jj = 1:Nm1o2
                t_idx = jj + sm1; %Index into the trailing buffer.
                t_idx = t_idx - Nm1o2*(Nm1o2 - t_idx < 0); %Push it back if it overflows.
                x_diffusion = x_diffusion + trail_x(t_idx)*FD(jj) + ahead_x(jj)*FD(Np1o2 + jj);
                y_diffusion = y_diffusion + trail_y(t_idx)*FD(jj) + ahead_y(jj)*FD(Np1o2 + jj);
            end%jj
            %Apply the smoothing operation InPlace
            this.XY(this.sides,1) = this.XY(this.sides,1) + lambda*x_diffusion;
            this.XY(this.sides,2) = this.XY(this.sides,2) + lambda*y_diffusion;
           
            clear ahead_x; %free
            clear ahead_y; %free
            clear trail_x; %free
            clear trail_y; %free
            
   
            this.Measure;
            this.ComputeNormals; %The distortion changes the normals.
            if this.AABB_present
                this.aabb.RedefineFromList(this.sides,this.XY);
            end%if
            if this.graphics_initialized
                UpdateRaw(this);
            end%if
        end%
        function SmoothOpen(this,lambda,N)
            if this.open
                warning('This routine does not support "closed" polygons.')
                return;
            end%if
            %If # of stencil points is not specified, default to 3.
            if nargin < 3
                N = 3;
            end%if
            if ~(exist('polynomial.m','file') == 2)
                error('This routine requires the file "polynomial.m"');
            end%if
            sequence = zeros(1,N);%Initialize a sequence buffer.
            
            %First, the finite difference scheme is applied to the left
            %points of the open curve.
            for ii = 1:N
                sequence(ii) = ii - 1;
            end%ii
            FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
            
            %The next step is to smooth all the points that are elligible
            %for a central difference.
            Nm1o2 = (N-1)*0.5; %Nm1o2 = "N minus 1 over 2"
            Np1o2 = Nm1o2 + 1; %This index corresponds to the middle of the central difference scheme.
            ahead_x = zeros(1,Nm1o2); %Stores the first (N-1)/2...
            ahead_y = zeros(1,Nm1o2); %... stencil points.
            trail_x  = zeros(1,Nm1o2); %Stores the trailing (N-1)/2...
            trail_y  = zeros(1,Nm1o2); %... stencil points.
            for ii = 1:Nm1o2
                ahead_x(ii) = this.XY(ii,1);
                ahead_y(ii) = this.XY(ii,2);
                trail_x(ii) = this.XY(this.sides - Nm1o2 + ii,1);
                trail_y(ii) = this.XY(this.sides - Nm1o2 + ii,2);
            end%ii
            
        end%function
        %Move to low-level
        
        function SmoothHighOrder(this,lambda,N)
            %Applies discrete Laplacian Smoothing to the polygon according
            %to a finite difference scheme for approximating the 2nd
            %derivative.
            if N < 3
                error(['The minimum # of points used is 3! input was',num2str(N)])
            end%if
            if mod(N,2) == 0 && N > 4
                N = N - 1;
                warning(['No central scheme uses an even number of stencil points. Input was:',num2str(N)]);
            end%if
            if N > this.sides
                N = this.sides; %This will at minimum be 3 (as per polygon rules).
                warning(['Polygon has less sides (',num2str(this.sides),') than points queried for smoothing (',num2str(N),')!']);
            end%if
            
            
            
            access2polynomial = (exist('polynomial.m','file') == 2);
            if this.open 
                if ~access2polynomial
                    warning('Cannot apply high-order smoothing on a "open" polygon without polynomial.m');
                    return;
                end%if
            end%if
            
            if N < 9
            elseif exist('polynomial.m','file') == 2
            end%if
            
        end%function
        
        function DetectSharpTurns(this,Rf)
            %This routine applies a "rolling" circle type of operation to
            %flag turns that exceed the radius of the tighest turn that the
            %polygon should make.
            bad_points = zeros(1,this.sides); %Allocate memory for an array of flags.
            %Sort the polygon by Coordinates
            [X_sorted,X_perm] = sort(this.XY(:,1));
            
            %For this to work, the polygon must be oriented inwards.
            if ~this.orientation
                ReverseNormals(this);
            end%if
            
            %Need to first find the where the first point ranks in the
            %sorted list.
            needle = 1;
            while X_perm(needle) ~= 1
                needle = needle + 1;
            end%while
            old_needle = needle;

            %Initialize the circle's center
            %XC1 = this.XY(ii,1) + Rf*this.nxy(1,1);
            %YC1 = this.XY(ii,2) + Rf*this.nxy(ii,2);
            for ii = 1:(this.sides - 1)
                iip1 = ii + 1;
                XC1 = this.XY(ii,1) + Rf*this.nxy(ii,1);
                YC1 = this.XY(ii,2) + Rf*this.nxy(ii,2);
                XC2 = this.XY(iip1,1) + Rf*this.nxy(iip1,1);
                YC2 = this.XY(iip1,2) + Rf*this.nxy(iip1,2);
                dx = XC2 - XC1; %Change in the x-coordinate of the polygon as it rolled.
                dy = YC2 - YC1; %Change in the y-coordinate of the polygon as it rolled.
                
                %Inspect the sign of "dx"
                if dx < 0
                    while X_perm(needle) ~= ii
                        %needle = (needle - 1)*(needle > 1) + this.sides*(needle == 1);
                        needle = needle - 1;
                        if needle == 0 %Keep index within bounds.
                            needle = poly.sides;
                        end%if
                    end%while
                    old_needle = needle;
                    needleL = needle - 1; %Need a counter to subtract from the first needle.
                    while X_sorted(needleL) > (XC2 - Rf) %&& X_sorted(1) > (F.xc - F.R)
                        dxc = this.XY(X_perm(needleL),1) - XC2; %Distance from center to leftward point.
                        dyc = this.XY(X_perm(needleL),2) - yc2; %
                        bad_points(X_perm(needleL)) = ((dxc*dxc + dyc*dyc) < Rf); 
                        needleL = needleL -1;
                        if needleL == 0 %If even the minimum X-value is to the left, break.
                            %needleL = poly.sides;
                            break;
                        end%if
                    end%while
                end%if
                if dx == 0
                    
                end%if
                if dx > 0
                    while X_perm(needle) ~= ii
                        needle = needle + 1;
                        if needle == poly.sides + 1 %Keep index within bounds.
                            needle = 1;
                        end%if
                    end%while
                    old_needle = needle;
                end%if
                                
                
                
                
                
                
                
                needle = old_needle;
            end%ii
            
            
            
            
            
            
            %The algorithm for an open polygon is slightly different.
            if this.open
                
            else %Closed polygon.
                
                
                
                
                
                for kk = 1:3
                    kkis1 = kk == 1;
                    kkis2 = kk == 2;
                    kkis3 = kk == 3;
                    %{
                    start  = kkis1*this.sides;
                    step   = kkis1*(2 - this.sides);
                    finish = kkis1*2;
                    %}
                    
                    start  = kkis1*1;
                    trail  = kkis1*(this.sides - 1);
                    ahead  = kkis1*1;
                    finish = kkis1*1;
                    
                    
                    for ii = start:head:finish
                        iim1 = start + trail;
                        iip1 = start + ahead;
                        
                    end%ii
                end%kk
            end%if
            
            
        end%function
        
        function AddNoise(this,mag)
            for ii = 1:this.sides
                this.XY(ii,1) = this.XY(ii,1) + rand*mag;
                this.XY(ii,2) = this.XY(ii,2) + rand*mag;
            end%ii
            this.Measure;
            this.ComputeNormals; %The distortion changes the normals.
            if this.AABB_present
                this.aabb.RedefineFromList(this.sides,this.XY)
            end%if
            if this.graphics_initialized
                UpdateRaw(this);
            end%if
        end%function
            
        %Orientation
        function ComputeNormals(this)
            if ~this.valid
                warning('This polygon is flagged as invalid.');
                return;
            end%if
            if isempty(this.XY)
                warning('Polygon has no associated coordinates!');
                return;
            end%if
            this.nxy = polygon.Normals2D(this.sides,this.XY,this.nxy); %Magic happens here.
            if this.graphics_initialized %Update the quiver plot.
                this.sketches.Normals.UData = this.nxy(:,1);
                this.sketches.Normals.VData = this.nxy(:,2);
            end%if
            this.normals_computed = true;
        end%function
        function ReverseNormals(this)
            if ~this.normals_computed
                warning('Normal vectors NOT computed. Cannot reverse.')
                return;
            end%if
            for ii = this.sides
                this.nxy(ii,1) = (-1)*this.nxy(ii,1);
                this.nxy(ii,2) = (-1)*this.nxy(ii,2);
            end%ii
            this.orientation = (-1)*this.orientation;
            if this.graphics_initialized
                for ii = 1:this.sides
                    this.sketches.Normals.UData(ii) = (-1)*this.sketches.Normals.UData(ii);
                    this.sketches.Normals.VData(ii) = (-1)*this.sketches.Normals.VData(ii);
                end%ii
            end%if
        end%function
        
        %Measure metric properties about this polygon
        function Measure(this)
            %Quantifyable metrics.
            this.xc = 0;
            this.yc =0;
            this.perimeter = 0; %Initialize a running sum for the perimeter.
            this.area = 0; %Initialize a running sum for the area
            
            min_x = this.XY(1,1);
            min_y = this.XY(1,2);
            sign_xchange = 0;
            sign_ychange = 0;
            sign_xold = 0;
            sign_yold = 0;
            sign_xnew = 0;
            sign_ynew = 0;
            
            for ii = 1:(this.sides - 1)
                iip1 = ii + 1;
                
                %Shoelace formula for area and centroids.
                dA = this.XY(ii,1)*this.XY(iip1,2) - this.XY(iip1,1)*this.XY(ii,2);
                this.area = this.area + dA;
                this.xc = this.xc + dA*(this.XY(iip1,1) + this.XY(ii,1));
                this.yc = this.yc + dA*(this.XY(iip1,2) + this.XY(ii,2));
                
                %Perimeter.
                dx = this.XY(iip1,1) - this.XY(ii,1); %Change in x-coordinate.
                dy = this.XY(iip1,2) - this.XY(ii,2); %Change in y-coordinate.
                this.perimeter = this.perimeter + sqrt(dx*dx + dy*dy);
                
                %Used for determination of convexity.
                sign_xnew = (dx >= 0) - (dx < 0);
                sign_ynew = (dy >= 0) - (dy < 0);
                sign_xchange = sign_xchange + (sign_xnew ~= sign_xold);
                sign_ychange = sign_ychange + (sign_ynew ~= sign_yold);
                sign_xold = sign_xnew;
                sign_yold = sign_ynew;
                min_x = min_x*(min_x <= this.XY(ii,1)) + this.XY(ii,1)*(min_x > this.XY(ii,1));
                min_y = min_y*(min_y <= this.XY(ii,2)) + this.XY(ii,2)*(min_y > this.XY(ii,2));
            end%ii
            
            %Area and Centroid: Last segment contribution.
            dA = this.XY(this.sides,1)*this.XY(1,2) - this.XY(1,1)*this.XY(this.sides,2);
            this.area = this.area + dA;
            this.xc = this.xc + dA*(this.XY(1,1) + this.XY(this.sides,1));
            this.yc = this.yc + dA*(this.XY(1,2) + this.XY(this.sides,2));
            
            %Area and Centroid: Correction factors.
            this.area = this.area*0.5;
            this.xc = this.xc/(6*this.area);
            this.yc = this.yc/(6*this.area);
            
            %Perimeter: Last segment contribution.
            dx = this.XY(1,1) - this.XY(this.sides,1); %Change in x-coordinate.
            dy = this.XY(1,2) - this.XY(this.sides,2); %Change in y-coordinate.
            this.perimeter = this.perimeter + sqrt(dx*dx + dy*dy);
            
            this.orientation = (this.area >= 0) - (this.area < 0);%If signed area is positive, polygon is oriented inwards.
            this.area = this.area*((this.area > 0) - (this.area < 0));%Branchless way of taking absolute value.
            
            %Convexity inference: If the number of sign changes in the
            %X-coordinate plus the number of sign changes in the
            %Y-coordinate is exactly 2 then the polygon is convex (assuming
            %that the first point of the polygon has the minimum X and Y
            %values for coordinates.
            this.convex = (sign_xchange + sign_ychange - (min_x ~= this.XY(1,1)) - (min_y ~= this.XY(1,2)) == 2);
        end%function
        function kappa = Curvature(this)
            kappa = zeros(this.sides,1);
            last_x = this.XY(this.sides,1);
            last_y = this.XY(this.sides,2);
            for kk = 1:2
                kkis1 = kk == 1;
                kkis2 = kk == 2;
                start = 1*kkis1 + this.sides*kkis2;
                step = 1*kkis1 + (1 - this.sides)*kkis2;
                finish = (this.sides - 1)*kkis1 + 2*kkis2;
                for ii = start:step:finish
                    iip1 = ii + step;
                    next_x = this.XY(iip1,1);
                    next_y = this.XY(iip1,2);

                    %fit 3-point circle.
                    [XC,YC] = circle.Center3P(...
                        last_x,last_y,...
                        this.XY(ii,1),this.XY(ii,2),...
                        next_x,next_y);
                    
                    %Compute the curvature of the circle.
                    kappa(ii) = 1/sqrt((last_x - XC)^2 + (last_y - YC)^2);
                    
                    %Update the buffer to remember "old" values.
                    last_x = this.XY(ii,1);
                    last_y = this.XY(ii,2);
                end%ii
            end%kk
        end%function
        
        %Is point inside?
        function inside = IsInside(this,XY_in)
        end%function
        
        %Authenticator: (Corrects as much as it can before having to
        %error).
        function authenticate(this)
            this.valid = isempty(this.xc) + isempty(this.yc) + isempty(this.XY);
            this.valid = ~this.valid;
        end%function
        
        %Reset function: (Sets all validity-related flags to false. This is
        %to be called at the beginning of functions that redefine an
        %already-existing data structure).
        function reset(this)
            this.normals_computed = false;
            this.valid = false;
            this.simple = [];
            this.convex = [];
            this.orientation = [];
        end%function.
        
        
    end%methods
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
            
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.Show;
            end%if
            
        end%function
        function stamp = Imprint(this)
            if ~this.canvas_set
                warning([this.name, ' does not have a canvas set. Cannot stamp!']);
                return;
            end%if
            stamp = line(...
                'Parent',this.canvas,...
                'XData',[this.XY(:,1);this.XY(1,1)],...
                'YData',[this.XY(:,2);this.XY(1,2)],...
                'Color',this.color,...
                'LineStyle',this.sketches.Curve.LineStyle,...
                'DisplayName',['Stamp of ',this.name]);
        end%function
        
        %Graphical setting functions
        function SetName(this,string)
            this.name = string;
            if this.graphics_initialized
                this.sketches.Curve.DisplayName = string;
                this.sketches.Normals.DisplayName = [string, 'Normals'];
                this.sketches.Centroid.DisplayName = [string, 'Centroid'];
            end%if
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.SetName(string)
            end%if
        end%function        
        function SetColor(this,RGB)
            this.color = RGB;
            if this.graphics_initialized
                this.sketches.Curve.Color = RGB;
                this.sketches.Normals.Color = RGB;
                this.sketches.Centroid.MarkerFaceColor = RGB;
                for ii = 1:this.sides
                    this.sketches.Vertex_Labels(ii).Color = RGB;
                    this.sketches.Segment_Labels(ii).Color = RGB;
                end%ii
            end%if
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.SetColor(RGB);
            end%if
            
            
        end%function
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.AABB_present
                this.aabb.SetCanvas(ax)
            end%if
            if this.graphics_initialized
                %If AABB.m is in the same folder as polygon.m
                if this.AABB_present
                    this.aabb.sketches.Box.Parent = ax;
                    this.aabb.sketches.Center.Parent = ax;
                end%if
                this.sketches.Curve.Parent = ax;
                this.sketches.Centroid.Parent = ax;
                this.sketches.Normals.Parent = ax;
                for ii = 1:this.sides
                    this.sketches.Vertex_Labels(ii) = ax;
                    this.sketches.Segment_Labels(ii) = ax;
                end%ii
            end%if
        end%function
        
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
         
            %Draw the polygon's centroid.
            this.sketches.Centroid = line(...
                this.canvas,...
                this.xc,...
                this.yc,...
                'Marker','x',...
                'MarkerFaceColor',this.color,...
                'MarkerEdgeColor',[0,0,0],...
                'Visible','off',...
                'DisplayName',[this.name,' centroid']);
            
            %Draw the normal vectors at the midpoint of the edges of the
            %polygon.
            this.sketches.Normals = quiver(...
                this.canvas,...
                zeros(this.sides,1),...
                zeros(this.sides,1),...
                zeros(this.sides,1),... %Normals must be precomputed prior to... 
                zeros(this.sides,1),... %... drawing this quiver object.
                'AutoScale','on',...
                'Color',this.color,...
                'Visible','off',...
                'DisplayName',[this.name,' Orientation']);            
            if this.normals_computed
                %The positions of the arrows are then computed.
                for ii = 1:(this.sides - 1)
                    iip1 = ii + 1;
                    this.sketches.Normals.XData(ii) = 0.5*(this.XY(iip1,1) + this.XY(ii,1));
                    this.sketches.Normals.YData(ii) = 0.5*(this.XY(iip1,2) + this.XY(ii,2));
                    this.sketches.Normals.UData(ii) = this.nxy(ii,1);
                    this.sketches.Normals.VData(ii) = this.nxy(ii,2);
                end%ii
                this.sketches.Normals.XData(this.sides) = 0.5*(this.XY(1,1) + this.XY(this.sides,1));
                this.sketches.Normals.YData(this.sides) = 0.5*(this.XY(1,2) + this.XY(this.sides,2));
                this.sketches.Normals.UData(this.sides) = this.nxy(this.sides,1);
                this.sketches.Normals.VData(this.sides) = this.nxy(this.sides,2);
            end%if
                   
            %Create the labels for Vertices
            this.sketches.Vertex_Labels = text(...
                this.XY(:,1),...
                this.XY(:,2),...
                '');
            
            for ii = 1:this.sides
                this.sketches.Vertex_Labels(ii).Color = this.color;
                this.sketches.Vertex_Labels(ii).Visible = 'off';
                this.sketches.Vertex_Labels(ii).Interpreter = 'latex';
                this.sketches.Vertex_Labels(ii).String = ['V$_{',num2str(ii),'}$'];
            end%ii
            
            %Create the segment labels.
            this.sketches.Segment_Labels = text(...
                zeros(this.sides,1),...
                zeros(this.sides,1),...
                '');
            for ii = 1:(this.sides - 1)
                iip1 = ii + 1;
                this.sketches.Segment_Labels(ii).Visible = 'off';
                this.sketches.Segment_Labels(ii).Interpreter = 'latex';
                this.sketches.Segment_Labels(ii).Color = this.color;
                this.sketches.Segment_Labels(ii).String = ['$S_{',num2str(ii),'}$'];
                this.sketches.Segment_Labels(ii).Position(1) = (this.XY(iip1,1) + this.XY(ii,1))*0.5;
                this.sketches.Segment_Labels(ii).Position(2) = (this.XY(iip1,2) + this.XY(ii,2))*0.5;
            end%ii
            this.sketches.Segment_Labels(this.sides).Visible = 'off';
            this.sketches.Segment_Labels(this.sides).Interpreter = 'latex';
            this.sketches.Segment_Labels(this.sides).String = ['S$_{',num2str(this.sides),'}$'];
            this.sketches.Segment_Labels(this.sides).Position(1) = (this.XY(1,1) + this.XY(this.sides,1))*0.5;
            this.sketches.Segment_Labels(this.sides).Position(2) = (this.XY(1,2) + this.XY(this.sides,2))*0.5;
            
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.InitializeGraphics;
                this.aabb.ToggleBox;
            end%if
            
            this.graphics_initialized  = true;
        end%function.           
        function TerminateGraphics(this)
            %Release system resources used to render graphics.
            delete(this.sketches.Curve);
            delete(this.sketches.Normals);
            delete(this.sketches.Centroid);
            for ii = this.sides:-1:1
                delete(this.sketches.Segment_Labels(ii));
                delete(this.sketches.Segment_Labels(ii));
            end%ii
            
            %Update graphics flag.
            this.graphics_initialized = false;
        end%function
        
        %Visibility toggling functions
        function Toggle(this,varargin)
            %NOTE: It is possible to avoid the if-else "ladder" by using
            %MATLAB's notation struct.(*string*). C-programming does not
            %allow this, so here I settle with the "ladder" until I find
            %something more reminiscent of what is allowed in C.
            %NOTE: GNU has created something called "gperf" which may help
            %with this issue when porting to C and OpenGL.
            for ii = 1:(nargin - 1)
                if strcmpi(varargin{ii},'Curve') == 1
                    VisibilityToggle(this,this.sketches.Curve);
                elseif strcmpi(varargin{ii},'Centroid') == 1
                    VisibilityToggle(this,this.sketches.Centroid);
                elseif strcmpi(varargin{ii},'Normals') == 1
                    VisibilityToggle(this,this.sketches.Normals);
                elseif strcmpi(varargin{ii},'Vertex_Labels') == 1
                    for jj = 1:this.sides
                        VisibilityToggle(this,this.sketches.Vertex_Labels(jj));
                    end%jj
                elseif strcmpi(varargin{ii},'Segment_Labels') == 1
                    for jj = 1:this.sides
                        VisibilityToggle(this,this.sketches.Segment_Labels(jj));
                    end%jj
                else
                    warning('Unrecognizable graphics option.');
                end%if
            end
        end%function
        function VisibilityToggle(this,handle)
            if ~this.graphics_initialized
                warning(['Cannot toggle ',handle.DisplayName,' of polygon "',this.name,'" (graphics uninitialized).']);
                return;
            end%if
            if ~isvalid(handle)
                warning(['Cannot toggle ',handle.DisplayName,' of polygon "',this.name,'" (deleted handle).']);
                return;
            end%if
            if strcmp(handle.Visible,'on') == 1
                handle.Visible = 'off';
            else
                handle.Visible = 'on';
            end%if
        end%function
        function ToggleAABB(this)
            if ~this.AABB_present
                warning('File AABB.m is not present in the same folder as polygon.m.')
                return;
            end%if
            this.aabb.ToggleBox;
        end%function
        function ToggleAABBCenter(this)
            if ~this.AABB_present
                warning('File AABB.m is not present in the same folder as polygon.m.')
                return;
            end%if
            this.aabb.ToggleCenter;
        end%function
        
        %Graphics updating functions.
        function UpdateByTranslation(this,dx,dy)
            %Recall that the plot needs a repeated point at the end to
            %display correctly. This is manually displaced here.
            this.sketches.Curve.XData(this.sides + 1) = this.sketches.Curve.XData(this.sides + 1) + dx;
            this.sketches.Curve.YData(this.sides + 1) = this.sketches.Curve.YData(this.sides + 1) + dy;
            
            %Displace the remainder of the graphics.
            for ii = 1:this.sides
                %Displace the curve.
                this.sketches.Curve.XData(ii) = this.sketches.Curve.XData(ii) + dx;
                this.sketches.Curve.YData(ii) = this.sketches.Curve.YData(ii) + dy;

                %Displace the normal vectors.
                this.sketches.Normals.XData(ii) = this.sketches.Normals.XData(ii) + dx;
                this.sketches.Normals.YData(ii) = this.sketches.Normals.YData(ii) + dy;
                
                %Displace the Vertex Labels
                this.sketches.Vertex_Labels(ii).Position(1) = this.sketches.Vertex_Labels(ii).Position(1) + dx;
                this.sketches.Vertex_Labels(ii).Position(2) = this.sketches.Vertex_Labels(ii).Position(2) + dy;

                %Displace the Segment Labels
                this.sketches.Segment_Labels(ii).Position(1) = this.sketches.Segment_Labels(ii).Position(1) + dx;
                this.sketches.Segment_Labels(ii).Position(2) = this.sketches.Segment_Labels(ii).Position(2) + dy;
            end%ii
            if this.AABB_present
                
            end%if
        end%function
        function UpdateRaw(this)
            %When ALL aspects of the polygon had to be recomputed
            
            %Update the centroid's sketch.
            this.sketches.Centroid.XData = this.xc;
            this.sketches.Centroid.YData = this.yc;
            
            for kk = 1:2
                %Branchless conditionals.
                kkis1 = kk == 1;
                kkis2 = kk == 2;
                start = 1*kkis1 + this.sides*kkis2;
                step = 1*kkis1 + (1 - this.sides)*kkis2;
                finish = (this.sides - 1)*kkis1 + 2*kkis2;
                
                %Actual loop.
                for ii = start:step:finish                    
                    %Update the curve.
                    this.sketches.Curve.XData(ii) = this.XY(ii,1);
                    this.sketches.Curve.YData(ii) = this.XY(ii,2);
                    
                    %Update the vertex_labels;
                    this.sketches.Vertex_Labels(ii).Position(1) = this.XY(ii,1);
                    this.sketches.Vertex_Labels(ii).Position(2) = this.XY(ii,2);
                    
                    %Update the X-coordinates of the normals.
                    this.sketches.Normals.UData(ii) = this.nxy(ii,1);
                    this.sketches.Normals.VData(ii) = this.nxy(ii,2);
                    
                    %Update the normals and the segment labels.
                    iip1 = ii + step;
                    this.sketches.Normals.XData(ii) = 0.5*(this.XY(iip1,1) + this.XY(ii,1));
                    this.sketches.Normals.YData(ii) = 0.5*(this.XY(iip1,2) + this.XY(ii,2));
                    this.sketches.Segment_Labels(ii).Position(1) = this.sketches.Normals.XData(ii);
                    this.sketches.Segment_Labels(ii).Position(2) = this.sketches.Normals.YData(ii);
                end%ii
            end%kk        
            
            %Need to update the last repeated point in the plotting handle.
            this.sketches.Curve.XData(this.sides + 1) = this.sketches.Curve.XData(1);
            this.sketches.Curve.YData(this.sides + 1) = this.sketches.Curve.YData(1);
            if this.AABB_present
                this.aabb.UpdateRaw;
            end%if
        end%function
        
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
        
        %Create instances of regular polygons, name them, and show their
        %normals.
        function [ax,polygons] = TestRegularPolygons
            %Clear slate of variables.
            clc
            clear
            close all
            ax = custom_axis;
            ax.Color = [0,0,0]; %Make background black for constrast.
            axis(ax,'equal');
            
            shapes = 4;
            names = {'Triangle','Square','Pentagon','Hexagon'};
            sides = [3,4,5,6];
            color = {[1,0,0],[0,1,0],[0,0,1],[1,1,0]};
            center = {[0,0],[0,3],[3,3],[3,0]};
            polygons = cell(shapes,1);
            for ii = 1:shapes
                polygons{ii} = polygon.CreateRegularByLength(...
                    center{ii}(1),...
                    center{ii}(2),...
                    sides(ii),...
                    2);
                polygons{ii}.SetCanvas(ax)
                polygons{ii}.SetName(names{ii});
                polygons{ii}.SetColor(color{ii});
                polygons{ii}.Show;
                polygons{ii}.ToggleAABB;

                %polygons{ii}.Toggle('Normals');
            end%ii
            title(ax,'Now showing polygons.');
            
            input('Next Feature?');
            title(ax,'Now showing polygon orientation.');
            for ii = 1:shapes
                polygons{ii}.Toggle('Normals');
            end%ii
            
            input('Next Feature?');
            title(ax,'Now showing polygon Axis-Aligned Bounding Boxes.');
            for ii = 1:shapes
                polygons{ii}.ToggleAABB;
            end%ii
            
            input('Next Feature?');
            title(ax,'Now reversing polygon''s orientations.');
            for ii = 1:shapes
                polygons{ii}.ReverseNormals;
            end%ii
            
            %{
            input('Next Feature?');
            title(ax,'Now displacing polygons.');
            
            %Spread the polygons out.
            ex = [-1,-1,+1,+1]/sqrt(2);
            ey = [-1,+1,+1,-1]/sqrt(2);
            steps = 100;
            ds = 1/steps;
            for jj = 1:steps
                for ii = 1:shapes
                    polygons{ii}.Displace(ex(ii)*ds,ey(ii)*ds);
                end%ii
                pause(0.00001);
            end%jj
            
            %Permute the positions of the polygons.
            ex = [ 0,+1, 0,-1]/sqrt(2);
            ey = [+1, 0,-1, 0]/sqrt(2);
            steps = 100;
            ds = 4/steps;
            for jj = 1:steps
                for ii = 1:shapes
                    polygons{ii}.Displace(ex(ii)*ds,ey(ii)*ds);
                end%ii
                pause(0.00001);
            end%jj
            %}
            input('Next Feature?');
            title(ax,'Now rotating polygons.');
            steps = 100;
            ds = 2*pi/steps;
            for jj = 1:steps
                for ii = 1:shapes
                    polygons{ii}.Rotate(ds);
                end%ii
                pause(0.00001);
            end%jj
            
            input('Conclude Demo?');

            
            %lgd = legend(ax,'show',...
            %    'location','northeastoutside');
            
            lgd = legend(ax,...
                [polygons{1}.sketches.Curve,...
                polygons{2}.sketches.Curve,...
                polygons{3}.sketches.Curve,...
                polygons{4}.sketches.Curve,...
                polygons{1}.sketches.Normals,...
                polygons{2}.sketches.Normals,...
                polygons{3}.sketches.Normals,...
                polygons{4}.sketches.Normals],...
                'location','northeastoutside');
            title(lgd,'LEGEND');
            lgd.Color = [1,1,1];
        end%function
        
        %Create instances of regular concave polygons and offset them.
        function [ax,poly] = TestStarPolygons
            %Clean slate of variables.
            clc;
            clear;
            close all;
            
            ax = custom_axis;
            axis(ax,'equal');
            x0 = 1;
            y0 = 1;
            pairs = 2;
            ri = 0.5;
            ro = 2;
            poly = polygon.CreateSimpleStar(x0,y0,pairs,ri,ro);
            poly.SetCanvas(ax);
            poly.Show;
        end%function
        
        %Demonstrate Measurement of metric properties.
        function [ax,poly] = TestMeasurements
            clc;
            clear;
            close all;
            ax = custom_axis;
            axis(ax,'equal')
            input('Ready?')

            load('fiber_copy.mat');
            
            fiber_number = 19;
            fiber = fiber_copy{fiber_number};
            fiber(:,1) = fiber(:,1) - min(fiber(:,1));
            fiber(:,2) = fiber(:,2) - min(fiber(:,2));
            %fiber(:,1) = flip(fiber(:,1));
            %fiber(:,2) = flip(fiber(:,2));

            
            poly = polygon.CreateFromList(length(fiber),fiber);
            %poly = polygon.CreateSimpleStar(0,0,10,3,4);
            
            %polygon.CreateSimpleStar(x0,y0,pairs,ri,ro)
            
            poly.SetCanvas(ax);
            poly.Show;
            
            Atxt = text(...
                ax,...
                NaN,...
                NaN,...
                'A',...
                'interpreter','latex');
            Ltxt = text(...
                ax,...
                NaN,...
                NaN,...
                'L',...
                'interpreter','latex');
            %Showcase area measurements.
            Ltot = 0;
            Atot = 0;
            Cx = 0;
            Cy = 0;
            
            segment = line(...
                ax,...
                [NaN,NaN],...
                [NaN,NaN],...
                'Color',[0,0,1],...
                'Marker','o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[0,0,1],...
                'LineStyle','-'); 
            
            triangle = line(...
                ax,...
                [0,NaN,NaN,0],...
                [0,NaN,NaN,0],...
                'LineStyle','-');
            
            centroid = line(...
                ax,...
                NaN,...
                NaN,...
                'Color',[0,0,0],...
                'Marker','s',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[0,0,0]);
            
            centroid_label = text(...
                ax,...
                NaN,...
                NaN,...
                'C',...
                'Interpreter','latex');
            lgd = legend(ax,[triangle,segment,centroid],...
                {'Triangle','Segment','Centroid'},...
                'Location','northeastoutside');
            
            
            for kk = 1:2
                kkis1 = kk == 1;
                kkis2 = kk == 2;
                start = 1*kkis1 + (poly.sides)*kkis2;
                step = 1*kkis1 + (1 - poly.sides)*kkis2;
                finish = (poly.sides - 1)*kkis1 + 2*kkis2;
                for ii = start:step:finish
                    iip1 = ii + step;
                    dx = poly.XY(iip1,1) - poly.XY(ii,1);
                    dy = poly.XY(iip1,2) - poly.XY(ii,2);
                    L = sqrt(dx*dx + dy*dy);
                    Ltot = Ltot + L;
                    
                    
                    %Update the length label.
                    Ltxt.Position(1) = 0.5*(poly.XY(ii,1) + poly.XY(iip1,1));
                    Ltxt.Position(2) = 0.5*(poly.XY(ii,2) + poly.XY(iip1,2));
                    Ltxt.String = ['L = ',num2str(L)];
                
                    %Update the leader lines.
                    segment.XData = [poly.XY(ii,1),poly.XY(iip1,1)];
                    segment.YData = [poly.XY(ii,2),poly.XY(iip1,2)];
                    
                    %Shoelace formula for area and centroids.
                    dA = poly.XY(ii,1)*poly.XY(iip1,2) - poly.XY(iip1,1)*poly.XY(ii,2);
                    Atot = Atot + dA;
                    Cx = Cx + dA*(poly.XY(iip1,1) + poly.XY(ii,1));
                    Cy = Cy + dA*(poly.XY(iip1,2) + poly.XY(ii,2));
                    
                    %Update Triangle render.
                    triangle.XData(2:3) = [poly.XY(ii,1),poly.XY(iip1,1)];
                    triangle.YData(2:3) = [poly.XY(ii,2),poly.XY(iip1,2)];
                    
                    %Update area text label.
                    Atxt.Position(1) = (poly.XY(iip1,1) + poly.XY(ii,1))/3;
                    Atxt.Position(2) = (poly.XY(iip1,2) + poly.XY(ii,2))/3;
                    Atxt.String = ['A = ',num2str(dA)];
                    
                    %Update centroid & label
                    centroid.XData = Cx;
                    centroid.YData = Cy;
                    centroid_label.Position(1) = Cx;
                    centroid_label.Position(2) = Cy;
                    centroid_label.String = ['(',num2str(Cx),',',num2str(Cy),')'];
                    
                    pause(0.005);
                end%ii
            end%kk
            Atot = 0.5*Atot;
            title(ax,['Signed Area = ',num2str(Atot)]);
            input('Sign Correction.')
            
            Cx = Cx/(6*Atot);
            Cy = Cy/(6*Atot);
            centroid.XData = Cx;
            centroid.YData = Cy;
            centroid_label.Position(1) = Cx;
            centroid_label.Position(2) = Cy;
            centroid_label.String = ['(',num2str(Cx),',',num2str(Cy),')'];
            
            poly.Measure;
            poly.Toggle('Normals');
        end%function
        
        %Test Laplacian smoothing on one of the Shaolin shards.
        
        %Create a regular polygon. Add randomized noise to its coordinates
        %and then apply several passes of Laplacian Smoothing.
        function [ax,poly1] = TestNoiseSmoothing
            clc;
            clear;
            close all;
            
            ax = custom_axis;
            ax.Color = [0,0,0];
            axis(ax,'equal');
            
            %Create a regular polygon of many sides and show it.
            poly1 = polygon.CreateRegularByLength(0,0,80,2);
            poly1.SetCanvas(ax);
            poly1.SetColor([0,1,0]);
            poly1.Show;
            %poly1.Toggle('Normals');
            
            set(ax,'XLim',[min(poly1.XY(:,1)),max(poly1.XY(:,1))]);
            set(ax,'YLim',[min(poly1.XY(:,2)),max(poly1.XY(:,2))]);
            %Add random noise to the regular polygon
            title(ax,'Proceed to distort the Polygon?');
            input('Add Noise?');
            distortions = 2;
            smoothings = 100;
            
            spacer = '----------';
            fprintf('%10s\t%10s\t%10s\t%10s \n','Area','Perimeter','Centroid-X','Centroid-Y');
            fprintf('%10s\t%10s\t%10s\t%10s \n',spacer,spacer,spacer,spacer);

            for ii = 1:distortions
                title(ax,['#Distortions = ',num2str(ii)]);
                poly1.AddNoise(0.75);
                fprintf('%10f\t%10f\t%10f\t%10f \n',poly1.area,poly1.perimeter,poly1.xc,poly1.yc);
                pause(0.5);
            end%ii
            
            title(ax,'Proceed to smooth the Polygon?');
            input('Smooth-out the polygon?');
            for ii = 1:smoothings
                title(ax,['#Smoothings = ',num2str(ii)]);
                %poly1.Smooth(0.5);
                poly1.Smooth(0.1,5);

                fprintf('%10f\t%10f\t%10f\t%10f \n',poly1.area,poly1.perimeter,poly1.xc,poly1.yc);
                pause(0.001);
            end%ii
        end%function
        
        function [ax,poly] = TestStamping
            clc;
            clear;
            close all;
            %This test demonstrates the stamping feature along with some of
            %the affine transformations. First, a regular polygon is
            %defined. It is first stamped, then, it will be subject to
            %concurrent translation and rotation transformations while it
            %generates stamps.
            
            ax = custom_axis;
            ax.XLim = [-2,6];
            ax.YLim = [-6,2];

            axis(ax,'equal');
            
            input('ready?')
            sides = 11;
            length = 2;
            %poly = polygon.CreateRegularByLength(0,0,sides,length);
            poly = polygon.CreateComplexStar(0,0,sides,length,3);
            poly.SetCanvas(ax)
            poly.SetColor([1,0,0]);
            poly.SetName([num2str(sides),'-sided polygon']);
            poly.Show;
            stamp = poly.Imprint;
            stamp.Color = [0,0,1];
            
            delta_x = 4;
            steps = 100;
            revs = 2;
            dx = delta_x/steps;
            d_theta = revs*2*pi/steps;
            for ii = 1:steps
                poly.Displace(dx,0);
                poly.Rotate(d_theta)
                pause(0.001);
            end%ii
            stamp = poly.Imprint;
            stamp.Color = [0,0,1];

            delta_y = -4;
            steps = 100;
            revs = 2;
            dy = delta_y/steps;
            d_theta = revs*2*pi/steps;
            for ii = 1:steps
                poly.Displace(0,dy);
                poly.Rotate(d_theta)
                pause(0.001);
            end%ii
            stamp = poly.Imprint;
            stamp.Color = [0,0,1];
            
            delta_x = -4;
            steps = 100;
            revs = 2;
            dx = delta_x/steps;
            d_theta = revs*2*pi/steps;
            for ii = 1:steps
                poly.Displace(dx,0);
                poly.Rotate(d_theta)
                pause(0.001);
            end%ii
            stamp = poly.Imprint;
            stamp.Color = [0,0,1];
        end%function
        
        function [ax,poly] = TestHighOrderSmoothing
            
        end%function.
        
    end%methods (Demonstrations)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Standalone computation of the polygon's area.
        function A = Area(sides,XY)
            A = 0; %Initialize running buffer for the area.
            for ii = 1:(sides - 1) %Loop computes area as if elements were quads.
                iip1 = ii + 1;
                A = A + XY(ii,1)*XY(iip1,2) - XY(iip1,1)*XY(ii,2);
            end%for            
            A = A/2; %Halve the area to compute as triangles.
        end%function
        
        %Standalone computation of the polygon's perimeter.
        function P = Perimeter(sides,XY)
            P = 0; %Initialize running counter for the perimeter.
            for ii = 1:(sides - 1) %This loop is valid for all sides but one.
                iip1 = ii + 1;
                P = P + sqrt((XY(iip1,1) - XY(ii,1))^2 +...
                    (XY(iip1,2) - XY(ii,2))^2);
            end%ii
            %Must compute the length of the segment between first and last
            %point manually.
            P = P + sqrt((XY(1,1) - XY(sides,1))^2 +...
                    (XY(1,2) - XY(sides,2))^2);
        end%function
        
        %Compute the Area and Perimeter simultaneously.
        function [A,P] = AreaPerimeter(sides,XY)
            %Since they share the same loop structure....
            P = 0; %Initialize running counter for the perimeter.
            A = 0; %Initialize running buffer for the area.
            for ii = 1:(sides - 1) %Loop computes area as if elements were quads.
                iip1 = ii + 1;
                
                P = P + sqrt((XY(iip1,1) - XY(ii,1))^2 +...
                    (XY(iip1,2) - XY(ii,2))^2);
                A = A + XY(ii,1)*XY(iip1,2) - XY(iip1,1)*XY(ii,2);
            end%for   
            A = A/2; %Halve the area to compute as triangles.
            P = P + sqrt((XY(1,1) - XY(sides,1))^2 +...
                    (XY(1,2) - XY(sides,2))^2);
        end%function
        
        %Self-Intersections
        function [SI,XY_SI,seg] = DetectSelfIntersections(sides,XY)
            if sides < 3
                error('Polygon with less than 3 sides?');
            end%if
            SI = 0;
            XY_SI = [];
            seg = [];
            if sides == 3
                return;
            end%if

            p_ii1 = zeros(1,2); %First point on segment "ii".
            p_ii2 = zeros(1,2); %Second point on segment "ii".
            p_jj1 = zeros(1,2); %First point on segment "jj".
            p_jj2 = zeros(1,2); %Second point on segment "jj".
            box_ii = AABB.Create2P(2,p_ii1,p_ii2); %Define a bounding box for segment "ii."
            box_jj = AABB.Create2P(2,p_jj1,p_jj2); %Define a bounding box for segment "jj."
            
            SI = 0;%Initialize counter for the self-intersections.
            for ii = 1:(sides - 3) %Check all sides
                iip1 = ii + 1;
                
                %Starting point of "ii^th" segment.
                p_ii1(1) = XY(ii,1);
                p_ii1(2) = XY(ii,2);
                
                %Ending point of "ii^th" segment.
                p_ii2(1) = XY(iip1,1);
                p_ii2(2) = XY(iip1,2);
                
                box_ii.RedefinePoints(2,p_ii1,p_ii2);
                
                for jj = (ii + 2):(sides-1) %Skip the segment immediately ahead (cannot intersect that one).
                    %fprintf('Checking %i against %i.\n',ii,jj)
                    jjp1 = jj + 1;
                    
                    p_jj1(1) = XY(jj,1); %Starting point of "ii^th" segment.
                    p_jj1(2) = XY(jj,2);
                    p_jj2(1) = XY(jjp1,1); %Ending point of "ii^th" segment.
                    p_jj2(2) = XY(jjp1,2);
                    box_jj.RedefinePoints(2,p_jj1,p_jj2);
                    
                    if AABB.BoxesOverlap(2,box_ii.p1,box_ii.p2,box_jj.p1,box_jj.p2)
                        fprintf('overlap\n')
                        [X,t1,t2] = intersect2p2p(p_ii1,p_ii2,p_jj1,p_jj2);
                        if 0 < t1 && t1 < 1 && 0 < t2 && t2 < 1
                            SI = SI + 1;
                            XY_SI = [XY_SI;[X(1),X(2)]];
                            seg = [seg; [ii,jj]];
                        end%if
                    end%if 
                end%jj
            end%ii
            
            %The above loop checks ALL edges against each other except for
            %the last one.
            p_ii1(1) = XY(sides,1);%Starting point of "last" segment.
            p_ii1(2) = XY(sides,2);
            p_ii2(1) = XY(1,1);%Ending point of "last" segment.
            p_ii2(2) = XY(1,2);
            box_ii.RedefinePoints(2,p_ii1,p_ii2);

            for jj = 2:sides-2
                %fprintf('Checking %i against %i.\n',sides,jj)
                jjp1 = jj + 1;
                p_jj1(1) = XY(jj,1); %Starting point of "ii^th" segment.
                p_jj1(2) = XY(jj,2);
                p_jj2(1) = XY(jjp1,1); %Ending point of "ii^th" segment.
                p_jj2(2) = XY(jjp1,2);
                box_jj.RedefinePoints(2,p_jj1,p_jj2);
                if AABB.BoxesOverlap(2,box_ii.p1,box_ii.p2,box_jj.p1,box_jj.p2)
                    [X,t1,t2] = intersect2p2p(p_ii1,p_ii2,p_jj1,p_jj2);
                    if 0 < t1 && t1 < 1 && 0 < t2 && t2 < 1
                        SI = SI + 1;
                        XY_SI = [XY_SI;[X(1),X(2)]];
                        seg = [seg; [sides,jj]];
                    end%if
                end%if
                
                
            end%ii
            %}
            
        end%function
        
        %Offset points by a certain thickness according to positive (O = 1)
        %or negative (O = -1) right-handed coordinate system.
        function XY_out = Offset(sides,XY,offset)
            XY_out = zeros(sides,2);
            
            %The first point will be done manually as it cannot be
            %generalized as part of the loop.
            [XY_out(1,1),XY_out(1,2)] = polygon.OffsetPoint(...
                XY(sides,1),XY(sides,2),... %(x1,y1)Point "to the left."
                XY(1,1)    ,XY(1,2),...     %(x2,y2)The point being offset.
                XY(2,1)    ,XY(2,2),...     %(x3,y3)Point "to the right."
                offset);            
            %The remaining sides are handled by identical code inside the
            %loop.
            for ii = 2:(sides - 1)
                iip1 = ii + 1;
                iim1 = ii - 1;
                [XY_out(ii,1),XY_out(ii,2)] = polygon.OffsetPoint(...
                    XY(iim1,1),XY(iim1,2),... %Point "to the left."
                    XY(ii,1)  ,XY(ii,2),... %The point being offset.
                    XY(iip1,1),XY(iip1,2),... %Point "to the right."
                    offset);
            end%ii
            
            %Need to do the last pair of segments manualy as well.
            [XY_out(sides,1),XY_out(sides,2)] = polygon.OffsetPoint(...
                XY(sides - 1,1),XY(sides - 1,2),... %(x1,y1) Point "to the left."
                XY(sides,1)    ,XY(sides,2),...     %(x2,y2) The point being offset.
                XY(1,1),XY(1,2),...                 %(x3,y3) Point "to the right."
                offset);
        end%function
        
        %Offset a single point. This computes the normals.
        function [X_off,Y_off] = OffsetPoint(x1,y1,x2,y2,x3,y3,offset)
            dx1 = x2 - x1;
            dy1 = y2 - y1;
            dx2 = x3 - x2;
            dy2 = y3 - y2;
            [nx1,ny1] = polygon.UnitNormal2D(dx1,dy1);
            [nx2,ny2] = polygon.UnitNormal2D(dx2,dy2);
            [X_off,Y_off,~,~] = polygon.Intersect2p2d(...
                x1 + nx1*offset,y1 + ny1*offset,... %First point.
                dx1,dy1,... %Direction at first point.
                x2 + nx2*offset,y2 + ny2*offset,... %Second point.
                dx2,dy2); %Direction at the second point.
        end%function
        
        function [X_off,Y_off] = OffsetPointWithNormals(x1,y1,x2,y2,x3,y3,nx1,ny1,nx2,ny2,offset)
            dx1 = x2 - x1;
            dy1 = y2 - y1;
            dx2 = x3 - x2;
            dy2 = y3 - y2;
            [X_off,Y_off,~,~] = polygon.Intersect2p2d(...
                x1 + nx1*offset,y1 + ny1*offset,... %First point.
                dx1,dy1,... %Direction at first point.
                x2 + nx2*offset,y2 + ny2*offset,... %Second point.
                dx2,dy2); %Direction at the second point.
        end%function
        
        %Apply a Laplacian smoothing to input XY coordinates
        function [X,Y] = AppendWeightedCentralDifference(N,FD,x,y,lbx,lby,rbx,rby,lambda)
            %"N" = # stencil points (including the point being edited.
            %"FD" = array of finite difference coefficients.
            %"x" = x-coordinate being smoothed.
            %"y" = y-coordinate being smoothed.
            %"lbx" = left buffer with (N-1)/2 entries of surrounding
            %x-coordinates.
            %"lby" = left buffer with (N-1)/2 entries of surrounding
            %y-coordinates.
            %"rbx" = left buffer with (N-1)/2 entries of surrounding
            %x-coordinates.
            %"rby" = left buffer with (N-1)/2 entries of surrounding
            %y-coordinates.
            %"lambda" = weighting factor
            
            
            
        end%function
        function [X,Y] = ApplyAsymetricDifference(sides,XY,N,FD,shift)
            %Same as "ApplyCentralDifference" but with the addition of the
            %"shift" variable. This is intended to apply forward and
            %backward differences. If the FD coefficients are "lopsided",
            %the user must specify the shift parameter correctly. For
            %backward difference, shift is "-(N+1)/2". Forward differences
            %should have shift = +(N+1)/2. Lopsided differences should
            %have: -(N+1)/2 < shift < +(N+1)/2.
           
        end%function
        
        %Compute the normal vectors of a polygon's sides. Needed mainly for
        %offseting operations.
        function nxy = Normals2D(sides,XY,nxy)
            if isempty(nxy) || nargin < 3
                nxy = zeros(sides,2);%Allocate output memory.
            end%if
            for kk = 1:2
                kkis1 = kk == 1;
                kkis2 = kk == 2;
                start = 1*kkis1 + sides*kkis2;
                step = 1*kkis1 + (1 - sides)*kkis2;
                finish = (sides - 1)*kkis1 + 2*kkis2;
                for ii = start:step:finish
                    iip1 = ii + step;
                    [nxy(ii,1),nxy(ii,2)] = polygon.UnitNormal2D(...
                        XY(iip1,1) - XY(ii,1),...
                        XY(iip1,2) - XY(ii,2));
                end%jj
            end%ii
        end%function
        
    end%methods (Static)
    %Elementary low-level functions that are not unique in application to
    %the class.
    methods (Static)
        
        %Find the unit normal vector to some direction according to a
        %right-handed coordinate system.
        function [nx,ny] = UnitNormal2D(dx,dy)
            mag = sqrt(dx*dx + dy*dy); %Magnitude of the input direction.
            A = dx/mag; %Normalized X-component of the direction.
            B = dy/mag; %Normalized Y-component of the direction.
            detA = A*A + B*B; %Determinant of the matrix.
            nx = -B/detA;
            ny = A/detA;
        end%function
        
        %Intersect two lines defined by two segments in two dimensions.
        function [x,y,t1,t2] = Intersect2p2p(x1,y1,x2,y2,x3,y3,x4,y4)
            dx12 = x2 - x1;
            dy12 = y2 - y1;
            dx34 = x4 - x3;
            dy34 = y4 - y3;
            [x,y,t1,t2] = polygon.Intersect2p2d(x1,y1,dx12,dy12,x3,y3,dx34,dy34);
        end%function
        
        %Intersect two points and 2 directions in two dimensions. This does
        %not check if the input directions are parallel!
        function [x,y,t1,t2] = Intersect2p2d(x1,y1,dx1,dy1,x2,y2,dx2,dy2)
            t2 = (dx1*(y2 - y1) - dy1*(x2 - x1))/(dx2*dy1 - dy2*dx1);
            t1 = (x2 - x1 + t2*dx2)/dx1;
            x = x1 + t1*dx1;
            y = y1 + t1*dy1;
        end%function
        
    end%methods (Static)
end%classdef
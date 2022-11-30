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
%%
% Discrete Curvature (Based on area)
%%
% $\kappa_{i} = \frac{2\alpha_{k}\pi}{L_{k-1} + L_{k}}$
%%
% Discrete Curvature (based on angles an 3-point circle)
%%
% $\kappa_{i} = \frac{1}{R} \sqrt{\frac{\alpha_{i} + \beta_{i}}{\sin{\alpha_{i}}+\sin{\beta_{i}}}}$
%% Class definition
classdef polygon < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        linewidth %(of the curve and last edge).
        name %Name of the polygon's graphics as they appear on axes legend.
        color %Color of polygon's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the polygon.
        sketches %Structure containing handles to the object's graphics.
        refresh_rate %A parameter that determines how long (ideally) calls to Refresh shoul take.
    end%properties
    %DEFINING DATA variables
    properties (SetAccess = private)
        sides %Number of sides (and coordinates).
        XY %XY coordinates (1st column is for x-values, 2nd for y-values.)
        open %Flag to denote whether the polygon should be thought of as an open curve.
    end%properties
    %METRIC variables
    properties (SetAccess = protected)
        aabb %Axis-Aligned bounding box of this polygon.
        orientation %Whether the normals point inward (true) or outward (false).
        nxy %XY components of the normal vectors at each segment. %DEPRECATE
        area %Area enclosed by the polygon (GLOBAL). 
        perimeter %Perimeter subtended by polyline (GLOBAL).
        simple %Flag to denote whether the polygon self-intersects (GLOBAL).
        convex %Flag to denote whether the polygon is convex (GLOBAL). 
        regular %Flag to denote whether the polygon is regular (GLOBAL). 
        xc %X-coordinate of the centroid.
        yc %Y-coordinate of the centroid. 
        SI %#of self-intersections.
        xSI %X-coordinates of self-intersections (SI by 1 array).
        ySI %Y-coordinates of self-intersections (SI by 1 array).
        sSI %Segments involved in self intersection (SI by 2 array).
        
        metL %"Local" vertex-wise metrics.
        
    end%properties    
    %FLAG and STATE variables.
    properties (Hidden = true)
        %Graphics-related flags.
        refresh %Structure holding boolean values that control which graphics refresh.
        updated %Structure holding boolean values that track graphics state of date.
        generated %Structure holding boolean values that track whether certain graphics are generated.
        canvas_set %Whether an axes object has been initialized.
        
        metL_alloc %Flags to indicate whether buffers for metL have been allocated.        
        
        normals_computed %Whether the inward/outward normals have been computed.
        AABB_present
        valid %Whether this instance has a valid definition.
    end%properties (Hidden)    
    %High-level instance CREATION routines.
    methods (Static)
        
        %Custom creation routines.
        function poly = CreateRegularByLength(x0,y0,sides,length)
            if sides < 3
                error('Need atleast 3 sides for a polygon.');
            end%if
            poly = polygon;
            poly.sides = sides;
            poly.Calloc(sides);
            poly.XY = CircularPolarLoop(...
                poly.XY,... %The polyline's buffer after calling "Calloc."
                poly.sides,... %The size of the buffer.
                length/(2*sin(d_theta)),... %This is the formula for a circle's segment.
                0,... %Start arc from 0 and...
                2*pi); %... end it at 2*pi to get the full circle.
            
            %Maybe replace this displacement loop with a function later on?
            for ii = 1:poly.sides
                poly.XY(ii,1) = poly.XY(ii,1) + x0; 
                poly.XY(ii,2) = poly.XY(ii,2) + y0;
            end%ii
            %Measure properties manually.
            poly.perimeter = poly.sides*length;
            poly.area = 0.25*sides*length*length/cot(pi/sides); %Formula from circle.
            poly.xc = x0;
            poly.yc = y0;
            
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            
            %Set Flags
            poly.open = false;
            poly.orientation = true; %This is oriented inwards.
            poly.regular = true;
            poly.convex = true;
            poly.simple = true;
            poly.valid = true;
        end%function
        function poly = CreateRegularByRadius(R,N,x0,y0,th1,th2)
            %Creates a polyline of "N" edges equal in length as implied by
            %a circular radius. The arc (open polyline) or circle (closed
            %polyline) are centered around the origin by default. The user
            %may specify custom centers by inputting "x0" and "y0". By
            %default, this routine creates circles. One may specify arcs by
            %inputting initial and final angles "th1" and "th2"
            %respectively.
            if R == 0
                error('Cannot input a radius of zero.');
            end%if
            %Default polar angle limits.
            if nargin < 5
                th1 = 0;
                th2 = 2*pi;
            end%if
            if th1 == 0 && th2 == 2*pi
                if N < 3
                    error(['A polygon needs at least three (3) sides! Input was:',num2str(N)]);
                end%id
                poly.open = false;
            else
                poly.open = true;
            end%if
            poly = polygon;
            poly.sides = N;
            poly.Calloc(N);
            %Create the circle centered around the origin.
            poly.XY = polygon.CircularPolarLoop(...
                poly.XY,... %Buffer where to store generated coordinates.
                poly.sides,... %Number of edges.
                R,... %Radius
                th1,... %Initial angle (radians).
                th2); %Final angle (radians).
            %Displace the coordinates so that the circle is centered around
            %user's origin.
            if nargin > 2
                poly.XY = polygon.DisplaceCoordinates(...
                    poly.XY,...
                    poly.sides,...
                    2,...
                    [x0,y0]);
            end%if
            %Set Flags
            poly.regular = true;
            poly.convex = true;
            poly.orientation = +1; %This is oriented inwards.
            poly.simple = true;
            poly.valid = true;
        end%function

        function poly = CreateSimpleStar(x0,y0,pairs,ri,ro)
            if pairs < 2
                error('Need at least two pairs for a simple star!');
            end%if
            poly = polygon; %Call constructor. 
            poly.sides = pairs*2;
            poly.Calloc(sides);
            
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
            poly.Calloc(sides);
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
            poly.Calloc(poly.sides);
            poly.XY = progenitor.Offset(offset);
            

            
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.Measure;
                        
            %Set flags.
            poly.regular = progenitor.regular;
            poly.open = progenitor.open;
            poly.valid = true;
        end%function
        function poly = CreateGranulatedCopy(progenitor,granules)
            %Given a parent polygon, creates another instance of a polygon
            %with evenly spaced "granules" along the sides of the polygon. 
            poly = polygon;
            poly.sides = progenitor.sides + (progenitor.sides - 1)*granules;
            poly.Calloc(poly.sides);
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
            polygon.Calloc(sides);
            for ii = 1:sides
                poly.XY(ii,1) = rand;
                poly.XY(ii,2) =rand;
            end%ii
            poly.open = false;
            poly.valid = true;
        end%function
        function poly = CreateFromList(points,XY,bool)
            %THIS REALLY NEEDS TO BE ABLE TO REMOVE REPEATED COORDINATES. AT
            %BARE MINIMUM IT SHOULD CHECK THAT NO TWO CONSECUTIVE COORDINATES
            %ARE IDENTICAL.
            
            %Create a polygon from a list of XY-coordinates. The user is
            %free to specify whether this is a closed or open discrete
            %curve. Certain operations like offsetting and smoothing change
            %their behavior depending on whether the curve is open or
            %closed. For open curves, set the "bool" flag to true. By
            %default, the routine assumes that the curves are closed. If
            %the flag is set to true but the XY list contains repeated
            %points at the beginning and end, the "bool" flag is overriden.
            poly = polygon;
            %Assume that it is closed by default.
            if nargin < 3
                bool = false;
            end%if
            poly.open = bool;
            %A repeated point will override bool.
            if XY(1,1) == XY(points,1) && XY(1,2) == XY(points,2)
                points = points - 1;
                poly.open = false;
            else
                poly.open = true;
            end%if
            if points < 3
                error('Need at least 3 sides for a polygon.')
            end%if
            poly.sides = points;
            poly.Calloc(poly.sides)
            for ii = 1:poly.sides
                poly.XY(ii,1) = XY(ii,1);
                poly.XY(ii,2) = XY(ii,2);
            end%ii
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            poly.Measure;
            poly.regular = false;
            poly.valid = true;
        end%function
    end%methods
    %High-level instance MODIFICATION and QUERY routines.
    methods
        %Constructor/Destructor
        function this = polygon(varargin)
            %Constructor
            this.XY = [];
            this.simple = [];
            this.convex = [];
            this.open = [];
            this.area = [];
            this.perimeter = [];
            this.xc = [];
            this.yc = [];
            
            this.metL = struct(...
                'areas',[],...
                'kappa',[],...
                'angles',[],...
                'lengths',[],...
                'nxy',[]);
            this.metL_alloc = struct(...
                'areas',false,...
                'kappa',false,...
                'angles',false,...
                'lengths',false,...
                'nxy',false);        
            
            %Self-Intersections
            this.SI = 0; %KEEP
            this.xSI = []; %KEEP
            this.ySI = []; %KEEP
            this.sSI = []; %KEEP
            
            %Graphics related.
            this.sketches = struct(...
                'Curve',[],...
                'LastEdge',[],...
                'Normals',[],...
                'Centroid',[],...
                'SelfIntersect',[],...
                'VertexLabels',[],...
                'SegmentLabels',[]);
            %{
            struct(...
                'Handle',[],...
                'Generated',false,...
                'Generator',@polygon.GenerateCurve,...
                'Updated',false,...
                'Refresh',false,...
                'Refresher',@polygon.RefreshCurve)
                %}
            this.refresh = struct(...
                'Curve',true,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',true,...
                'VertexLabels',false,...
                'SegmentLabels',false);
            this.updated = struct(...
                'Curve',false,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',false,...
                'VertexLabels',false,...
                'SegmentLabels',false);
            this.generated = struct(...
                'Curve',false,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',false,...
                'VertexLabels',false,...
                'SegmentLabels',false);
            
            %Default customization options.
            this.refresh_rate = 60;%Desired refresh rate in Hz;
            this.name = 'Polygon';
            this.color = [0,0,0];
            this.linewidth = 1;

            %Flag and state variables.
            this.normals_computed = false;
            this.orientation = [];
            this.AABB_present = (exist('AABB','file') == 2);
        end%function
        function delete(this)
            %Graphics Cleanup
            if this.generated.Curve
                delete(this.sketches.Curve);
                delete(this.sketches.LastEdge);
            end%if
            if this.generated.Normals
                delete(this.sketches.Normals);
            end%if
            if this.generated.Centroid
                delete(this.sketches.Centroid);
            end%if
            if this.generated.VertexLabels
                for ii = 1:this.sides
                    delete(this.sketches.VertexLabels(ii))
                end%ii
            end%if
            if this.generated.SegmentLabels
                for ii = 1:this.sides
                    delete(this.sketches.SegmentLabels(ii))
                end%ii
            end%if
        end%function
        
        %Memory allocation
        function Calloc(this,sides)
            this.XY = zeros(sides,2); %Buffer for XY coordinates.
            this.nxy = zeros(sides,2); %Buffer for 2D components of normals.
            %kappa = zeros(sides,1); %Buffer for discrete curvature measure.
            %L = zeros(sides,1);
            %angles = zeros(sides,1)
        end%function
        function ReCalloc(this,new_sides)
            %This routine changes the size of various storage buffers. It
            %is meant to be used as part of the "Redefine" routines.
            if new_sides == this.sides
                return;
            elseif new_sides < this.sides
                last = 1 + new_sides;
                %Eliminate excess storage.
                this.XY(last:this.sides,:) = [];
                this.nxy(last:this.sides,:) = [];                
                %Do the same for the curves.
                if this.generated.Curve 
                    this.sketches.Curve.XData(last:end) = [];
                    this.sketches.Curve.YData(last:end) = [];
                end%if
                if this.generated.Normals 
                    this.sketches.Normals.XData(last:end) = [];
                    this.sketches.Normals.YData(last:end) = [];
                    this.sketches.Normals.UData(last:end) = [];
                    this.sketches.Normals.VData(last:end) = [];
                end%if                     
            elseif new_sides > this.sides
                %Increase storage.
                delta = this.sides - new_sides;
                this.XY = [this.XY; zeros(delta,2)];
                this.nxy = [this.XY; zeros(delta,2)];
                %TerminateGraphics(this);
                %InitializeGraphics(this);
                if this.generated.Curve 
                    this.sketches.Curve.XData = [this.sketches.Curve.XData,zeros(1,delta)];
                    this.sketches.Curve.YData = [this.sketches.Curve.YData,zeros(1,delta)];
                    for ii = 1:this.sides
                        this.sketches.Curve.XData(ii) = 0;
                        this.sketches.Curve.YData(ii) = 0;
                    end%ii
                end%if
                if this.generated.Normals 
                    this.sketches.Normals.XData = [this.sketches.Normals.XData,zeros(1,delta)];
                    this.sketches.Normals.YData = [this.sketches.Normals.YData,zeros(1,delta)];
                    this.sketches.Normals.UData = [this.sketches.Normals.UData,zeros(1,delta)];
                    this.sketches.Normals.VData = [this.sketches.Normals.VData,zeros(1,delta)];
                end%if      
            end %if
            this.sides = new_sides;
        end%function
        
        %Redefinition routines (modification)
        function RedefineFromList(this,points,XY,bool)
            if nargin < 3
                bool = false;
            end%if
            this.open = bool;
            %A repeated point will override bool.
            if XY(1,1) == XY(points,1) && XY(1,2) == XY(points,2)
                points = points - 1;
                this.open = false;
            end%if
            if points < 3
                error('Need at least 3 sides for a polygon.')
            end%if
            %If redefining reset the validity flags.
            this.ReCalloc(points);
            this.sides
            for ii = 1:this.sides
                this.XY(ii,1) = XY(ii,1);
                this.XY(ii,2) = XY(ii,2);
            end%ii
            Measure(this);           
            this.valid = true;
            Refresh(this);
        end%function
        function RedefineAsWeightedAverage(this,poly1,w1,poly2,w2)
            if poly1.sides ~= poly2.sides
                error('Input polygons do not have equal sides')
            end%if
            if poly1.open ~= poly2.open
                poly.open = false;
                warning('Input polygons have different "open" flags. Defaulting to closed.');
            else
                this.open = poly1.open;
            end%if
            this.ReCalloc(poly1.sides);
            w = w1 + w2;%Total weight.
            for ii = 1:poly1.sides
                this.XY(ii,1) = (poly1.XY(ii,1)*w1 + poly2.XY(ii,1)*w2)/w;
                this.XY(ii,2) = (poly1.XY(ii,2)*w1 + poly2.XY(ii,2)*w2)/w;
            end%ii
            this.Measure;
            this.Refresh;
            this.valid = true;
        end%function
        function RedefineAsRegular(this)
            length = this.perimeter/this.sides;
            
        end%function

        %DEPRECATE AFTER tensor.m matures
        %Affine Transformations (modification) MOVE to tensor.m
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
            Refresh(this);
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
            Refresh(this)
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
            this.Measure;
            Refresh(this);
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
            Refresh(this);
        end%function
        %DEPRECATE AFTER tensor.m matures

        
        %Custom Transformations (modification)
        %DEPRECATE ALL ONCE tensor.m matures.
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
                if ~(exist('polynomial.m','file') == 2)
                    error('Smoothing of open curves requires the file "polynomial.m"');
                end%if
                this.SmoothOpen(lambda,N);
            end%if
            
            this.Measure;            
            this.Refresh;            
        end%function
        function Disperse(this,lambda,N)
        end%function
        function SnipeSelfIntersections(this,sniper)
            %If the polygon has any self intersections, this routine will
            %apply blending techniques to correct such a flaw. A blending
            %curve is produced and discretely evaluated as many times as
            %there are points in the faulty region.
            if this.simple
                warning('This polygon was flagged as being "simple." Are you sure it has Self-Intersections?')
                return;
            end%if
            
            %Select your torpedo.
            if strcmpi(sniper,'ellipse') == 1
                %Fit a minimum eccentricity fillet through the affected
                %region.
                for ii = 1:this.SI
                    
                end%ii
            elseif strcmpi(sniper,'circle') == 1
                %Fit a circular fillet through the affected region
                
                
                %WIP

            elseif strcmpi(sniper,'bezier') == 1
                %Fit a quadratic Bezier curve through the affected region.
                
                %WIP
            elseif strcmpi(sniper,'chamfer') == 1
                %Connects a line between the point immediately before the
                %flaw and the point immediately after the flaw. An evenly
                %spaced number of points is evaluated in between this line.
                
                %WIP
            end%if
            
            %Update flags.
            this.refresh.SelfIntersect = false;
            this.simple = true;
        end%function
        
        %Move to low-level
        % replace "this" with "sides", and "XY".
        function SmoothClose(this,lambda,N)
            %Append a weighted "pass" of finite difference scheme onto a
            %sequence of 2D coordinates. The finite difference scheme is
            %for some discrete measure of the sequence's derivatives.
            %Even-numbered derivatives result in a diffusion effect.
            %Odd-numbered derivatives result in a dispersion effect. This
            %routine assumes that the sequence represents a closed polygon.
            
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
                %Hardcoded Central finite difference schema.
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
                       sequence(1) = -(N-1)*0.5;
                       for ii = 2:N
                           sequence(ii) = sequence(ii-1) + 1;
                       end%ii
                       FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
                       clear sequence;%free
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
           
            
            clear FD;%free
            clear ahead_x; %free
            clear ahead_y; %free
            clear trail_x; %free
            clear trail_y; %free
            
            if this.AABB_present
                this.aabb.RedefineFromList(this.sides,this.XY);
            end%if
        end%
        function SmoothOpen(this,lambda,N)
            %This routine is meant to be identical to its "closed"
            %counterpart except that it will treat the sequence as if it
            %represents an open curve. The finite difference scheme is
            %first "forward", and then switches to "central once enough
            %stencil points are available on both sides for central
            %difference to occur. Once the other end is reached, the scheme
            %switches to a backward difference.
            skip_endpoints = true;
            
            
            %If # of stencil points is not specified, default to 3.
            if nargin < 3
                N = 3;
            end%if
            
            %Like with open curve smoothing, a minimal copy of "trailing"
            %values is needed. The "ahead" values are not needed though.
            Nm1 = N - 1;
            Nm1o2 = Nm1*0.5; %Nm1o2 = "N minus 1 over 2"
            Np1o2 = Nm1o2 + 1; %This index corresponds to the middle of the central difference scheme.
            trail_x  = zeros(1,Nm1o2); %Stores the trailing (N-1)...
            trail_y  = zeros(1,Nm1o2); %... stencil points.
            
            

            for ii = 1:Nm1o2 %This time, more points need to be remembered...
                trail_x(ii) = this.XY(ii,1); %...this is because towards...
                trail_y(ii) = this.XY(ii,2); %...the end backward schemes...
            end%ii                            ...are applied.
            
            trail2_x = zeros(1,N);
            trail2_y = zeros(1,N);
            for ii = 1:N
                trail2_x(ii) = this.XY(this.sides-N+ii,1);
                trail2_y(ii) = this.XY(this.sides-N+ii,2); 
            end%ii
            
            %Derive the first forward difference scheme.
            sequence = zeros(1,N);%Initialize a sequence buffer.
            for ii = 1:N
                sequence(ii) = ii - 1; %Set to [0, +1, +2, +3,....,+N]
            end%ii
            FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
            
            %This first loop will smooth the first "(N-1)/2" points. The
            %first one is smoothed according to a forward difference (st by
            %"FD". The next "(N-1)/2 - 1" points are smoothed according to
            %shifted forward differences.
            for ii = 1:Nm1o2
                x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
                y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
                
                %Reference the points behind (stored in trail_x/trail_y
                %buffers).
                for jj = 1:(ii - 1) 
                    x_diffusion = x_diffusion + trail_x(jj)*FD(jj);
                    y_diffusion = y_diffusion + trail_y(jj)*FD(jj);
                end%jj
                %Reference the points ahead
                for jj = ii:N
                    x_diffusion = x_diffusion + this.XY(jj,1)*FD(jj);
                    y_diffusion = y_diffusion + this.XY(jj,2)*FD(jj);
                end%%jj
       
                %*No need to update the trail_x and trail_y buffers here.*
                
                %Spend the original points.
                if ~skip_endpoints
                    this.XY(ii,1) = this.XY(ii,1) + lambda*x_diffusion;
                    this.XY(ii,2) = this.XY(ii,2) + lambda*y_diffusion;
                end
                %Update Lagrange interpolation sequence.
                for jj = 1:N
                    sequence(jj) = sequence(jj) - 1;
                end%jj
                %Derive the new forward difference scheme.
                FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
            end%ii
            
            
            %By the end of the first loop, "FD" has the coefficients for a
            %central difference. We re-introduce the "spanner" and "sm1"
            %variables from the open smoothing version to reduce the
            %overhead of copying.
            spanner = 1; 
            sm1 = spanner - 1;
            for ii = Np1o2:(this.sides - Nm1o2)
                x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
                y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
                
                %Contribution of the current point.
                x_diffusion = x_diffusion + this.XY(ii,1)*FD(Np1o2);
                y_diffusion = y_diffusion + this.XY(ii,2)*FD(Np1o2);
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
            
            %UNCOMMENT WITH "%{
            
            
            %At this stage, all points elligible for central difference
            %have been spent. Backward difference on the last "(N-1)/2"
            %points now ensus.   
            for ii = (this.sides - Nm1o2):this.sides
                %Update Lagrange sequence to derive backward difference
                %schema.
                for jj = 1:N
                    sequence(jj) = sequence(jj) - 1; %Set to [-N, 1-N, 2-N, 3-N,....,0]
                end%ii
                FD = polynomial.LagrangeFiniteDifference(N,sequence,2);
                
                x_diffusion = 0;%Set/Reset the X-coordinate smoothing.
                y_diffusion = 0;%Set/Reset the Y-coordinate smoothing.
                %Reference the points behind (stored in trail_x/trail_y
                %buffers).
                for jj = 1:(ii - (this.sides - Nm1o2))
                %for jj = 1:(Np1o2 + ii -this.sides + Nm1o2)
                    %t_idx = jj + sm1; %Index into the trailing buffer.
                    %t_idx = t_idx - Nm1*(Nm1 - t_idx < 0); %Push it back if it overflows.                    
                    x_diffusion = x_diffusion + trail2_x(jj)*FD(jj);
                    y_diffusion = y_diffusion + trail2_y(jj)*FD(jj);
                end%jj

                %Reference the points ahead
                %for jj = ((Np1o2 + ii -this.sides + Nm1o2+1)):N
                for jj =(1 + ii - (this.sides - Nm1o2)):N
                    idx = ii + sequence(jj) + 1;
                    x_diffusion = x_diffusion + this.XY(idx,1)*FD(jj);
                    y_diffusion = y_diffusion + this.XY(idx,2)*FD(jj);
                end%%jj                                
                
                
                %Update the permuted buffer before expending the original
                %point.
                %trail_x(spanner) = this.XY(ii,1);
                %trail_y(spanner) = this.XY(ii,2);
                spanner = (spanner + 1)*(spanner < Nm1) + 1*(spanner >= Nm1);
                %spanner = (spanner + 1)*(spanner < Nm1o2) + 1*(spanner >= Nm1o2);
                sm1 = spanner - 1;

                %Apply the smoothing operation InPlace (spend the original
                %point).
                if ~skip_endpoints
                this.XY(ii,1) = this.XY(ii,1) + lambda*x_diffusion;
                this.XY(ii,2) = this.XY(ii,2) + lambda*y_diffusion;     
                end
            end%ii
            %}
            
            clear FD;%free
            clear sequence; %free
            clear trail_x; %free
            clear trail_y; %free

            if this.AABB_present
                this.aabb.RedefineFromList(this.sides,this.XY);
            end%if
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

        %NEEDS TO BE COMPLETED
        function tight = Rollercoaster(this,Rf,sink)
                        
            sink = 1;
            
            %This routine applies a "rolling" circle type of operation to
            %flag turns that exceed the radius of the tighest turn that the
            %polygon should make.
            
            tight = ones(1,this.sides); %Allocate memory for an array of flags.
            
            %For this to work, the polygon must be oriented inwards.

            %mult*
            
            

            %Sort either X or Y coordinates.
            [X_sorted,X_perm] = sort(this.XY(:,1));
            %Need to first find the where the first point ranks in the
            %sorted list.
            needle = 1;
            while X_perm(needle) ~= 1
                needle = needle + 1;
                if needle == length(X_perm)
                    %needle = 1;
                    break;
                end%if
            end%while
            old_needle = needle;

            flags = scatter3(...
                this.canvas,...
                this.XY(:,1),...
                this.XY(:,2),...
                NaN(this.sides,1),...
                'Marker','s',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,0]);
                        
            %Definition of the "Rolling Circle."                      
            F = circle.CreateXYR(...
                0.5*(this.XY(1,1) + this.XY(2,1)) + Rf*this.nxy(1,1),... %X-center.
                0.5*(this.XY(1,2) + this.XY(2,2)) + Rf*this.nxy(1,2),... %Y-center
                abs(Rf));%Radius
            F.SetCanvas(this.canvas);
            F.Show;
            F.sketches.Curve.Color = [1,1,0];
            F.sketches.Orientation.Visible = 'on';
            F.sketches.Orientation.Color = [1,1,0];
            for kk = 1:3
                switch kk
                    case 1
                        start = 1;
                        step1 = 1;
                        step2 = 2;
                        finish = this.sides - 2;
                    case 2
                        if this.open
                            break;
                        else
                            start = this.sides - 1;
                            step1 = 1;
                            step2 = 2 - this.sides;
                            finish = this.sides - 1;
                        end%if
                    case 3
                        if this.open
                            break;
                        else
                            start = this.sides;
                            step1 = 1 - this.sides;
                            step2 = 2 - this.sides;
                            finish = this.sides;
                        end%if
                end%switch
                for ii = start:step1:finish
                    iip1 = ii + step1;
                    iip2 = ii + step2;
                    xc2 = 0.5*(this.XY(iip1,1) + this.XY(iip2,1)) + sink*Rf*this.nxy(iip1,1);
                    yc2 = 0.5*(this.XY(iip1,2) + this.XY(iip2,2)) + sink*Rf*this.nxy(iip1,2);
                    
                    %Compute Arc length traversed. (ANIMATION + FUNCTIONALITY)
                    
                    dsx = this.XY(iip1,1) - this.XY(ii,1);
                    dsy = this.XY(iip1,2) - this.XY(ii,2);
                    ds = sqrt(dsx*dsx + dsy*dsy);
                    F.RotateByArc(ds);
                    %}
                    
                    %Compute distance moved by the center (ANIMATION + FUNCTIONALITY)
                    dxc = xc2 - F.xc;
                    dyc = yc2 - F.yc;
                    F.Displace(dxc,dyc);

                    %needle = old_needle;
                    if dxc < 0
                        %Locate the needle of the next point. Look in the negative sort
                        %direction as the change in coordinate is negative.
                        while X_perm(needle) ~= ii
                            %{
                            needle = needle - 1;
                            if needle == 0 %Keep index within bounds.
                                needle = this.sides;
                            end%if
                            %}
                            needle = (needle - 1)*(needle > 1) + this.sides*(needle == 1);
                        end%while
                        old_needle = needle;
                    elseif dxc > 0
                        %Locate needle of the next point. Since the change is positive,
                        %look forwards.
                        while X_perm(needle) ~= ii
                            %{
                            needle = needle + 1;
                            if needle == this.sides + 1 %Keep index within bounds.
                                needle = 1;
                            end%if
                            %}
                            needle = (needle + 1)*(needle < this.sides) + 1*(needle == this.sides);
                        end%while
                        old_needle = needle;
                    else %if dxc = 0
                        fprintf('Cannot handle vertical lines yet! (Murphy''s Law)\n');
                    end%if
                    
                    %Needle found, now, check all the points in the negative sort
                    %direction that are within xc - R.
                    nL = needle; %Need a counter to subtract from the first needle.
                    while X_sorted(nL) > (F.xc - F.R) %&& X_sorted(1) > (F.xc - F.R)
                        if F.ContainsPoint(this.XY(X_perm(nL),1),this.XY(X_perm(nL),2))
                            flags.ZData(X_perm(nL)) = 0;%"Un-NaN" the marker.
                            tight(X_perm(nL)) = -1;
                        end%if
                        nL = nL - 1;
                        if nL == 0 %If even the minimum X-value is to the left, break.
                            %needleL = poly.sides;
                            break;
                        end%if
                    end%while
                    
                    %Cheack nearby "right" points.
                    nR = needle;
                    while X_sorted(nR) < (F.xc + F.R) %&& X_sorted(poly.sides) < (F.xc + F.R)
                        if F.ContainsPoint(this.XY(X_perm(nR),1),this.XY(X_perm(nR),2))
                            flags.ZData(X_perm(nR)) = 0;%"Un-NaN" the marker.
                            tight(X_perm(nR)) = -1;
                        end%if
                        nR = nR + 1;
                        if nR == (this.sides + 1)
                            %needleR = 1;
                            break;
                        end%if
                    end%while
                    needle = old_needle;
                    %input('Next?');
                    pause(0.00001);
                end%ii
            end%kk
            delete(F); %Deallocate the circle

            %Filter any ghost points
            for kk = 1:3
                switch kk
                    case 1
                        start = 2;
                        back = 1; %back
                        next = 1;
                        finish = this.sides - 1;
                    case 2
                        start = this.sides;
                        back = 1;
                        next = 1 - this.sides;
                        finish = 2;
                    case 3
                        start = 1;
                        next = 1;
                        back = 1 - this.sides;
                        finish = 1;
                end%switch
                for ii = start:next:finish
                    if tight(ii - back) == -1 && tight(ii + next) == -1
                        tight(ii) = -1;
                        flags.ZData(ii) = 0;
                    end%if
                end%ii
                if this.open
                    break;
                end%if
            end%for
        end%function
        function RedefineAsComplexStar(this)
            if mod(this.sides,2) == 0
                warning('Support for Complex star generation is limited to odd-sides polygons.');
                return;
            end%if
            
            length = this.perimeter/this.sides;
        end%function
        %NEEDS TO BE COMPLETED

        %NEEDS TO BE WRITTEN
        function RedefineAsSimpleStar(this)
            
        end%function
        function targets = SpotHighCurvature
        end%function
        function targets = SpotSelfIntersections
        end%function
        function authenticate(this)
            %Authenticator: (Corrects as much as it can before having to
            %error).
            this.valid = isempty(this.xc) + isempty(this.yc) + isempty(this.XY);
            this.valid = ~this.valid;
        end%function
        function inside = PointInPolygon(this,XY_in)
            %This function treats the polyline as if it were a polygon
            
            %First, determine whether the point is even inside the
            %polygon's AABB
        end%function
        %NEEDS TO BE WRITTEN

        
        function DetectSharpTurns(this,Rf)

            bad_points = this.Rollercoaster(Rf,1);
            %Extract bands of continuous bad flags.
            [L,idx] = polygon.DetectSignChanges(this.sides,bad_points);
            idx
            
            %Uncomment to enable blending
            
            %If the number of sign changes is even, then the the beginning
            %points do not form a tight turn. If the number of sign changes
            %is odd, it is possi
            targets = floor(0.5*L) + (bad_points(1) ~= bad_points(this.sides));
            bands = 0;
            kk = 1;
            %while bands < targets
            while kk < L
                if bad_points(idx(kk)) == +1
                    %Skip points with a "+1."
                    kk = kk + 1;
                else %A point with "-1" is found.
                    
                    %{
                    %DOT PRODUCT EXCLUDER:
                    dx1 = this.XY(idx(kk),1) - this.XY(idx(kk)-1,1);
                    dy1 = this.XY(idx(kk),2) - this.XY(idx(kk)-1,2);
                    dx2 = this.XY(idx(kk+1)+1,1) - this.XY(idx(kk+1)-1,1);
                    dy2 = this.XY(idx(kk+1)+1,2) - this.XY(idx(kk+1)-1,2);
                    if dx1*dx2 + dy1*dy2 > 0
                        kk = kk  + 2;
                        continue;
                    end%if
                    %}
                    
                    N = idx(kk + 1) - idx(kk);
                    %{
                    [this.XY,status,circ] = polygon.CircularFilletBlend(...
                        this.sides,...
                        this.XY,...
                        N,... %Number of points to put on Filet.
                        (idx(kk):(idx(kk+1)-1)),... %Indices to fillet.
                        abs(Rf),... %Fillet Radius.
                        this.XY(idx(kk)-1,1)    ,this.XY(idx(kk)-1,2),...   %x1,y1
                        this.XY(idx(kk),1)      ,this.XY(idx(kk),2),...     %x2,y2
                        this.XY(idx(kk + 1),1)  ,this.XY(idx(kk + 1),2),... %x3,y3
                        this.XY(idx(kk + 1)+1,1),this.XY(idx(kk + 1)+1,2)); %x4,y4
                    %}
                    status = false;
                    if ~status
                        fprintf('Circular fillet failed! Switching to LEE blend.\n');
                        
                        idx_m = round(0.5*(idx(kk) + (idx(kk+1)-1))); %index of some middle point.
                        
                        idx1 = idx(kk)-1;
                        idx2 = idx(kk);
                        idx3 = idx(kk + 1);
                        idx4 = idx(kk + 1)+1;
                        
                        if idx4 > this.sides
                            break;
                        end%if
                        
                        [this.XY,status,elli] = polygon.EllipticalFilletBlend(...
                            this.sides,...
                            this.XY,...
                            N,(idx(kk):(idx(kk+1)-1)),...
                            1,...
                            this.XY(idx1,1),this.XY(idx1,2),... %x1,y1
                            this.XY(idx2,1),this.XY(idx2,2),... %x2,y2
                            this.XY(idx3,1),this.XY(idx3,2),... %x3,y3
                            this.XY(idx4,1),this.XY(idx4,2),... %x4,y4
                            idx_m); 
                        elli.SetCanvas(this.canvas);
                        elli.SetColor(this.color);
                        elli.Toggle('Curve','MajorRadius','MinorRadius');
                        elli.sketches.Curve.LineStyle = '--';
                        
                        %{
                        elli.Draw(this.canvas,100);
                        elli.DrawMajorRadius(this.canvas);
                        elli.DrawMinorRadius(this.canvas)
                        elli.sketch.Color = this.color;
                        elli.sketch.LineStyle = '--';
                        %}
                    else
                        circ.SetCanvas(this.canvas);
                        circ.SetColor([1,1,1]);
                        circ.Show;
                    end                    
                    kk = kk + 2;
                    bands = bands + 1;
                end%if
            end%while
            
            %If the first sign change is from "-1" to "+1" there is a
            %possibility that a tight turn encompasses beginning and ending
            %indices.
            if ~this.open && L > 0
                if bad_points(idx(1)) == +1 && bad_points(idx(L)) == -1
                %If the curve is open, ignore it. Else, look into it.
                N = this.sides - idx(L) + idx(1);
               % bll = [idx(L):this.sides,1:(idx(1)-1)];
                idx_m = 1;
                [this.XY,status,elli] = polygon.EllipticalFilletBlend(...
                            this.sides,...
                            this.XY,...
                            N,[idx(L):this.sides,1:(idx(1)-1)],...
                            1,...
                            this.XY(idx(L)-1,1),this.XY(idx(L)-1,2),...   %x1,y1
                            this.XY(idx(L),1)  ,this.XY(idx(L),2),...     %x2,y2
                            this.XY(idx(1),1),this.XY(idx(1),2),... %x3,y3
                            this.XY(idx(1)+1,1)  ,this.XY(idx(1)+1,2),...
                            idx_m); %x4,y4
                        
                        
                        elli.SetCanvas(this.canvas);
                        elli.SetColor(this.color);
                        elli.Toggle('Curve','MajorRadius');
                        elli.sketches.Curve.LineStyle = '--';
                        
                        
                        %{
                        elli.Draw(this.canvas,100);
                        elli.DrawMajorRadius(this.canvas);
                        elli.DrawMinorRadius(this.canvas)
                        elli.sketch.Color = this.color;
                        elli.sketch.LineStyle = '--';
                    %}
                end%if
            end%if
            %}
            
            
            this.Measure;
            %%%%%%%%%%%%%SECOND ATTEMPT

            
  
            
            %delete(flags); %Deallocate the flags.
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
            Refresh(this);
        end%function
            
        %Orientation
        %DEPRECATE
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
            this.normals_computed = true;
        end%function
        %DEPRECATE
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
            RefreshNormals(this);
        end%function
        
        %Measure metric properties about this polyline.
        function Measure(this)
            %Quantifyable metrics.
            this.xc = 0; %Running sum for X-coordinate of the centroid.
            this.yc =0; %Running sum for Y-coordinate of the centroid.
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
            
            %old_dx = this.XY(1,1) - this.XY(this.sides,1);
            %old_dy = this.XY(1,2) - this.XY(this.sides,2);
            for kk = 1:2
                switch kk 
                    case 1
                        start = 1;
                        step = 1;
                        finish = this.sides - 1;
                    case 2
                        start = this.sides;
                        step = 1 - this.sides;
                        finish = 2;
                end%switch
                for ii = start:step:finish
                    iip1 = ii +step;
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
            end%kk
                        
            %Area and Centroid: Correction factors.
            this.area = this.area*0.5;
            this.xc = this.xc/(6*this.area);
            this.yc = this.yc/(6*this.area);
                        
            this.orientation = (this.area >= 0) - (this.area < 0);%If signed area is positive, polygon is oriented inwards.
            this.area = this.area*((this.area > 0) - (this.area < 0));%Branchless way of taking absolute value.
            
            %Convexity inference: If the number of sign changes in the
            %X-coordinate plus the number of sign changes in the
            %Y-coordinate is exactly 2 then the polygon is convex (assuming
            %that the first point of the polygon has the minimum X and Y
            %values for coordinates.
            this.convex = (sign_xchange + sign_ychange - (min_x ~= this.XY(1,1)) - (min_y ~= this.XY(1,2)) == 2);        
            ComputeNormals(this); %DEPRECATE
            
            if this.AABB_present
                this.aabb.RedefineFromList(this.sides,this.XY)
            end%if 
            
            %Check whether the polygon self-intersects.
            if this.sides == 3
                this.simple = true;
            else
                [this.xSI,this.ySI,this.sSI] = selfintersect(this.XY(:,1),this.XY(:,2));
                this.simple = isempty(this.xSI);
                %If SI's are detected, it is probably best to track them.
                if ~this.simple
                    this.SI = length(this.xSI);
                    this.refresh.SelfIntersect = true;
                end%if
            end%if
            
        end%function
        function MeasureShoelace(this)
            %Given a valid definition of a polygon this routine measures
            %the "Shoelace" quantities: Signed Area, Signed X-Centroid, and
            %Signed Y-Centroid.
            this.xc = 0; %Running sum for X-coordinate of the centroid.
            this.yc =0; %Running sum for Y-coordinate of the centroid.
            this.area = 0; %Initialize a running sum for the area
            
            %Allocate memory for the local metrics
            this.metL.areas = zeros(1,this.sides);
            this.metL_alloc.areas = true;
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        next = 1;
                        finish = this.sides;
                    case 2
                        start = this.sides;
                        next = 1 - this.sides;
                        finish = 2;
                end%
                for ii = start:next:finish
                    iip1 = ii + next;
                    %Make this a function
                    dA = this.XY(ii,1)*this.XY(iip1,2) - this.XY(iip1,1)*this.XY(ii,2);
                    %Make this a function
                    this.area = this.area + dA;
                    this.xc = this.xc + dA*(this.XY(ii,1) + this.XY(iip1,1));
                    this.yc = this.yc + dA*(this.XY(ii,2) + this.XY(iip1,2));
                end%ii
            end%kk
            %Apply correction factors to the shoelace formula.
            this.area = this.area*0.5;
            this.xc = this.xc/(6*this.area);
            this.yc = this.yc/(6*this.area);
        end%function
        function MeasureLAKN(this)
            %"LAKN" = Lengths, Angles, Kappa, and Normals.
            %Given a valid definition of a polygon, this routine measures
            %the lengths (L-2 norms) of the edges of the polygon's
            %boundaries and the interior angles of the vertices.
            
            %Allocate memory for the local metrics
            this.metL.lengths = zeros(1,this.sides);
            this.metL.angles = zeros(1,this.sides);
            this.metL.kappa = zeros(1,this.sides);
            this.metL.nxy = zeros(this.sides,2);

            %Update flag buffers.
            this.metL_alloc.lengths = true;
            this.metL_alloc.angles = true;
            this.metL_alloc.kappa = true;
            this.metL_alloc.nxy = true;
            
            %Manually compute the lengths and displacements of the last
            %segment. This avoids recomputing the "dx" and "dy" quantities
            %twice except for the last segment.
            dxiim1 = this.XY(1,1) - this.XY(this.sides,1);
            dyiim1 = this.XY(1,2) - this.XY(this.sides,2);
            for kk = 1:3
                switch kk
                    case 1 %First Vertex
                        start = 1;
                        next = 1;
                        back = 1 - this.sides;
                        finish = 1;
                    case 2 %Middle Vertices
                        start = 2;
                        next = 1;
                        back = 1;
                        finish = this.sides - 1;
                    case 3 %Last Vertex
                        start = this.sides;
                        next = 1 - this.sides;
                        back = 1;
                        finish = 3;
                end%switch
                for ii = start:next:finish
                    iip1 = ii + next;
                    iim1 = ii - back;
                    dxii = this.XY(iip1,1) - this.XY(ii,1);
                    dyii = this.XY(iip1,2) - this.XY(ii,2);
                    
                    %Compute the length;
                    this.metL.lengths(ii) = sqrt(dxii*dxii + dyii*dyii);
                    
                    %Compute the unit normal components using a hardcoded
                    %90 degree CCW rotation.
                    this.metL.nxy(ii,1) = -dyii/this.metL.lengths(ii);
                    this.metL.nxy(ii,2) = +dxii/this.metL.lengths(ii);

                    %Computation of internal angles and curvature.
                    this.metL.angles(ii) = acos((dxii*dxiim1 + dyii*dyiim1)/(this.metL.lengths(ii)*this.metL.lengths(iim1)));                    
                    if dxiim1*dyii - dxii*dyiim1 < 0 %Determine reflex vs. convex vertex from cross product.
                        this.metL.angles(ii) = 2*pi - this.metL.angles(ii);
                    end%if               
                    
                    %Computation of curvature.
                    this.metL.kappa(ii) = 2*(pi - this.metL.angles(ii))/(this.metL.lengths(ii) + this.metL.lengths(iim1));
                                        
                    %Remember the quantity computed here to avoid
                    %recomputation.
                    dxiim1 = dxii;
                    dyiim1 = dyii;
                end%ii
            end%for
        end%function
        function MeasureNormals
        end%function
        
        
        function FilterSelfIntersections(this)
            % THIS NEEDS A LOT MORE WORK!!!
            %This routine is called to help handle self-intersection
            %problems, namely "loop-within-loop" phenomenona. This is
            %because the correction algorithms need two blending edges.
            if this.simple
                return;
            end%if
            deltas = zeros(this.SI,1);
            for ii = 1:this.SI
                deltas(ii) = this.sSI(ii,2) - this.sSI(ii,1);
            end%ii
            [deltas,perm] = sort(deltas);
            for ii = this.SI
                
            end%ii
            
            
            new_SI = zeros(this.SI,2);
            keep = 1;
            kk = 0;
            while kk < this.SI
                kk = kk + 1;
                ww = kk + 1;
                new_SI(keep,1) = this.sSI(kk,1);
                new_SI(keep,2) = this.sSI(kk,2);
                                
                while this.sSI(kk,2) > this.sSI(ww,1) && this.sSI(ww,1) > this.sSI(kk,1) && this.sSI(kk,2) > this.sSI(ww,2) && this.sSI(ww,2) > this.sSI(kk,1)
                    ww = ww + 1;
                end%while
                keep = keep + 1;
                kk = ww;
            end%while
       
        end%function
        

        function XY_out = Offset(this,offset)
            if ~this.valid
                error('Polygon is flagged as invalid.');
            end%if
            if this.open
                XY_out = polygon.OffsetOpen(this.sides,this.XY,this.nxy,offset);
            else
                XY_out = polygon.OffsetClosed(this.sides,this.XY,this.nxy,offset);
            end%if
            
        end%function
        

       
        
    end%methods
    %Graphical setups.
    methods
        
        %Whether the polygon is to be drawn.       
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
        function SetRefreshRate(this,rate)
            %Set a target refresh rate for the Refresh routine.
            if rate == 0
                rate = 1;
            end%if
            this.refresh_rate = rate*(-1*(rate < 0));
        end%function
        function SetName(this,string)
            this.name = string;
            if this.generated.Curve
                this.sketches.Curve.DisplayName = string;
            end%if
            if this.generated.Normals
                this.sketches.Normals.DisplayName = [string,'''s Normals'];
            end%if
            if this.generated.Centroid
                this.sketches.Centroid.DisplayName = [string,'''s Centroid'];
            end%if
            
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.SetName(string)
            end%if
        end%function        
        function SetColor(this,RGB)
            this.color = RGB;
            if this.generated.Curve
                this.sketches.Curve.Color = RGB;
                this.sketches.LastEdge.Color = RGB;
            end%if
            if this.generated.Normals
                this.sketches.Normals.Color = RGB;
            end%if
            if this.generated.Centroid
                this.sketches.Centroid.MarkerFaceColor = RGB;
            end%if
            if this.generated.VertexLabels
                for ii = 1:this.sides
                    this.sketches.VertexLabels(ii).Color = RGB;
                end%ii
            end%if
            if this.generated.SegmentLabels
                for ii = 1:this.sides
                    this.sketches.SegmentLabels(ii).Color = RGB;
                end%if
            end %if

            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.SetColor(RGB);
            end%if
            
        end%function
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.generated.Curve
                this.sketches.Curve.Parent = ax;
                this.sketches.LastEdge.Parent = ax;
            end%if
            if this.generated.Centroid
                this.sketches.Centroid.Parent = ax;
            end%if
            if this.generated.Normals
                this.sketches.Normals.Parent = ax;
            end%if
            if this.generated.VertexLabels
                for ii = 1:this.sides
                    this.sketches.VertexLabels(ii).Parent = ax;
                end%ii
            end%if
            if this.generated.SegmentLabels
                for ii = 1:this.sides
                    this.sketches.SegmentLabels(ii).Parent = ax;
                end%ii
            end%if
            
            if this.AABB_present
          %      this.aabb.SetCanvas(ax)
            end%if
        end%function
        function SetLineWidth(this,width)
            this.linewidth = width;
            if this.generated.Curve
                this.sketches.Curve.LineWidth = width;
                this.sketches.LastEdge.LineWidth = width;
            end%if
        end%function
        
        %Create the Graphical objects.
        %DEPRECATE
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            this.GenerateCurve;
            this.RefreshCurve;
            
            this.GenerateCentroid;
            this.RefreshCentroid;

            this.GenerateNormals;
            this.RefreshNormals;
            
            this.GenerateSelfIntersections;
            this.RefreshSelfIntersect;

            
            %Bounding Box initialization.
            if this.AABB_present
                this.aabb.InitializeGraphics;
                this.aabb.ToggleBox;
            end%if
            
            this.graphics_initialized  = true;
        end%function.           
        function TerminateGraphics(this)
            %Release system resources used to render graphics.
            if ~isempty(this.sketches.Curve)
                delete(this.sketches.Curve);
                this.generated.Curve = false;
            end%if
            if ~isempty(this.sketches.Normals)
                delete(this.sketches.Normals);
                this.generated.Normals = false;
            end%if
            if ~isempty(this.sketches.Centroid)
                delete(this.sketches.Centroid);
                this.generated.Centroid = false;
            end%if            

            if ~isempty(this.sketches.VertexLabels)
                for ii = this.sides:-1:1
                    delete(this.sketches.VertexLabels(ii));
                end%ii
                this.generated.VertexLabels = false;
            end%if
            if ~isempty(this.sketches.SegmentLabels)
                for ii = this.sides:-1:1
                    delete(this.sketches.SegmentLabels(ii));
                end%ii
                this.generated.SegmentLabels = false;
            end%if
            %Update graphics flag.
            this.graphics_initialized = false;
        end%function
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
                        this.GenerateCurve;
                    end%if
                    if ~this.updated.Curve 
                        this.RefreshCurve;
                    end%if
                    this.ToggleVisible(this.sketches.Curve);
                    this.ToggleVisible(this.sketches.LastEdge);
                elseif strcmpi(varargin{ii},'LastEdge') == 1
                    this.ToggleVisible(this.sketches.LastEdge);

                elseif strcmpi(varargin{ii},'Centroid') == 1
                    if ~this.generated.Centroid
                        this.GenerateCentroid;
                    end%if
                    if ~this.updated.Centroid
                        this.RefreshCentroid;
                    end%if
                    this.ToggleVisible(this.sketches.Centroid);
                elseif strcmpi(varargin{ii},'Normals') == 1
                    if ~this.generated.Normals
                        this.GenerateNormals;
                    end%if
                    if ~this.updated.Normals
                        this.RefreshNormals;
                    end%if
                    this.ToggleVisible(this.sketches.Normals);
                elseif strcmpi(varargin{ii},'VertexLabels') == 1
                    if ~this.generated.VertexLabels
                        this.GenerateVertexLabels;
                    end%if
                    if ~this.updated.VertexLabels
                        this.RefreshVertexLabels;
                    end%if
                    for jj = 1:this.sides
                        this.ToggleVisible(this.sketches.VertexLabels(jj));
                    end%ii
                elseif strcmpi(varargin{ii},'SegmentLabels') == 1
                    if ~this.generated.SegmentLabels
                        this.GenerateSegmentLabels;
                    end%if
                    if ~this.updated.SegmentLabels
                        this.RefreshSegmentLabels;
                    end%if
                    this.ToggleVisible(this.sketches.SegmentLabels);    
                elseif strcmpi(varargin{ii},'AABB') == 1
                    this.aabb.Show;
                else
                    warning('Unrecognizable graphics option. Valid options are: "Curve", "Normals","Centroid", "VertexLabels", "SegmentLabels", and "AABB"');
                end%if
            end%ii
        end%function
        function ToggleVisible(this,handle)
            if strcmp(handle.Visible,'on') == 1
                handle.Visible = 'off';
            else
                handle.Visible = 'on';
            end%if
        end%function
        function ToggleRefresh(this,varargin)
            %Instruct the object instance which graphics to refresh
            %whenever the Refresh routine is called.
            for ii = 1:(nargin - 1)
                if strcmpi(varargin{ii},'Curve') == 1
                    this.refresh.Curve = ~(this.refresh.Curve == true);
                end%if
                if strcmpi(varargin{ii},'Normals') == 1
                    this.refresh.Normals = ~(this.refresh.Normals == true);
                end%if
                if strcmpi(varargin{ii},'Centroid') == 1
                    this.refresh.Centroid = ~(this.refresh.Centroid == true);
                end%if
                if strcmpi(varargin{ii},'VertexLabels') == 1
                    this.refresh.VertexLabels = ~(this.refresh.VertexLabels == true);
                end%if
                if strcmpi(varargin{ii},'SegmentLabels') == 1
                    this.refresh.SegmentLabels = ~(this.refresh.SegmentLabels == true);
                end%if
            end%ii
        end%function
        
        %Graphical Refresh routines.
        function Refresh(this,varargin)
            tic
            %Set all state of date flags to false. These are all set back
            %to true by the refresh routines if the refresh flag for each
            %curve is set to true.

            this.updated.Curve = false;
            this.updated.Normals = false;
            this.updated.Centroid = false;
            this.updated.VertexLabels = false;
            this.updated.SelfIntersect = false;
            this.updated.SegmentLabels = false;

            %Refresh routines will only execute if their respective
            %graphics' state of refresh is set to true.
            if this.refresh.Curve 
                RefreshCurve(this);
            end%if
            if this.refresh.Centroid
                RefreshCentroid(this);
            end%if
            if this.refresh.Normals
                RefreshNormals(this);
            end%if
            if this.refresh.VertexLabels
                RefreshVertexLabels(this);
            end%if
            if this.refresh.SegmentLabels
                RefreshSegmentLabels(this);
            end%if   
            if this.refresh.SelfIntersect
                RefreshSelfIntersect(this);
            end%if
            time = 1/this.refresh_rate - toc;
            pause(0 + time*(time > 0));
            drawnow
        end%function
        function RefreshCurve(this)
            %This refreshes the "main" curve. All internal edges.
            for ii = 1:this.sides
                this.sketches.Curve.XData(ii) = this.XY(ii,1);
                this.sketches.Curve.YData(ii) = this.XY(ii,2);
            end%ii
            
            %Update the last edge.
            this.sketches.LastEdge.XData(1) = this.XY(this.sides,1);
            this.sketches.LastEdge.YData(1) = this.XY(this.sides,2);
            this.sketches.LastEdge.XData(2) = this.XY(1,1);
            this.sketches.LastEdge.YData(2) = this.XY(1,2);
                       
            this.updated.Curve = true;
        end%function
        function RefreshCentroid(this)
            this.sketches.Cetroid.XData = this.xc;
            this.sketches.Cetroid.YData = this.yc;
            this.updated.Centroid = true;
        end%function
        function RefreshNormals(this)
            %Updates the graphics for the Normals with the current values
            %stored in the instance.
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        step = 1;
                        finish = this.sides - 1;
                    case 2
                        start = this.sides;
                        step = 1 - this.sides ;
                        finish = 2;
                end%switch
                for ii = start:step:finish
                    iip1 = ii + step;
                    this.sketches.Normals.UData(ii) = this.nxy(ii,1);
                    this.sketches.Normals.VData(ii) = this.nxy(ii,2);
                    this.sketches.Normals.XData(ii) = 0.5*(this.XY(iip1,1) + this.XY(ii,1));
                    this.sketches.Normals.YData(ii) = 0.5*(this.XY(iip1,2) + this.XY(ii,2));
                end%ii
            end%kk
            this.updated.Normals = true;
        end%function
        function RefreshSelfIntersect(this)
            if this.SI == 0
                return;
            end%if
            this.sketches.SelfIntersect.XData = this.xSI;
            this.sketches.SelfIntersect.YData = this.ySI;
            this.updated.SelfIntersect = true;
        end%function
        function RefreshVertexLabels(this)
            for ii = 1:this.sides
                this.sketches.VertexLabels(ii).Position(1) = this.XY(ii,1);
                this.sketches.VertexLabels(ii).Position(2) = this.XY(ii,2);
            end%ii
            this.updated.VertexLabels = true;
        end%function
        function RefreshSegmentLabels(this)

            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        step = 1;
                        finish = this.sides - 1;
                    case 2
                        start = this.sides;
                        step = 1 - this.sides;
                        finish = 2;
                end%switch
                for ii = start:step:finish
                    iip1 = 1 + step;
                    
                end%ii
            end%kk
            this.updated.SegmentLabels = true;
        end%function        
        
        %Graphics generation functions.
        function GenerateCurve(this)
            GenerateDefaultCanvas(this);
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',zeros(1,this.sides),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(1,this.sides),... %MATLAB won't allow ANY property access otherwise.
                'Linewidth',this.linewidth,...
                'LineStyle','-',...
                'DisplayName',this.name,...
                'Visible','off');
            this.sketches.LastEdge= line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',[this.XY(this.sides,1),this.XY(1,1)],...
                'YData',[this.XY(this.sides,2),this.XY(1,2)],...
                'Linewidth',this.linewidth,...
                'LineStyle','-',...
                'Visible','off');
            if this.open
                this.sketches.LastEdge.LineStyle = '--';
            end%if
            this.generated.Curve = true;
        end%function
        function GenerateCentroid(this)
            %Draw the polygon's centroid.
            GenerateDefaultCanvas(this);
            this.sketches.Centroid = line(...
                this.canvas,...
                this.xc,...
                this.yc,...
                'Marker','x',...
                'MarkerFaceColor',this.color,...
                'MarkerEdgeColor',[0,0,0],...
                'Visible','off',...
                'DisplayName',[this.name,' centroid'],...
                'Visible','off');
            this.generated.Centroid = true;
        end%function
        function GenerateNormals(this)
            %Draw the normal vectors at the midpoint of the edges of the
            %polygon.
            GenerateDefaultCanvas(this);
            this.sketches.Normals = quiver(...
                this.canvas,...
                zeros(this.sides,1),...
                zeros(this.sides,1),...
                zeros(this.sides,1),... %Normals must be precomputed prior to... 
                zeros(this.sides,1),... %... drawing this quiver object.
                'AutoScale','on',...
                'Color',this.color,...
                'Visible','off',...
                'DisplayName',[this.name,' Orientation'],...
                'Visible','off');     
            this.generated.Normals = true;
        end%function
        function GenerateSelfIntersections(this)
            %If the polyline is complex (has self-intersections) draw
            %markers where the intersections occur.
            GenerateDefaultCanvas(this);
            this.sketches.SelfIntersect = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',this.xSI,... 
                'YData',this.ySI,... 
                'Linewidth',1,...
                'LineStyle','-',...
                'Marker','d',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',this.color,...
                'DisplayName',[this.name, 'Self-Intersections'],...
                'Visible','off');
            this.generated.SelfIntersect = true;
        end%function
        function GenerateVertexLabels(this)
            %WARNING: THIS ROUTINE IS GRAPHICS AND MEMORY INTENSIVE WHEN
            %WORKING WITH LARGE POLYGONS. IF THIS ROUTINE IS DESIRED, IT IS
            %RECOMMENDED THAT IT BE USED ONLY AS A POST-PROCESSING
            %OPERATION. DO NOT USE IN THE MIDDLE OF OTHER INTENSIVE
            %CALCULATIONS OR YOU WILL GET TERRIBLE PERFORMANCE AND WAIT
            %TIMES WHILE RENDERING OR SAVING MATLAB FIGURES.
            
            %Prints labels for the vertices and their coordinates.
            this.sketches.VertexLabels = text(...
                this.XY(:,1),...
                this.XY(:,2),...
                '');
            for ii = 1:this.sides
                this.sketches.VertexLabels(ii).Color = this.color;
                this.sketches.VertexLabels(ii).Interpreter = 'latex';
                this.sketches.VertexLabels(ii).String =['V$_{',num2str(ii),'}$'];
                this.sketches.VertexLabels(ii).Visible = 'off';
            end%ii
            this.generated.VertexLabels = true;
        end%function
        function GenerateSegmentLabels(this)
            %WARNING: THIS ROUTINE IS GRAPHICS AND MEMORY INTENSIVE WHEN
            %WORKING WITH LARGE POLYGONS. IF THIS ROUTINE IS DESIRED, IT IS
            %RECOMMENDED THAT IT BE USED ONLY AS A POST-PROCESSING
            %OPERATION. DO NOT USE IN THE MIDDLE OF OTHER INTENSIVE
            %CALCULATIONS OR YOU WILL GET TERRIBLE PERFORMANCE AND WAIT
            %TIMES WHILE RENDERING OR SAVING MATLAB FIGURES.
            
            %Prints labels for the segments.
            this.sketches.SegmentLabels = text(...
                zeros(this.sides,1),...
                zeros(this.sides,1),...
                '');
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        step = 1;
                        finish = this.sides - 1;
                    case 2
                        start = this.sides;
                        step = 1 - this.sides;
                        finish = 2;
                end%switch
                for ii = start:step:finish
                    iip1 = ii + step;
                    this.sketches.SegmentLabels(ii).Interpreter = 'latex';
                    this.sketches.SegmentLabels(ii).Color = this.color;
                    this.sketches.SegmentLabels(ii).String = ['$S_{',num2str(ii),'}$'];
                    this.sketches.SegmentLabels(ii).Position(1) = (this.XY(iip1,1) + this.XY(ii,1))*0.5;
                    this.sketches.SegmentLabels(ii).Position(2) = (this.XY(iip1,2) + this.XY(ii,2))*0.5;
                    this.sketches.SegmentLabels(ii).Visible = 'off';
                end%ii
            end%kk
            this.generated.VertexLabels = true;
        end%function
        function GenerateSelfIntersect(this)
            if this.simple
                warning('This polygon is flagged as simple. Are you sure it self intersects?');
                return;
            end%if
            this.sketches.SelfIntersect = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',this.xSI,... 
                'YData',this.ySI,... 
                'Linewidth',1,...
                'LineStyle','-',...
                'Marker','d',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',this.color,...
                'DisplayName',[this.name, 'Self-Intersections'],...
                'Visible','off');
            this.generated.SelfIntersect = true;
        end%function
        function GenerateDefaultCanvas(this)
            %A routine that is called by the above graphics generation
            %routines to create a canvas in case one has not been preset
            %already.
            if ~this.canvas_set
                warning('Default canvas created!');
                ax = custom_axis;
                this.SetCanvas(ax);
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
        

    end%methods (Graphics)
    %Unit Tests
    methods (Static)
        function pass = Test0
            %Tests the area measurement of the shoelace routine. This test
            %is a replica of "Mathologer"'s video on the subject.
            %The answer is 55.
            polygon.CleanSlate;
            XY = [...
                +4,+4;...
                +0,+1;...
                -2,+5;...
                -6,+0;...
                -1,-4;...
                +5,-2];
            sides = 6;
            poly = polygon.CreateFromList(sides,XY);
            pass = poly.area == 55;
        end%function
        function pass = Test1
            polygon.CleanSlate;
            XY = [...
                -4,+6;...
                +0,+2;...
                +2,+5;...
                +7,+0;...
                +5,-6;...
                +3,+3;...
                +0,-5;...
                -6,+0;...
                -2,+1];
            %XY(:,1) = flip(XY(:,1));
            %XY(:,2) = flip(XY(:,2));

            sides = 9;
            
            ax = custom_axis;
            axis(ax,'equal');
            poly = polygon.CreateFromList(sides,XY);

            A = zeros(1,poly.sides);
            A = polygon.Angles(poly.sides,poly.XY,A,poly.orientation);
            
            R = 0.4;
            circles = cell(sides,1);
            sectors = cell(sides,1);
            
            for kk = 1:2
                switch kk
                    %case 1 %First Vertex
                    %    start = 1;
                    %    next = 1 - sides;
                    %    finish = 2;
                    case 1 %Middle Vertices
                        start = 1;
                        next = 1;
                        finish = sides - 1;
                    case 2 %Last Vertex
                        start = sides;
                        next = 1 - sides;
                        finish = 2;
                end%switch
                for ii = start:next:finish
                    iip1 = ii + next;
                    [ii,iip1,kk]
                    circles{ii} = circle.CreateXYR(XY(ii,1),XY(ii,2),R);
                    circles{ii}.SetCanvas(ax);
                    %circles{ii}.Show;
                    dx = poly.XY(iip1,1) - poly.XY(ii,1);
                    dy = poly.XY(iip1,2) - poly.XY(ii,2);
                    th1 = atan2(dx,dy);
                    th2 = th1 + A(ii);
                    sectors{ii} = CircularSector.CreateByCuttingCircle(circles{ii},th1,th2);
                    sectors{ii}.InheritFromCircle;
                    sectors{ii}.GenerateArcPolyline;
                    sectors{ii}.Polyline.Toggle('LastEdge');
                    text(ax,XY(ii,1),XY(ii,2),[num2str(A(ii)*180/pi),'$^{\circ}$'],'interpreter','latex');
                end%ii
            end%kk
                    
            poly.SetCanvas(ax);
            poly.Toggle('Curve');
            poly.Toggle('VertexLabels');

            xprt = polygon.ExportToGmsh(poly.sides,poly.XY,false,'Test');
            
            pass = false;
        end%function
        
        function [ax,poly,poly2] = TestOpenPolygon
            polygon.CleanSlate
            ax = custom_axis;
            axis(ax,'equal');
            load('fiber_copy.mat');
            fiber = 25;
            XY = fiber_copy{fiber};

            
            %Open fibers = 26,27,28
            %Closed fibers = 20,21,22,23,24,25
            
            poly = polygon.CreateFromList(length(XY),XY,true);
            poly.SetRefreshRate(4);
            poly.SetCanvas(ax);
            poly.Toggle('Curve');
            poly.Toggle('Normals');

            lambda = 0.25;
            N = 7;
            passes = 20;
            input('Ready?')
            for ii = 1:passes
                poly.Smooth(lambda,N);
                title(ax,['Pass #',num2str(ii)])
            end%ii
            
            poly2 = polygon.CreateFromOffset(poly,0.09);
            poly2.SetCanvas(ax);
            poly2.SetColor([1,0,0]);
            poly2.Toggle('Curve');
            
            poly3 = polygon.CreateFromOffset(poly,-0.09);
            poly3.SetCanvas(ax);
            poly2.SetColor([0,0,1]);
            poly3.Toggle('Curve');

            p1 = poly.sketches.Curve;
            p2 = poly2.sketches.Curve;
            p3 = poly3.sketches.Curve;
            lgd = legend(ax,[p1,p2,p3],{'Progenitor','Inward Offset','Outward Offset'},'interpreter','latex','location','northeast');
            title(lgd,'LEGEND','interpreter','latex');
            
            
            drawnow
            %poly2.
    
            %set(ax,'XLim',[12.2,14.2]);
            %set(ax,'YLim',[1.8,3.6]);
            %set(ax,'XLim',[2,5]);
            %set(ax,'YLim',[1.5,5.5]);

        end%function        
        function [ax,polygons] = TestRegularPolygons
            %Create instances of regular polygons, name them, and show their
            %normals.            
            polygon.CleanSlate;
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
                polygons{ii}.open = true;
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
        function [ax,poly] = TestStarPolygons
            %Create instances of regular concave polygons and offset them.
            polygon.CleanSlate
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
        function [ax,poly] = TestMeasurements
            %Demonstrate Measurement of metric properties.
            polygon.CleanSlate;
            ax = custom_axis;
            axis(ax,'equal')
            input('Ready?')

            load('fiber_copy.mat');
            
            fiber_number = 28;
            fiber = fiber_copy{fiber_number};
            fiber(:,1) = fiber(:,1) - min(fiber(:,1));
            fiber(:,2) = fiber(:,2) - min(fiber(:,2));
            
            %Reverse orientation
            fiber(:,1) = flip(fiber(:,1));
            fiber(:,2) = flip(fiber(:,2));

            %ax.XLim = [-0.1,1]
            %ax.YLim = [-0.2,0.6]
            poly = polygon.CreateFromList(length(fiber),fiber);
            %poly = polygon.CreateSimpleStar(0,0,10,3,4);
            
            %polygon.CreateSimpleStar(x0,y0,pairs,ri,ro)
            
            poly.SetCanvas(ax);
            poly.Show;
            poly.GenerateVertexLabels;
            poly.Toggle('VertexLabels');
            %ax.XLim = [-0.0712,0.4309];
            %ax.YLim = [-0.0868,0.2853];
            poly.sketches.Curve.Marker = 'o';
            poly.sketches.Curve.MarkerFaceColor = [0,0,0];
            poly.sketches.Curve.MarkerSize = 4;

            
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
                'Color',[0,0,1],...
                'Linewidth',1,...
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
                    
                    pause(0.001);
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
            %poly.GenerateVertexLabels;
            %poly.Toggle('VertexLabels')
        end%function
        function [ax,poly] = TestNoiseSmoothing
            %Create a regular polygon. Add randomized noise to its coordinates
            %and then apply several passes of Laplacian Smoothing.
            polygon.CleanSlate;
            
            
            lambda = 0.51; %Scaling factor.
            stencil = 3; %Number of stencil points.
            %{
            lambda = 0.35; %Scaling factor.
            stencil = 3; %Number of stencil points.
            
            %}
            
            ax = custom_axis;
            ax.Color = [1,1,1];
            axis(ax,'equal');
            ax.XTick = [];
            ax.YTick = [];

            %Create a regular polygon of many sides and show it.
            load('fiber_copy.mat');
            fiber = 13
            %fiber = 13;
            %fiber = 12;
            %fiber = 43;
            %fiber = 38
            %fiber = 34;
            XY = fiber_copy{fiber};
            
            %{
            load('noise.mat');
            XY = noise;
            %}
            poly = polygon.CreateFromList(length(XY),XY);

            %poly = polygon.CreateRegularByLength(0,0,80,2);
            poly.SetCanvas(ax);
            poly.SetColor([0,0,1]);
            poly.SetRefreshRate(20);
            poly.Show;
            poly.sketches.Curve.LineWidth = 2;
            poly.Toggle('Normals');
            %poly.GenerateVertexLabels;
            %poly.Toggle('VertexLabels');

            %poly.Toggle('Centroid');
            poly.ToggleAABB;
            set(ax,'XLim',[min(poly.XY(:,1)),max(poly.XY(:,1))]);
            set(ax,'YLim',[min(poly.XY(:,2)),max(poly.XY(:,2))]);
            
            
            %Add random noise to the regular polygon
            %title(ax,'Proceed to distort the Polygon?');
            input('Add Noise?');
            distortions = 2;
            smoothings = 100;
            
            spacer = '----------';
            fprintf('%10s\t%10s\t%10s\t%10s\t%10s\t%10s \n','Area','Perimeter','Centroid-X','Centroid-Y','Var X', 'Var Y');
            fprintf('%10s\t%10s\t%10s\t%10s\t%10s\t%10s \n',spacer,spacer,spacer,spacer,spacer,spacer);

            %{
            for ii = 1:distortions
                title(ax,['#Distortions = ',num2str(ii)]);
                poly.AddNoise(0.75);
                fprintf('%10f\t%10f\t%10f\t%10f \n',poly.area,poly.perimeter,poly.xc,poly.yc);
            end%ii
            %}
            
            
            
            %title(ax,'Proceed to smooth the Polygon?');
            input('Smooth-out the polygon?');
            kk = 0;
            stamp = Imprint(poly);

            for ii = 1:smoothings
                kk = kk + 1;
                summary = describe(poly.XY(2:end,:)- poly.XY(1:end-1,:));
                
                input('Next?');
                title(ax,['#Smoothings = ',num2str(ii)]);
                poly.Smooth(lambda,stencil);
                if mod(kk,10) == 0
                    stamp = Imprint(poly);
                    %{
                    text(ax,...
                        (max(poly.XY(:,1))),...
                        (max(poly.XY(:,2)) + min(poly.XY(:,2)))*0.5,...
                        num2str(ii));
                    %}
                end
                fprintf('%10f\t%10f\t%10f\t%10f\t%10f\t%10f \n',poly.area,poly.perimeter,poly.xc,poly.yc, summary.var(1),summary.var(2));
                poly.Refresh
            end%ii
        end%function
        function [ax,poly] = TestStamping
            %This test demonstrates the stamping feature along with some of
            %the affine transformations. First, a regular polygon is
            %defined. It is first stamped, then, it will be subject to
            %concurrent translation and rotation transformations while it
            %generates stamps.
            polygon.CleanSlate;
            
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
        function [ax,poly] = TestPointInclusion
            %This works only for rectangles.
            polygon.CleanSlate
            
            ax = custom_axis;
            poly = polygon.CreateRegularByLength(...
                0,...
                0,...
                4,...
                2);
            poly.SetCanvas(ax)
            poly.Show;
            
            xp = 2*rand;
            yp = 2*rand;
            line(...
                'Parent',ax,...
                'XData',xp,...
                'YData',yp,...
                'Marker','o',...
                'MarkerFaceColor',[0,0,1],...
                'MarkerEdgeColor',[0,0,0]);
            if polygon.ContainmentConvexPolygonvsPoint(poly.sides,poly.XY,xp,yp)
                fprintf('Point inside polygon!\n');
            else
                fprintf('Point outside polygon!\n');
            end%if
            
        end%function
        function [ax,poly,poly2] = TestOffset
            polygon.CleanSlate;
            ax = custom_axis;
            poly =  polygon.CreateRegularByLength(...
                0,...
                0,...
                10,...
                2);
            poly.SetCanvas(ax);
            poly.Show;
            
            offset = 10;
            
            poly2 = polygon.CreateFromList(...
                poly.sides,...
                poly.Offset(offset));
            
            %poly2 = polygon.CreateFromOffset(poly,2);
            poly2.SetCanvas(ax);
            poly2.SetColor([1,0,0])
            poly2.Show;
        end%function.
        
    end%methods (Demonstrations)
    %Low-level SPECIALIZED routines SPECIFIC TO THE CLASS.
    methods (Static)
         function xprt = ExportToGmsh(N,XY,open,filename)
            xprt = fopen([filename,'.msh'],'w');
            %Mesh Format header
            fprintf(xprt,'$MeshFormat\n');
            fprintf(xprt,'4.1 0 8\n');
            fprintf(xprt,'$EndMeshFormat\n');
            %Entities header
            fprintf(xprt,'$Entities\n');
            edges = N*(~open) + (N-1)*open;
            faces = ~open + open;
            fprintf(xprt,'%i %i %i %i\n',N,edges,faces,0);
            %Print CAD vertices.
            %<ID> <X> <Y> <Z> <n_tags> <tag 1> ... <tag_n>
            for ii = 1:N
                fprintf(xprt,'%i %f %f %f %i\n',ii,XY(ii,1),XY(ii,2),0.0,0);
            end%ii
            %Printf CAD edges
            %<ID> <minX> <minY> <minZ> <maxX> <maxY> <maxZ> <n_tags> <> ...
            %... <> <n_points = 2> <point 1> <point 2>
            
            for ii = 1:edges
                iip1 = ii + 1;
                if ii == edges
                    iip1 = 1;
                end%if
                %edgetag = ii                    %(1)
                minX = min(XY(ii,1),XY(iip1,1)); %(2)
                minY = min(XY(ii,2),XY(iip1,2)); %(3)
                %minZ = 0;                        %(4)
                maxX = max(XY(ii,1),XY(iip1,1)); %(5)
                maxY = max(XY(ii,2),XY(iip1,2)); %(6)
                %maxZ = 0;                        %(7)
                %ntags = 0;                       %(8)
                %npoints = 2;                     %(9)
                %pointID1 = ii                   %(10)
                %pointID2 = iip1                 %(11)
                %              1  2  3  4  5  6  7  8  9 10 11
                fprintf(xprt,'%i %f %f %f %f %f %f %i %i %i %i\n',ii,minX,minY,0,maxX,maxY,0,0,2,ii,iip1);
            end%ii
            if ~open
                %surface_tag = 1;    %(1)
                minX = min(XY(:,1)); %(2)
                minY = min(XY(:,2)); %(3)
                %minZ = 0;           %(4)
                maxX = max(XY(:,1)); %(5)
                maxY = max(XY(:,2)); %(6)
                %maxZ = 0;           %(7)
                %ntags = 0;          %(8)
                %ncurves = edges;    %(9)
                %              1  2  3  4  5  6  7  8  9
                fprintf(xprt,'%i %f %f %f %f %f %f %i %i ',1,minX,minY,0,maxX,maxY,0,0,edges);
                for ii = 1:edges
                    fprintf(xprt,'%i ',ii);
                end%ii
                fprintf(xprt,'\n');
            end%if
            fprintf(xprt,'$EndEntities');
            fclose(xprt);
        end%function
   
        %Blending functions (take two segments and blend points in between)
        function [XY,status,circ] = CircularFilletBlend(sides,XY,N,idx,Rf,x1,y1,x2,y2,x3,y3,x4,y4)
            %INPUTS:
            %sides: #of sides in the polygon.
            %XY: coordinates of the polygon.
            %x1,y1: coordinates of the start point of segment 1.
            %x2,y2: coordinates of the end point of segment 1.
            %x3,y3: coordinates of the start point of segment 2.
            %x4,y4: coordinates of the end point of segment 2.
            %Rf: fillet radius.
            %N: %# of discrete points to evaluate on the fillet.
            %idx: Indices into XY where to replace the points.
            
            %OUTPUTS:
            %XY: Coordinates with edited values.
            %status: Whether the fillet was successful (if it failed, the
            %data is not edited).
            
            status = true; %If it fails, it gets set to false later.
            dx1 = x2 - x1;
            dy1 = y2 - y1;
            dx2 = x4 - x3;
            dy2 = y4 - y3;
            
            %Must know if the segments are oriented roughly in the same
            %direction. If so, flip the orientation of one of them.
            if dx1*dx2 + dy1*dy2 > 0 
                [x1,x2] = polygon.swap(x1,x2);
                [y1,y2] = polygon.swap(y1,y2);
            end%if
            
            %Compute an orientation.
            ex = 0.5*(x2 - x1 - (x4 - x3));
            ey = 0.5*(y2 - y1 - (y4 - y3));
            [xc,yc,t1,t2,th1,th2] = circle.CenterFilletSegments(...
                x1,y1,...
                x2,y2,...
                x3,y3,...
                x4,y4,...
                Rf,ex,ey);
            if t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1
                %If this first fillet attempt failed, it could be because the
                %first guess for "ex" and "ey" was wrong. The fillet is
                %recomputed with orientation reversed.
                ex = -ex;
                ey = -ey;
                [xc,yc,t1,t2,th1,th2] = circle.CenterFilletSegments(...
                    x1,y1,...
                    x2,y2,...
                    x3,y3,...
                    x4,y4,...
                    Rf,ex,ey);
                if t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1
                    %If it failed again, there is no solution.
                    status = false;
                end%if
            end%if
            if nargout > 2
                circ = circle.CreateXYR(xc,yc,Rf);
            end%if
            if status %edit only if successful.
                XY(idx,:) = circle.polar_loop(xc,yc,Rf,th1,th2,N);
                
                %See if the las two segments intersect.
                [~,~,t1,t2] = polygon.Intersect2p2p(...
                    XY(idx(1)-1,1) ,XY(idx(1)-1,2),...
                    XY(idx(1),1)  ,XY(idx(1),2),...
                    XY(idx(N),1)  ,XY(idx(N),2),...
                    XY(idx(N)+1,1),XY(idx(N)+1,2));
                if t1 < 1 && t2 < 1 && t1 > 0 && t2 > 0
                    %Segments intersect, therefore, the blending points were
                    %assigned in reverse order. Proceed to flip the order.
                    XY(idx,1) = flip(XY(idx,1));
                    XY(idx,2) = flip(XY(idx,2));
                end%if
            end%if
            
        end%function
        function [XY,status,elli] = EllipticalFilletBlend(sides,XY,N,idx,K,x1,y1,x2,y2,x3,y3,x4,y4,idx_m)
            %idx_m index of a point towards center of the tight turn.
                        
             elli = ellipse.CreateLEE2p2d(...
                x1,y1,...
                x2-x1,y2-y1,...
                x3,y3,...
                x4-x3,y4-y3);
            %{
            elli = ellipse.CreateLeastEccentricFillet(...
                x1,y1,...
                x2-x1,y2-y1,...
                x3,y3,...
                x4-x3,y4-y3);
            %}
            status = true;
            
            %idx_m = round(0.5*(idx(1) + idx(N))); %index of some middle point.
            dxc = elli.xc - XY(idx_m,1);
            dyc = elli.yc - XY(idx_m,2);
            if dxc*elli.eax + dyc*elli.eay > 0
                %elli.ReverseMajorRadius;
                %elli.flip_a;
            end%if

            %OLD
            %{
            A = [XY(idx(1),1),XY(idx(1),2)];
            E = [XY(idx(N),1),XY(idx(N),2)];
            XY(idx,:) = elli.sector2P(A,E,N);
            %}
            
            start = idx(1)-1;
            finish = idx(N)+1;
            %A = [XY(start,1),XY(start,2)];
            %E = [XY(finish,1),XY(finish,2)];
            %temp = elli.sector2P(A,E,N+2);
            temp = elli.CreatePolarSector2P(...
                XY(start,1),XY(start,2),...
                XY(finish,1),XY(finish,2),...
                N+2);
      
            
            %XY([start,idx,finish],:) = elli.sector2P(A,E,N+2);
            XY(idx,:) = temp(2:end-1,:);
  
            %See if the last two segments intersect.
            [~,~,t1,t2] = polygon.Intersect2p2p(...
                XY(idx(1)-1,1),XY(idx(1)-1,2),...
                XY(idx(1),1)  ,XY(idx(1),2),...
                XY(idx(N),1)  ,XY(idx(N),2),...
                XY(idx(N)+1,1),XY(idx(N)+1,2));
            if t1 < 1 && t2 < 1 && t1 > 0 && t2 > 0
                %Segments intersect, therefore, the blending points were
                %assigned in reverse order. Proceed to flip the order.
                XY(idx,1) = flip(XY(idx,1));
                XY(idx,2) = flip(XY(idx,2));
            end%if
            
        end%function
        
        %Containment tests (Whether polygon encloses another object).
        function contains = ContainmentConvexPolygonvsPoint(sides,XY,xp,yp)
            %Assumes the polygon abstraction is used to represent a closed
            %curve. It also assumes that the polygon is convex, for which
            %an easier inclusion test than raycasting is available.
            contains = true; %"Innocent until proven guilty."
            for kk = 1:3
                switch kk
                    case 1
                        start = 1;
                        back = 1 - sides;
                        next = 1;
                        finish = 1;
                    case 2
                        start = 2;
                        back = 1;
                        next = 1;
                        finish = sides - 1;
                    case 3
                        start = sides;
                        back = 1;
                        next = 1-sides;
                        finish = 2;
                end%switch
                for ii = start:next:finish
                    %This sucks, because MATLAB cannot reuse the pointer
                    %for the edge that can be reused. Some edges are
                    %computed twice unnecessarily.
                    iim1 = ii - back;
                    iip1 = ii + next;
                    
                    %Vector from corner of polygon to next corner.
                    dxf = XY(iip1,1) - XY(ii,1);
                    dyf = XY(iip1,2) - XY(ii,2);
                    
                    %Vector from corner of polygon to previous corner.
                    dxb = XY(iim1,1) - XY(ii,1);
                    dyb = XY(iim1,2) - XY(ii,2);
                    
                    %Vector from corner of polygon to point.
                    dxp = xp - XY(ii,1);
                    dyp = yp - XY(ii,2);

                    if dxf*dxp + dyf*dyp < 0 || dxb*dxp + dyb*dyp < 0
                        contains = false;
                        break;
                    end%if
                end%ii
            end%kk
        end%function
        function contains = ContainmentConcavePolygonvsPoint(sides,XY,xp,yp)
        end%function
        
        %Standalone measuring functions.
        function A = Area(sides,XY)
            %Standalone computation of the polygon's area.
            A = 0; %Initialize running buffer for the area.
            for ii = 1:(sides - 1) %Loop computes area as if elements were quads.
                iip1 = ii + 1;
                A = A + XY(ii,1)*XY(iip1,2) - XY(iip1,1)*XY(ii,2);
            end%for            
            A = 0.5*A; %Halve the area to compute as triangles.
        end%function
        function P = Perimeter(sides,XY)
            %Standalone computation of the polygon's perimeter.
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
        function [A,P] = AreaPerimeter(sides,XY)
            %Compute the Area and Perimeter simultaneously.
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
        function L = Lengths(sides,XY,L)
            %Standalone routine to measure the lengths of edges. Useful for
            %when one needs to keep track of individual lengths of segments
            %of the polyline. Must provide a properly allocated buffer "L"
            %and buffer "XY".
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        next = 1;
                        finish = sides - 1;
                    case 2
                        start = sides;
                        next = 1 - sides;
                        finish = 0;
                end%switch
                for ii = start:next:finish
                    iip1 = ii + next;
                    dx = XY(iip1,1) - XY(ii,1);
                    dy = XY(iip1,2) - XY(ii,2);
                    L(ii) = sqrt(dx*dx + dy*dy);
                end%ii
            end%kk
        end%function
        function N = Normals(sides,XY,N)
            %Given a buffer "N" that is column-oriented, computes the unit
            %normals from scratch.
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        next = 1;
                        finish = sides - 1;
                    case 2
                        start = sides;
                        next = 1 - sides;
                        finish = 2;
                end%kk
                for ii = start:next:finish
                    iip1 = ii + next;
                    [N(ii,1),N(ii,2)] = UnitNormal2D(...
                        XY(iip1,1) - XY(ii,1),...%Change in X
                        XY(iip1,2) - XY(ii,2)... %Change in Y
                        );
                end%ii
            end%kk
        end%function
        function N = NormalsReuseL(sides,XY,N,L)
            %Computes the unit normals just like before, but utilizes the
            %precomputed lengths to save effort.
            for kk = 1:2
                switch kk
                    case 1
                        start = 1;
                        next = 1;
                        finish = sides - 1;
                    case 2
                        start = sides;
                        next = 1 - sides;
                        finish = 2;
                end%kk
                for ii = start:next:finish
                    iip1 = ii + next;
                    N(ii,1) = (XY(ii,2) - XY(iip1,2))/L(ii);
                    N(ii,2) = (XY(iip1,1) - XY(ii,1))/L(ii);
                end%ii
            end%kk
        end%function
        function A = Angles(N,XY,A,winding)
            %Standalone internal angle measurement routine. Provide buffers
            %"XY" of size N by 2, and "A" of size "N". Computes the Angles
            %from scratch. By default, the routine assumes that the winding
            %order is Clockwise (CW). Counterclockwise (CCW) windings are
            %permitted. Specify CW and CCW orderings with "+1" and "-1"
            %respectively.
            
            %WARNING: Winding can only be "+1" or "-1", however, no error
            %checking is  provided here.
            if nargin < 4
                winding = +1;
            end%if
            %winding = +1 = CCW.
            %winding = -1 = CW.
            dxiim1 = XY(1,1) - XY(N,1);
            dyiim1 = XY(1,2) - XY(N,2);
            Liim1 = sqrt(polygon.Dot1d1d(dxiim1,dyiim1,dxiim1,dyiim1));
            
            ccw = winding == +1;
            cw = winding == -1;
            for kk = 1:3
                switch kk
                    case 1 %First Vertex (cw) or Last Vertex (ccw).
                        start =      1*cw + ccw*N;
                        next =       1*cw - ccw*1;
                        finish =     1*cw + ccw*2;
                    case 2 %
                        start =      2*cw + ccw*(N-1);
                        next =       1*cw - ccw*1;
                        finish = (N-1)*cw + ccw*2;
                    case 3%Last Vertex (cw) or First vertex (ccw).
                        %fprintf('kk = 3\n');
                        start =      N*cw + ccw*1;
                        next =   (1-N)*cw + ccw*(N-1);
                        finish =     2*cw + ccw*(N-2);
                end%switch
                for ii = start:next:finish
                    if kk == 1
                        fprintf('kk = 1\n');
                    end
                    iip1 = ii + next;
                    dxii = XY(iip1,1) - XY(ii,1);
                    dyii = XY(iip1,2) - XY(ii,2);
                    Lii = sqrt(polygon.Dot1d1d(dxii,dyii,dxii,dyii));
                    A(ii) = pi - acos(polygon.Dot1d1d(dxii,dyii,dxiim1,dyiim1)/(Lii*Liim1));
                    if polygon.Cross1d1d(dxii,dyii,dxiim1,dyiim1) < 0
                        A(ii) = 2*pi - A(ii);
                    end%if
                    dxiim1 = dxii;
                    dyiim1 = dyii;
                    Liim1 = Lii;
                end%ii
            end%kk
        end%function
        function A = AnglesReuseL(sides,XY,A,L)
        end%function
        function K = Curvature(N,XY,K)
            %Standalone vertex-wise area-based curvature measure. Computes
            %an approximate discrete curvature from scratch.
            
            %This is based on the paper by Julia Cufi, Agusti Revento, and
            %Carlos J Rodriz title: "Curvature for Polygons."
            %Published on December 13, 2017 on the American Mathematical
            %Monthly journal.
            dxiim1 = XY(1,1) - XY(N,1);
            dyiim1 = XY(1,2) - XY(N,2);
            Liim1 = sqrt(polygon.Dot1d1d(dxiim1,dyiim1,dxiim1,dyiim1));
            for kk = 1:2
                switch kk
                    case 1 %First Vertex
                        start = 1;
                        next = 1;
                        finish = N - 1;
                    case 2 %Last Vertex
                        start = N;
                        next = 1-N;
                        finish = 2;                        
                end%switch
                for ii = start:next:finish
                    iip1 = ii + next;
                    dxii = XY(iip1,1) - XY(ii,1);
                    dyii = XY(iip1,2) - XY(ii,2);
                    Lii = sqrt(polygon.Dot1d1d(dxii,dyii,dxii,dyii));
                    %Compute the external angle;
                    angle = acos(polygon.Dot1d1d(dxiim1,dyiim1,dxii,dyii)/(Liim1*Lii));
                    %Compute the curvature.
                    K(ii) = 2*angle/(Lii + Liim1);
                    
                    dyiim1 = dyii;
                    dxiim1 = dxii;
                    Liim1 = Lii;
                end%ii
            end%for
        end%function
        function K = CurvatureReuseAL(N,XY,K,A,L)
            %Computes the area-based curvature with the assistance of
            %precomputed Internal Angles and lengths.
        end%function
        
        
        function [SI,XY_SI,seg] = DetectSelfIntersections(sides,XY)
            %Inspect the polygon for self-intersections.
            %WARNING: THIS IS MY OWN ATTEMPT AT MAKING A SELF-INTERSECTION
            %DETECTOR. IT IS SLOW O(n^2), NOT THE BEST.
            %MUST UPGRADE TO SORT-ASSISTED SEARCH.
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
        
        %Offset points according to a sequence of normal vectors.
        function XY_out = OffsetClosed(sides,XY,nxy,offset)
            %WARNING: This routine does not check if the normal vectors
            %supplied by nxy are in fact unitary.
            XY_out = zeros(sides,2);
            for kk = 1:3
                switch kk
                    case 1 %First offset is coupled to last point.
                        start = 1;
                        back = sides -1;
                        next = 1;
                        finish = 1;
                    case 2 %Intermediate offsets require their immediate neighbors.
                        start = 2;
                        back = -1;
                        next = 1;
                        finish = (sides - 1);
                    case 3 %Last offset is coupled to first point.
                        start = sides;
                        back = -1;
                        next = 1 - sides;
                        finish = 2;
                end%switch
                for ii = start:next:finish
                    iim1 = ii + back;
                    iip1 = ii + next;
                    [XY_out(ii,1),XY_out(ii,2),~,~] = polygon.Intersect2p2d(...
                        XY(iim1,1) + nxy(iim1,1)*offset, XY(iim1,2) + nxy(iim1,2)*offset,... %First point.
                        XY(ii,1)   - XY(iim1,1),         XY(ii,2) - XY(iim1,2),... %Direction at first point.
                        XY(ii,1)   + nxy(ii,1)*offset,   XY(ii,2) + nxy(ii,2)*offset,... %Second point.
                        XY(iip1,1) - XY(ii,1),           XY(iip1,2) - XY(ii,2)); %Direction at the second point.
                end%ii
            end%kk
        end%function
        function XY_out = OffsetOpen(sides,XY,nxy,offset)
            %WARNING: This routine does not check if the normal vectors
            %supplied by nxy are in fact unitary.
            XY_out = zeros(sides,2);
            
            %The first point is manually offset according to a single
            %normal.
            [nx,ny] = polygon.UnitNormal2D(...
                XY(2,1) - XY(1,1),...
                XY(2,2) - XY(1,2));
            XY_out(1,1) = XY(1,1) + offset*nx;
            XY_out(1,2) = XY(1,2) + offset*ny;
            
            %Middle points are offset just like before (according to
            %intersections of pairs of offset lines).
            for ii = 2:(sides - 1)
                iim1 = ii - 1;
                iip1 = ii + 1;
                [XY_out(ii,1),XY_out(ii,2),~,~] = polygon.Intersect2p2d(...
                    XY(iim1,1) + nxy(iim1,1)*offset, XY(iim1,2) + nxy(iim1,2)*offset,... %First point.
                    XY(ii,1)   - XY(iim1,1),         XY(ii,2) - XY(iim1,2),... %Direction at first point.
                    XY(ii,1)   + nxy(ii,1)*offset,   XY(ii,2) + nxy(ii,2)*offset,... %Second point.
                    XY(iip1,1) - XY(ii,1),           XY(iip1,2) - XY(ii,2)); %Direction at the second point.
            end%ii
            
            %Finally, the last point is also offset according to a single
            %normal.
            [nx,ny] = polygon.UnitNormal2D(...
                XY(sides,1) - XY(sides - 1,1),...
                XY(sides,2) - XY(sides - 1,2));
            XY_out(sides,1) = XY(sides,1) + offset*nx;
            XY_out(sides,2) = XY(sides,2) + offset*ny;
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
        %DEPRECATE
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
        %DEPRECATE
        
        %Coordinate Manipulation routines
        function XY = DisplaceCoordinates(XY,sides,dim,dx)
            %A general Coordinate displacement routine. Given a coordinate
            %buffer "XY", its size (in the form of the "sides" integer), a
            %dimension "dim" to represent "N" as in "R^N" space, and an
            %array of displacement components "dx."
            for ii = 1:sides
                for jj = 1:dim
                    XY(ii,jj) = XY(ii,jj) + dx(jj);
                end%jj
            end%ii
        end%function
    end%methods (Static)
    %Low-level ELEMENTARY routines
    methods (Static)
        %Elementary 2D geometry routines.
        function XY = CircularPolarLoop(XY,sides,R,th1,th2)
            %Simple function that evaluates discrete points on the arc of a
            %circular sector. It is assumed that th2 > th1. For spanning of
            %the whole circle, set th1 = 0 and th2 = 2*pi. That is, units
            %are assumed to be radians.
            d_theta = (th2 - th1)/(sides - 1);
            theta = th1;
            for ii = 1:sides
                XY(ii,1) = R*sin(theta);
                XY(ii,2) = R*cos(theta);
                theta = theta + d_theta;
            end%ii
        end%function 
        
        %Elementary 2D vector algebra routines.
        function [nx,ny] = UnitNormal2D(dx,dy)
            %Find the unit normal vector to some direction according to a
            %hardcoded 90-degree (pi/4 radian) CCW rotation.
            mag = sqrt(dx*dx + dy*dy); %Magnitude of the input direction.
            nx = -dy/mag;
            ny = +dx/mag;
        end%function
        function [x,y,t1,t2] = Intersect2p2p(x1,y1,x2,y2,x3,y3,x4,y4)
            %Intersect two lines defined by two segments in two dimensions.
            [x,y,t1,t2] = polygon.Intersect2p2d(...
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
        function dot = Dot2p2p(x1,y1,x2,y2,x3,y3,x4,y4)
            dot = polygon.Dot1d1d(...
                x2 - x1,...
                y2 - y1,...
                x4 - x3,...
                y4 - y3);
        end%function
        function dot = Dot1d1d(dx1,dy1,dx2,dy2)
            dot = dx1*dx2 + dy1*dy2;
        end%function
        function det = Cross2p2p(x1,y1,x2,y2,x3,y3,x4,y4)
            %Computes the magnitude of the vector cross product of two
            %vectors encoded as 2-point segments
            det = polygon.Cross1d1d(...
                x2 - x1,...
                y2 - y1,...
                x4 - x3,...
                y4 - y3);
        end%function
        function det = Cross1d1d(dx1,dy1,dx2,dy2)
            %Computes the magnitude of the vector cross product of two
            %vectors contained in the 2D plane encoded by their components.
            det = dx1*dy2 - dy1*dx2;
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
        function [L,idx] = DetectSignChanges(len,array)
            %Given an array, this routines finds all the indices where a
            %sign change occurs. The indices are for the element
            %immediately before another element of different sign.
            
            %INPUTS:
            %len: Array length.
            %array: The array.
            
            %OUTPUTS:
            %L: # of sign changes detected.
            %idx: The indices into "array" immediately before a sign
            %change.
            ii = 1;%Initialize an indexer into the array.
            sign1 = polygon.Signum(array(1)); %Get the sign of the first element.
            idx_chain = dchain;
            while ii < len
                ii = ii + 1;
                sign2 = polygon.Signum(array(ii));
                if sign2 ~= sign1 && sign2 ~= 0
                    idx_chain.append(dlink(ii),'right');
                    sign1 = sign2;
                end%if
            end%while
            
            
            L = idx_chain.length;
            idx = idx_chain.chain2array;
        end%function
    end%methods (Static)
end%classdef
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
% Discrete Curvature
%%
% $\kappa_{k} = \frac{2\alpha_{k}\pi}{L_{k-1} + L_{k}}$
%% Class definition
classdef polygon < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the polygon's graphics as they appear on axes legend.
        color %Color of polygon's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the polygon.
        sketches %Structure containing handles to the object's graphics.
        refresh %Structure holding boolean values that control which graphics refresh.
        updated %Structure holding boolean values that track graphics state of date.
        generated %Structure holding boolean values that track whether certain graphics are generated.
        refresh_rate %A parameter that determines how long (ideally) calls to Refresh shoul take.
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
        
        SI %#of self-intersections.
        xSI %X-coordinates of self-intersections (SI by 1 array).
        ySI %Y-coordinates of self-intersections (SI by 1 array).
        sSI %Segments involved in self intersection (SI by 2 array).
    end    
    %FLAG and STATE variables.
    properties (Hidden = true)
        %Graphics-related flags.
        graphics_initialized %Whether "line" and or "quiver" have been generated.
        normals_computed %Whether the inward/outward normals have been computed.
        canvas_set %Whether an axes object has been initialized.
        AABB_present
        vertex_labels_on
        segment_labels_on
        last_edge %This is the line that connects the first point to the last one.
        
        valid %Whether this instance has a valid definition.
    end%properties (Hidden)    
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = polygon(varargin)
            %Defining and metric data.
            this.XY = [];%KEEP
            this.simple = [];%KEEP
            this.convex = [];%SKIP
            this.open = [];%KEEP
            this.area = [];%KEEP
            this.perimeter = [];%KEEP
            this.xc = [];%SKIP
            this.yc = [];%SKIP
            
            %Self-Intersections
            this.SI = 0; %KEEP
            this.xSI = []; %KEEP
            this.ySI = []; %KEEP
            this.sSI = []; %KEEP
            
            %Graphics related.
            this.sketches = struct(...
                'Curve',[],...
                'Normals',[],...
                'Centroid',[],...
                'SelfIntersect',[],...
                'Vertex_Labels',[],...
                'Segment_Labels',[]);
            this.refresh = struct(...
                'Curve',true,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',true,...
                'Vertex_Labels',false,...
                'Segment_Labels',false);
            this.updated = struct(...
                'Curve',false,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',false,...
                'Vertex_Labels',false,...
                'Segment_Labels',false);
            this.generated = struct(...
                'Curve',false,...
                'Normals',false,...
                'Centroid',false,...
                'SelfIntersect',false,...
                'Vertex_Labels',false,...
                'Segment_Labels',false);
            
            
            this.refresh_rate = 1;%Desired refresh rate in Hz;
            
            this.name = 'Polygon';
            this.color = [0,0,0];
            
            %Flag and state variables.
            this.vertex_labels_on = false;
            this.segment_labels_on = false;
            this.graphics_initialized = false;
            this.normals_computed = false;
            this.orientation = [];
            this.AABB_present = (exist('AABB','file') == 2);
        end%function
        
        %Custom creation routines.
        function poly = CreateRegularByLength(x0,y0,sides,length)
            if sides < 3
                error('Need atleast 3 sides for a polygon.');
            end%if
            poly = polygon;
            poly.sides = sides;
            poly.Calloc(sides);
            
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
            poly.Calloc(sides);
            
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
            poly.Measure;
            
            if poly.AABB_present
                poly.aabb = AABB.CreateFromList(2,poly.sides,poly.XY);
            end%if
            
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

            %If # of stencil points is not specified, default to 3.
            if nargin < 3
                N = 3;
            end%if
            
            %Like with open curve smoothing, a minimal copy of "trailing"
            %values is needed. The "ahead" values are not needed though.
            Nm1o2 = (N-1)*0.5; %Nm1o2 = "N minus 1 over 2"
            Np1o2 = Nm1o2 + 1; %This index corresponds to the middle of the central difference scheme.
            trail_x  = zeros(1,N-1); %Stores the trailing (N-1)/2...
            trail_y  = zeros(1,N-1); %... stencil points.
            %trail_x  = zeros(1,Nm1o2); %Stores the trailing (N-1)/2...
            %trail_y  = zeros(1,Nm1o2); %... stencil points.
            
            for ii = 1:Nm1o2 %This time, more points need to be remembered...
                trail_x(ii) = this.XY(ii,1); %...this is because towards...
                trail_y(ii) = this.XY(ii,2); %...the end backward schemes...
            end%ii                            ...are applied.
                        
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
                %this.XY(ii,1) = this.XY(ii,1) + lambda*x_diffusion;
                %this.XY(ii,2) = this.XY(ii,2) + lambda*y_diffusion;
                %{
                dx = lambda*x_diffusion;
                dy = lambda*y_diffusion;
                this.XY(ii,1) = this.XY(ii,1) + dx*0.0001;
                this.XY(ii,2) = this.XY(ii,2) + dy*0.0001;
                
                
                if sqrt(dx*dx + dy*dy) < sqrt((this.XY(ii+1,1)-this.XY(ii,1))^2 + (this.XY(ii+1,2)-this.XY(ii,2))^2)
                   this.XY(ii,1) = this.XY(ii,1) + dx;
                   this.XY(ii,2) = this.XY(ii,2) + dy;
                else
                    fprintf('Delta greater than length averted!\n');
                end%if
                %}
                
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
            
            %{
            %At this stage, all points elligible for central difference
            %have been spent. Backward difference on the last "(N-1)/2"
            %points now ensus.
            trail2_x = zeros(1,N-1);
            trail2_y = zeros(1,N-1);
            
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
                %for jj = 1:(ii - (this.sides - Nm1o2))
                for jj = 1:(Np1o2 + ii -this.sides + Nm1o2)
                    t_idx = jj + sm1; %Index into the trailing buffer.
                    t_idx = t_idx - Nm1*(Nm1 - t_idx < 0); %Push it back if it overflows.                    
                    x_diffusion = x_diffusion + trail_x(t_idx)*FD(jj);
                    y_diffusion = y_diffusion + trail_y(t_idx)*FD(jj);
                end%jj

                %Reference the points ahead
                for jj = ((Np1o2 + ii -this.sides + Nm1o2+1)):N
                %for jj =(1 + ii - (this.sides - Nm1o2)):N
                    idx = ii + sequence(jj) + 1;
                    x_diffusion = x_diffusion + this.XY(idx,1)*FD(jj);
                    y_diffusion = y_diffusion + this.XY(idx,2)*FD(jj);
                end%%jj                                
                
                
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
        function DetectSharpTurns(this,Rf)
            %This routine applies a "rolling" circle type of operation to
            %flag turns that exceed the radius of the tighest turn that the
            %polygon should make.
            
            bad_points = ones(1,this.sides); %Allocate memory for an array of flags.
            
            %For this to work, the polygon must be oriented inwards.
            if this.orientation
                Rf = -Rf;
                %this.ReverseNormals;
            end%if
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
                abs(Rf));%Radiues
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
                    xc2 = 0.5*(this.XY(iip1,1) + this.XY(iip2,1)) + 0.95*Rf*this.nxy(iip1,1);
                    yc2 = 0.5*(this.XY(iip1,2) + this.XY(iip2,2)) + 0.95*Rf*this.nxy(iip1,2);
                    
                    %Compute Arc length traversed. (ANIMATION + FUNCTIONALITY)
                    %{
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
                    needleL = needle; %Need a counter to subtract from the first needle.
                    while X_sorted(needleL) > (F.xc - F.R) %&& X_sorted(1) > (F.xc - F.R)
                        if F.ContainsPoint(this.XY(X_perm(needleL),1),this.XY(X_perm(needleL),2))
                            flags.ZData(X_perm(needleL)) = 0;%"Un-NaN" the marker.
                            bad_points(X_perm(needleL)) = -1;
                        end%if
                        needleL = needleL - 1;
                        if needleL == 0 %If even the minimum X-value is to the left, break.
                            %needleL = poly.sides;
                            break;
                        end%if
                    end%while
                    
                    %Cheack nearby "right" points.
                    needleR = needle;
                    while X_sorted(needleR) < (F.xc + F.R) %&& X_sorted(poly.sides) < (F.xc + F.R)
                        if F.ContainsPoint(this.XY(X_perm(needleR),1),this.XY(X_perm(needleR),2))
                            flags.ZData(X_perm(needleR)) = 0;%"Un-NaN" the marker.
                            bad_points(X_perm(needleR)) = -1;
                        end%if
                        needleR = needleR + 1;
                        if needleR == (this.sides + 1)
                            %needleR = 1;
                            break;
                        end%if
                    end%while
                    needle = old_needle;
                    %input('Next?')
                    %pause(0.0001);
                end%ii
            end%kk

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
                    if bad_points(ii - back) == -1 && bad_points(ii + next) == -1
                        bad_points(ii) = -1;
                        flags.ZData(ii) = 0;
                    end%if
                end%ii
                if this.open
                    break;
                end%if
            end%for

            %%%%%%%%%%%%%SECOND ATTEMPT

            %Extract bands of continuous bad flags.
            [L,idx] = polygon.DetectSignChanges(this.sides,bad_points);
            idx
            
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
            
            
            
            this.Measure;
            %%%%%%%%%%%%%SECOND ATTEMPT

            
            %%%%%%%%%%%%%FIRST ATTEMPT
            %{
            if mod(L,2) == 1 %Odd-numbered sign changes can be are due to sharp turn
                for ii = L-1:-2:1
                    iim1 = ii - 1;
                    iip1 = ii + 1;
                    points = idx(ii) - idx(iim1); %Number of faulty points.

                end%ii
                %[XY,status,circ] = CircularFilletBlend(XY,N,idx,Rf,x1,y1,x2,y2,x3,y3,x4,y4)
            else %Even number of sign changes.
                for ii = L:-2:2
                    iim1 = ii - 1;
                    iip1 = ii + 1;
                    points = idx(ii) - idx(iim1); %Number of faulty points.
                    tgt1 = idx(ii);
                    tgt2 = idx(iim1) - 1;
                    if points > 2
                        %Make this into a function.
                        %Midpoint
                        m_x = 0.5*(this.XY(tgt2,1) + this.XY(tgt1,1));
                        m_y = 0.5*(this.XY(tgt2,2) + this.XY(tgt1,2));
                        ex = this.XY(round(0.5*(tgt1 + tgt2)),1) - m_x;
                        ey = this.XY(round(0.5*(tgt1 + tgt2)),2) - m_y;
                        quiver(this.canvas,m_x,m_y,ex,ey,'Color',[1,1,1],'linewidth',1);
                        
                        blender = circle.CreateFilletWithRadius(...
                            this.XY(idx(iim1) - 1,1),this.XY(idx(iim1) - 1,2),...
                            this.XY(idx(iim1),1)    ,this.XY(idx(iim1),2),...
                            this.XY(idx(ii) - 1,1)  ,this.XY(idx(ii) - 1,2),...
                            this.XY(idx(ii),1)      ,this.XY(idx(ii),2),...
                            Rf,ex,ey);
                        %Make this into a function
                        blender.SetColor([1,1,1]);
                        blender.sketches.Curve.Color = [1,0,0];
                        blender.SetCanvas(this.canvas);
                        blender.Show;
                    end%if
                end%ii
            end%if
            %}
            %%%%%%%%%%%%%FIRST ATTEMPT
  
            
            delete(F); %Deallocate the circle
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
        
        %Measure metric properties about this polygon
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
        
        
        %Query routines
        function inside = IsInside(this,XY_in)
        
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
        function SetRefreshRate(this,rate)
            %Set a target refresh rate for the Refresh routine.
            if rate == 0
                rate = 1;
            end%if
            this.refresh_rate = rate*(-1*(rate < 0));
        end%function
        function SetName(this,string)
            this.name = string;
            if this.graphics_initialized
                this.sketches.Curve.DisplayName = string;
                this.sketches.Normals.DisplayName = [string,'''s Normals'];
                this.sketches.Centroid.DisplayName = [string,'''s Centroid'];
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
                this.last_edge.Color = RGB;
                this.sketches.Normals.Color = RGB;
                this.sketches.Centroid.MarkerFaceColor = RGB;
                if this.generated.Vertex_Labels
                    for ii = 1:this.sides
                        this.sketches.Vertex_Labels(ii).Color = RGB;
                    end%ii
                end%if
                if this.generated.Segment_Labels
                    for ii = 1:this.sides
                        this.sketches.Segment_Labels(ii).Color = RGB;
                    end%ii
                end%if
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
                'XData',zeros(1,this.sides + 1),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(1,this.sides + 1),... %MATLAB won't allow ANY property access otherwise.
                'Linewidth',1,...
                'LineStyle','-',...
                'DisplayName',this.name);
            this.last_edge = line(...
                'Parent',this.canvas,...
                'Color',this.color,...
                'XData',[this.XY(this.sides,1),this.XY(1,1)],...
                'YData',[this.XY(this.sides,2),this.XY(1,2)],...
                'Linewidth',1,...
                'LineStyle','-');
            if this.open
                this.last_edge.LineStyle = '--';
            end%if
            this.generated.Curve = true;
            RefreshCurve(this);
            
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
            this.generated.Centroid = true;

            
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
            this.generated.Normals = true;
            RefreshNormals(this);
            
            %Generate markers to flag self intersections if any.
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
                'DisplayName',[this.name, 'Self-Intersections']);
            RefreshSelfIntersect(this);
            this.generated.SelfIntersect = true;
            
            
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

            if ~isempty(this.sketches.Vertex_Labels)
                for ii = this.sides:-1:1
                    delete(this.sketches.Vertex_Labels(ii));
                end%ii
                this.generated.Vertex_Labels = false;
            end%if
            if ~isempty(this.sketches.Segment_Labels)
                for ii = this.sides:-1:1
                    delete(this.sketches.Segment_Labels(ii));
                end%ii
                this.generated.Segment_Labels = false;
            end%if
            %Update graphics flag.
            this.graphics_initialized = false;
        end%function
        
        %Annotations.
        function GenerateVertexLabels(this)
            %WARNING: THIS ROUTINE IS GRAPHICS AND MEMORY INTENSIVE WHEN
            %WORKING WITH LARGE POLYGONS. IF THIS ROUTINE IS DESIRED, IT IS
            %RECOMMENDED THAT IT BE USED ONLY AS A POST-PROCESSING
            %OPERATION. DO NOT USE IN THE MIDDLE OF OTHER INTENSIVE
            %CALCULATIONS OR YOU WILL GET TERRIBLE PERFORMANCE AND WAIT
            %TIMES WHILE RENDERING OR SAVING MATLAB FIGURES.
            
            %Prints labels for the vertices and their coordinates.
            this.sketches.Vertex_Labels = text(...
                this.XY(:,1),...
                this.XY(:,2),...
                '');
            for ii = 1:this.sides
                this.sketches.Vertex_Labels(ii).Color = this.color;
                this.sketches.Vertex_Labels(ii).Visible = 'off';
                this.sketches.Vertex_Labels(ii).Interpreter = 'latex';
                this.sketches.Vertex_Labels(ii).String = {...
                    ['V$_{',num2str(ii),'}$'],['(',num2str(this.XY(ii,1)),';',num2str(this.XY(ii,2)),')']};
            end%ii
            this.vertex_labels_on = true;            
        end%function
        function GenerateSegmentLabels(this)
            %WARNING: THIS ROUTINE IS GRAPHICS AND MEMORY INTENSIVE WHEN
            %WORKING WITH LARGE POLYGONS. IF THIS ROUTINE IS DESIRED, IT IS
            %RECOMMENDED THAT IT BE USED ONLY AS A POST-PROCESSING
            %OPERATION. DO NOT USE IN THE MIDDLE OF OTHER INTENSIVE
            %CALCULATIONS OR YOU WILL GET TERRIBLE PERFORMANCE AND WAIT
            %TIMES WHILE RENDERING OR SAVING MATLAB FIGURES.
            
            %Prints labels for the segments.
            this.sketches.Segment_Labels = text(...
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
                    this.sketches.Segment_Labels(ii).Visible = 'off';
                    this.sketches.Segment_Labels(ii).Interpreter = 'latex';
                    this.sketches.Segment_Labels(ii).Color = this.color;
                    this.sketches.Segment_Labels(ii).String = ['$S_{',num2str(ii),'}$'];
                    this.sketches.Segment_Labels(ii).Position(1) = (this.XY(iip1,1) + this.XY(ii,1))*0.5;
                    this.sketches.Segment_Labels(ii).Position(2) = (this.XY(iip1,2) + this.XY(ii,2))*0.5;
                end%ii
            end%kk
            this.segment_labels_on = true;
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
                'DisplayName',[this.name, 'Self-Intersections']);
            this.generated.SelfIntersect = true;
            
        end%function
        
        %Visibility toggling functions
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
                if strcmpi(varargin{ii},'Vertex_Labels') == 1
                    this.refresh.Vertex_Labels = ~(this.refresh.Vertex_Labels == true);
                end%if
                if strcmpi(varargin{ii},'Segment_Labels') == 1
                    this.refresh.Segment_Labels = ~(this.refresh.Segment_Labels == true);
                end%if
            end%ii
        end%function
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
                    VisibilityToggle(this,this.sketches.Curve);
                    VisibilityToggle(this,this.last_edge);
                    if ~this.updated.Curve && strcmp(this.sketches.Curve.Visible,'on') == 1
                        RefreshCurve(this);
                    end%if
                elseif strcmpi(varargin{ii},'Centroid') == 1
                    VisibilityToggle(this,this.sketches.Centroid);
                    if ~this.updated.Centroid
                        RefreshCentroid(this);
                    end%if
                elseif strcmpi(varargin{ii},'Normals') == 1
                    VisibilityToggle(this,this.sketches.Normals);
                    if ~this.updated.Normals && strcmp(this.sketches.Normals.Visible,'on') == 1
                        RefreshNormals(this);
                    end%if
                elseif strcmpi(varargin{ii},'Vertex_Labels') == 1
                    if isempty(this.sketches.Vertex_Labels)
                        return;
                    end%if
                    if ~this.updated.Vertex_Labels
                        RefreshVertexLabels(this);
                    end%if
                    for jj = 1:this.sides
                        VisibilityToggle(this,this.sketches.Vertex_Labels(jj));
                    end%jj
                elseif strcmpi(varargin{ii},'Segment_Labels') == 1
                    if isempty(this.sketches.Segment_Labels)
                        return;
                    end%if
                    if ~this.updated.Segment_Labels
                        RefreshSegmentLabels(this);
                    end%if
                    for jj = 1:this.sides
                        VisibilityToggle(this,this.sketches.Segment_Labels(jj));
                    end%jj
                elseif strcmpi(varargin{ii},'SelfIntersect') == 1
                    if isempty(this.sketches.SelfIntersect)
                        return;
                    end%if
                    if ~this.updated.Segment_Labels
                        RefreshSelfIntersect(this);
                    end%if
                else
                    warning('Unrecognizable graphics option. Valid options are: "Curve","Normals","Centroid", "Vertex_Labels", and "Segment_Labels"');
                end%if
            end
        end%function
        %NEEDS TO BE UPDATED
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
        %NEEDS TO BE UPDATED
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
        
        %Graphical Refresh routines.
        function Refresh(this,varargin)
            tic
            %Set all state of date flags to false. These are all set back
            %to true by the refresh routines if the refresh flag for each
            %curve is set to true.
            this.updated.Curve = false;
            this.updated.Normals = false;
            this.updated.Centroid = false;
            this.updated.Vertex_Labels = false;
            this.updated.Segment_Labels = false;

            %Refresh routines will only execute if their respective
            %graphics' state of refresh is set to true.
            if this.refresh.Curve 
                RefreshCurve(this);
            else
                this.Toggle('Curve');
            end%if
            if this.refresh.Centroid
                RefreshCentroid(this);
            else
                this.Toggle('Centroid');
            end%if
            if this.refresh.Normals
                RefreshNormals(this);
            else
                this.Toggle('Normals');
            end%if
            if this.refresh.Vertex_Labels
                RefreshVertexLabels(this);
            else
                this.Toggle('Vertex_Labels');
            end%if
            if this.refresh.Segment_Labels
                RefreshSegmentLabels(this);
            else
                this.Toggle('Segment_Labels');
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
            if this.open
                this.sketches.Curve.XData(this.sides + 1) = NaN;
                this.sketches.Curve.YData(this.sides + 1) = NaN;
            else
                this.sketches.Curve.XData(this.sides + 1) = this.XY(1,1);
                this.sketches.Curve.YData(this.sides + 1) = this.XY(1,2);
            end%if
            
            %Update the last edge.
            this.last_edge.XData = [this.XY(this.sides,1),this.XY(1,1)];
            this.last_edge.YData = [this.XY(this.sides,2),this.XY(1,2)];            
            this.updated.Curve = true;
        end%function
        function RefreshCentroid(this)
            
            this.sketches.Cetroid.XData = this.xc;
            this.sketches.Cetroid.YData = this.yc;
            this.updated.Centroid = true;
        end%function
        function RefreshNormals(this)
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
                this.sketches.Vertex_Labels(ii).Position(1) = this.XY(ii,1);
                this.sketches.Vertex_Labels(ii).Position(2) = this.XY(ii,2);
            end%ii
            this.updated.Vertex_Labels = true;
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
            this.updated.Segment_Labels = true;
        end%function
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
        
        function [ax,poly,poly2] = TestOpenPolygon
            polygon.CleanSlate
            ax = custom_axis;
            axis(ax,'equal');
            load('fiber_copy.mat');
            fiber = 28;
            XY = fiber_copy{fiber};

            
            %Open fibers = 26,27,28
            %Closed fibers = 20,21,22,23,24,25
            
            poly = polygon.CreateFromList(length(XY),XY,true);
            poly.SetRefreshRate(4);
            poly.SetCanvas(ax);
            poly.Show;
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
            poly2.Show;
   
            
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
            
            fiber_number = 37;
            fiber = fiber_copy{fiber_number};
            fiber(:,1) = fiber(:,1) - min(fiber(:,1));
            fiber(:,2) = fiber(:,2) - min(fiber(:,2));
            
            %Reverse orientation
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
            %poly.Toggle('Vertex_Labels')
        end%function
        function [ax,poly] = TestSmoothing

        end%function
        function [ax,poly] = TestNoiseSmoothing
            %Create a regular polygon. Add randomized noise to its coordinates
            %and then apply several passes of Laplacian Smoothing.
            polygon.CleanSlate;
            
            lambda = 0.35; %Scaling factor.
            stencil = 5; %Number of stencil points.
            
            
            ax = custom_axis;
            ax.Color = [0,0,0];
            axis(ax,'equal');
            
            %Create a regular polygon of many sides and show it.
            poly = polygon.CreateRegularByLength(0,0,80,2);
            poly.SetCanvas(ax);
            poly.SetColor([0,1,0]);
            poly.Show;
            poly.Toggle('Normals');
            
            set(ax,'XLim',[min(poly.XY(:,1)),max(poly.XY(:,1))]);
            set(ax,'YLim',[min(poly.XY(:,2)),max(poly.XY(:,2))]);
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
                poly.AddNoise(0.75);
                fprintf('%10f\t%10f\t%10f\t%10f \n',poly.area,poly.perimeter,poly.xc,poly.yc);
                pause(0.05);
            end%ii
            
            title(ax,'Proceed to smooth the Polygon?');
            input('Smooth-out the polygon?');
            for ii = 1:smoothings
                title(ax,['#Smoothings = ',num2str(ii)]);
                poly.Smooth(lambda,stencil);
                fprintf('%10f\t%10f\t%10f\t%10f \n',poly.area,poly.perimeter,poly.xc,poly.yc);
            end%ii
        end%function
        function [ax,poly] = TestDispersion
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

    end%methods (Demonstrations)
    %Low-level SPECIALIZED routines SPECIFIC TO THE CLASS.
    methods (Static)
        
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
        
        %Standalone measuring functions.
        function A = Area(sides,XY)
            %Standalone computation of the polygon's area.
            A = 0; %Initialize running buffer for the area.
            for ii = 1:(sides - 1) %Loop computes area as if elements were quads.
                iip1 = ii + 1;
                A = A + XY(ii,1)*XY(iip1,2) - XY(iip1,1)*XY(ii,2);
            end%for            
            A = A/2; %Halve the area to compute as triangles.
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
        function K = Curvature(sides,XY,K)
            %sides: Number of sides in the polygon.
            %XY: Coordinates of the polygon.
            %K: Discrete curvature measure.
            
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
    %Low-level ELEMENTARY routines
    methods (Static)
        %Basic Geometry features.
        function [nx,ny] = UnitNormal2D(dx,dy)
            %Find the unit normal vector to some direction according to a
            %right-handed coordinate system.
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
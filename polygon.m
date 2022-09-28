classdef polygon < handle
    properties (SetAccess = public)
        canvas
        sketches
        
        sketch %DEPRECATE
        sketch_center %DEPRECATE
        sketch_labels %DEPRECATE
    end%properties
    properties (SetAccess = private)
        sides %Number of sides (and coordinates).
        XY %XY coordinates (1st column is for x-values, 2nd for y-values.)
    end%properties
    properties (SetAccess = protected)
        aabb %Axis-Aligned bounding box of this polygon.
        nxy %XY components of the normal vectors at each segment.
        centroid %
        area %Area enclosed by the polygon.
        P %Perimeter subtended by polygon.
        simple %Flag to denote whether the polygon self-intersects.
        convex %Flag to denote whether the polygon is convex.
    end
    properties (Hidden = true)
        %Graphics-related flags.
        graphics_initialized
        normals_computed
        canvas_set
        
        updated
        valid
    end%properties (Hidden)
    
    %High-level functions that MODIFY specific instances of the class
    %object.
    methods
        %Constructor
        function this = polygon(varargin)
               this.XY = [];
               this.simple = [];
               this.convex = [];
               this.area = [];
               this.P = [];
               this.sketch = [];
               this.sketches = struct(...
                   'Curve',[],...
                   'Normals',[],...
                   'Vertex_Labels',[],...
                   'Segment_Labels',[]);
               this.graphics_initialized = false;
               normals_computed = false;
        end%function
        
        %Compute the normals of this polygon
        function ComputeNormals(this)
            if ~this.valid
                warning('This polygon is flagged as invalid.');
                return;
            end%if
            if isempty(this.XY)
                warning('Polygon has no associated coordinates!');
                return;
            end%if
            this.nxy = polygon.Normals2D(this.sides,this.XY); %Magic happens here.
            if this.graphics_initialized %Update the quiver plot.
                this.sketches.Normals.UData = this.nxy(:,1);
                this.sketches.Normals.VData = this.nxy(:,2);
            end%if
            this.normals_computed = true;
        end%function
        
        %Reverse the direction of the normal vectors.
        function ReverseNormals(this)
            if this.normals_computed
                for ii = this.sides
                    this.nxy(ii,1) = (-1)*this.nxy(ii,1);
                    this.nxy(ii,2) = (-1)*this.nxy(ii,2);
                end%ii
                if this.graphics_initialized
                    for ii = 1:this.sides
                        this.sketches.Normals.UData(ii) = -this.sketches.Normals.UData(ii);
                        this.sketches.Normals.VData(ii) = -this.sketches.Normals.VData(ii);
                    end%ii
                end%if
            else
                warning('Normal vectors NOT computed. Cannot reverse.')
                return;
            end
        end%function
        
        %Redefine a polygon from 
        function RedefineFromList(this,sides,XY_in)
            this.XY = XY_in;
            this.sides = sides;
            this.valid = true;
            this.ComputeNormals;
            if this.graphics_initialized
                this.sketches.Curve.XData = XY_in(:,1);
                this.sketches.Curve.YData = XY_in(:,2);
            end%if
        end%function
        
        %Measure metric properties about this polygon
        function Measure(this)
            
        end%function
        
        %Is point inside?
        function status = IsInside(this,XY_in)
        end%function
                
        %Draw
        function Draw(this,ax)
            if isempty(this.sketch) || ~isvalid(this.sketch)
                this.sketch = custom_line(...
                    'Parent',ax,...
                    'DisplayName','Polygon');
            end%if
            this.sketch.XData = [this.XY(:,1);this.XY(1,1)];
            this.sketch.YData = [this.XY(:,2);this.XY(1,2)];
        end%function
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
        end%function
        
        %Set the axes object onto which to draw the graphics.
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.graphics_initialized
                this.sketches.Curve.Parent = ax;
                this.sketches.Center_Label.Parent = ax;
                this.sketches.Orientation.Parent = ax;
            end%if
        end%function
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',[this.XY(:,1);this.XY(1,1)],... %Do not initalize with empty ("[]" ) because...
                'YData',[this.XY(:,2);this.XY(1,2)],... %MATLAB won't allow ANY property access otherwise.
                'LineStyle','-',...
                'DisplayName',[num2str(this.sides),'-sided polygon']);
         
            %Draw the normal vectors at the midpoint of the edges of the
            %polygon.
            this.sketches.Normals = quiver(...
                this.canvas,...
                zeros(this.sides,1),...
                zeros(this.sides,1),...
                zeros(this.sides,1),... %Normals must be precomputed prior to... 
                zeros(this.sides,1),... %... drawing this quiver object.
                'AutoScale','on',...
                'Color',[0,0,0],...
                'Visible','off',...
                'DisplayName','Polygon Orientation');
            
            if this.normals_computed
                this.sketches.Normals.UData = this.nxy(:,1);
                this.sketches.Normals.VData = this.nxy(:,2);
                %The positions of the arrows are then computed.
                for ii = 1:(this.sides - 1)
                    iip1 = ii + 1;
                    this.sketches.Normals.XData(ii) = (this.XY(iip1,1) + this.XY(ii,1))*0.5;
                    this.sketches.Normals.YData(ii) = (this.XY(iip1,2) + this.XY(ii,2))*0.5;
                end%ii
                this.sketches.Normals.XData(ii) = (this.XY(1,1) + this.XY(this.sides,1))*0.5;
                this.sketches.Normals.YData(ii) = (this.XY(1,2) + this.XY(this.sides,2))*0.5;
            end%if
                   
            %Create the labels for Vertices
            this.sketches.Vertex_Labels = text(...
                this.XY(:,1),...
                this.XY(:,2),...
                '');
            for ii = 1:this.sides
                this.sketches.Vertex_Labels(ii).Visible = 'off';
                this.sketches.Vertex_Labels(ii).Interpreter = 'latex';
                this.sketches.Vertex_Labels(ii).String = ['$V_{',num2str(ii),'}$'];
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
                this.sketches.Segment_Labels(ii).String = ['$S_{',num2str(ii),'}$'];
                this.sketches.Segment_Labels(ii).Position(1) = (this.XY(iip1,1) + this.XY(ii,1))*0.5;
                this.sketches.Segment_Labels(ii).Position(2) = (this.XY(iip1,2) + this.XY(ii,2))*0.5;
            end%ii
            this.sketches.Segment_Labels(this.sides).Visible = 'off';
            this.sketches.Segment_Labels(this.sides).Interpreter = 'latex';
            this.sketches.Segment_Labels(this.sides).String = ['$S_{',num2str(this.sides),'}$'];
            this.sketches.Segment_Labels(this.sides).Position(1) = (this.XY(1,1) + this.XY(this.sides,1))*0.5;
            this.sketches.Segment_Labels(this.sides).Position(2) = (this.XY(1,2) + this.XY(this.sides,2))*0.5;
            
            
            
            this.graphics_initialized  = true;
        end%function.           
        
    end%methods (Graphics)
    
    %High-level functions that CREATE instances of the polynomial objects
    %from low-level code.Error checking involved.
    methods (Static)
        
        %Create a regular polygon
        function poly = CreateRegular(x0,y0,sides,length)
            if sides < 3
                error('Need atleast 3 sides for a polygon.');
            end%if
            poly = polygon;
            poly.sides = sides;
            poly.XY = zeros(sides,2);
            
            %Can this be replace with a "polar loop"?
            theta = 0;
            delta_theta = 2*pi/sides;
            R = length/(2*sin(delta_theta));
            for ii = 1:sides
                poly.XY(ii,1) = x0 + R*cos(theta);
                poly.XY(ii,2) = y0 + R*sin(theta);
                theta = theta + delta_theta;
            end%ii
            poly.convex = true;
            poly.simple = true;
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
            stops = 0;
            stop = []; %IN THE FUTURE, REPLACE THIS WITH THE DOUBLY-LINKED LIST ONCE THAT MATURES.
            for ii = 1:(sides - 1)
                iip1 = ii + 1;
                if XY(ii,1) == XY(iip1,1) && XY(ii,2) == XY(iip1,2)
                    stops = stops + 1;
                    stop = [stop,ii];
                end%if
            end%ii
            
            
            %Initially, take the inputs for granted.
            poly.sides = sides;
            poly.XY = XY;
            
            if XY(1,1) == XY(sides,1)
                 poly.sides = sides - 1;
                 poly.XY = XY(1:end-1,:);
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
            
            poly.valid = true;
            poly.ComputeNormals;
        end%function
        
        %Create a polygon from another polygon's offset.
        function poly = CreateFromOffset(poly2off,thick)
            if ~poly2off.valid
                error('Input polygon is flagged as invalid.')
            end%if
            poly = polygon;
            poly.sides = poly2off.sides;
            poly.XY = polygon.Offset(...
                poly2off.sides,... %Sides of the original polygon.
                poly2off.XY,...%Coordinates of original polygon.
                thick); %Thickness
            poly.valid = true;
            poly.ComputeNormals;
        end%function
        
        %Create a polygon from random coordinates
        function poly = CreateRandom(sides)
            poly = polygon;
            poly.sides = sides;
            poly.XY = rand(sides,2);
            poly.valid = true;
            polt.ComputeNormals;
        end%function
        
    end%methods
    
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
                offset)
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
        
        %Apply a Laplacian smoothing to input XY coordinates
        function XY_out = LaplacianSmooth(sides,XY_in,L)
            
        end%function
        function XY_out = LaplacianSmooth2(sides,XY_in,L,P)
            if nargin < 3
                L = 0.5;
            end%if
            if nargin < 4
                P = 2;
            end%if
            %Use an absolutely disgusting amalgamation of nested for loops
            %and branchless conditionals for to derive the finite
            %difference scheme.
            hp = 0.5*P; %"Half of P."
            M = zeros(P,P);%Allocate memory for coefficient matrix.
            for ww = 1:2
                ii_1 = (ww == 1)*(1) + (ww == 2)*(hp + 1);
                ii_2 = (ww == 1)*(hp) + (ww == 2)*(P);
                dx = (ww == 1)*(-hp) + (ww == 2)*(1);
                for ii = ii_1:ii_2
                    for jj = 1:P
                        dx = dx + 1;
                        M(ii,jj) = power(dx,jj)/factorial(jj); %Optimize this, no need to call factorial.
                    end%jj
                end%ii
            end%ww
            Minv = M^-1;
            coeffs(:,hp + 1) = -sum(Minv,2); %Sum row-wise
            coeffs(:,[1:hp,hp+2:P+1]) = Minv;
            coeffs
            
            %This algorithm assumes that the list of XY coordinates entails
            %a loop.
            XY_out = zeros(sides,2); %Allocate output buffer.
            
            %Points whose "left" neighbors need array access to loop around
            %from back to front.
            for ii = 1:hp
                kk = sides - hp + ii - 1;
                for jj = 1:(P+1)
                    XY_out(ii) = XY_out(ii) + coeffs(jj,2)*XY_in(kk);
                    kk = kk + 1 - sides*(kk == sides);
                end%jj
            end%ii
            
            %Points whose "left" and "right" neighbors can be accessed
            %contiguously.
            for ii = (hp + 1):(sides - hp - 1)
                kk = ii - hp;
                for jj = 1:(P+1)
                    XY_out(ii) = XY_out(ii) + coeffs(jj,2)*XY_in(kk);
                    kk = kk + 1;
                end%jj
                %XY_out(ii)/(P+1)
            end%ii
            
            %Points whose "right" neighbors need array access to loop
            %around from front to back.
            for ii = (sides - hp):sides
                kk = ii;
                for jj = 1:(P+1)
                    XY_out(ii) = XY_out(ii) + coeffs(jj,2)*XY_in(kk);
                    kk = kk + 1 - sides*(kk == sides);
                end%jj
            end%ii
        end%function
        
        %Compute the normal vectors of a polygon's sides. Needed mainly for
        %offseting operations.
        function nxy = Normals2D(sides,XY)
            nxy = zeros(sides,2);%Allocate output memory.
            for ii = 1:(sides - 1)
                iip1 = ii + 1;
                dx = XY(iip1,1) - XY(ii,1);
                dy = XY(iip1,2) - XY(ii,2);
                [nxy(ii,1),nxy(ii,2)] = polygon.UnitNormal2D(dx,dy);
            end%ii
            dx = XY(1,1) - XY(sides,1);
            dy = XY(1,2) - XY(sides,2);
            [nxy(sides,1),nxy(sides,2)] = polygon.UnitNormal2D(dx,dy);
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
%% circle.m
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% A simple data structure used to represent circles in the 2D plane and to
% perform elementary geometric and CAD operations. Basic creation routines
% such as the three-point circle and fillet are available. Point inclusion
% detection.
%% Formulae
% Area
%%
% $A = \pi R^{2}$
%%
% Perimeter
%% 
% $P = 2\pi R$
%% 
% Three-point circle.
%% 
% $\begin{array}{c}
% A_{1} = x_{1} - x_{2}\\
% A_{2} = y_{1} - y_{2}\\
% A_{3} = x_{1} - x_{3}\\
% A_{4} = y_{1} - y_{3}\\
% \|A\| = A_{1}A_{4} - A_{3}A_{2}\\
% B_{1} = x_{1}^{2} - x_{2}^{2} + y_{1}^{2} - y_{2}^{2}\\
% B_{2} = x_{1}^{2} - x_{3}^{2} + y_{1}^{2} - y_{3}^{2}\\
% x_{c} = \frac{A_{4}B_{1} - A_{2}B_{2}}{\|A\|}\\
% y_{c} = \frac{A_{1}B_{2} - A_{3}B_{1}}{\|A\|}\\
% R = \sqrt{(x_{c} - x_{1})^{2} + (y_{c} - y_{1})^{2}}\end{array}$
%%
% Circular fillet (blend two straight lines)
%%
% $\begin{array}{c}
% A_{1} = x_{2} - x_{1}\\
% A_{2} = x_{3} - x_{4}\\
% A_{3} = y_{2} - y_{1}\\
% A_{4} = y_{3} - y_{4}\\
% \|A\| = A_{1}A_{4} - A_{3}A_{2}\\
% L_{1} = \sqrt{A_{1}^{2} + A_{3}^{2}}\\
% L_{2} = \sqrt{A_{2}^{2} + A_{4}^{2}}\\
% B_{1} = x_{3} - x_{1} + R(\frac{A_{3}}{L_{1}} + \frac{A_{4}}{L_{2}})\\
% B_{2} = y_{3} - y_{1} + R(\frac{A_{2}}{L_{2}} + \frac{A_{1}}{L_{1}})\\
% t_{1} = \frac{A_{4}B_{1} - A_{2}B_{2}}{\|A\|}\\
% t_{2} = \frac{A_{1}B_{2} - A_{3}B_{1}}{\|A\|}\\
% x_{c} = x_{1} + t_{1}A_{1} - R\frac{A_{3}}{L_{1}}\\
% y_{c} = y_{1} + t_{1}A_{3} + R\frac{A_{1}}{L_{1}}\\
% \end{array}$
%% Class definition
classdef circle < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the circle's graphics as they appear on axes legend.
        color %Color of circle's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the graphics.
        sketches %Structured array that contains all the plotting handles.
    end%properties (Public)
    %DEFINING DATA variables
    properties (SetAccess = private)
        xc %X-coordinate of center.
        yc %Y-coordinate of center.
        ex %X-component of the assigned orientation.
        ey %Y-component of the assigned orientation.
        R %Radius.
    end%properties (Private)
    %METRIC variables
    properties (SetAccess = protected)
        area
        perimeter
        kappa
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        valid
        
        %Graphics related flags.
        graphics_initialized %DEPRECATE
        canvas_set
        generated
        updated
        refresh
    end%properties (Hidden)
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = circle
            %Geometry-related
            this.xc = [];
            this.yc = [];
            this.ex = 1/sqrt(2);
            this.ey = this.ex;
            this.R = [];
            this.area = [];
            this.perimeter = [];
            this.valid = false;
            this.sketches = struct(...
                'Curve',[],... %Circle is plotted in here.
                'N',50,... %Points on the circle.
                'Center',[],... %Center is plotted. 
                'Radius_Label',[],...
                'Center_Label',[],...
                'Orientation',[]); %A user-specified direction.
            this.generated = struct(...
                'Curve',[],... %Circle is plotted in here.
                'Center',[],... %Center is plotted.
                'Radius_Label',[],...
                'Center_Label',[],...
                'Orientation',[]); %A user-specified direction.
            this.updated = struct(...
                'Curve',[],... %Circle is plotted in here.
                'Center',[],... %Center is plotted.
                'Radius_Label',[],...
                'Center_Label',[],...
                'Orientation',[]); %A user-specified direction.
            this.refresh = struct(...
                'Curve',[],... %Circle is plotted in here.
                'Center',[],... %Center is plotted.
                'Radius_Label',[],...
                'Center_Label',[],...
                'Orientation',[]); %A user-specified direction.

            %Graphics-related.
            this.canvas = [];
            this.color = [0,0,0];
            this.canvas_set = false;
            this.graphics_initialized = false;
        end%function        
        
        %Custom creation routines.
        function this = CreateXYR(xc,yc,R)
            %Define from a center and radius
            this = circle;
            if R < 0
                error('Cannot define a negative radius!');
            end%if
            this.R = R;
            this.xc = xc;
            this.yc = yc;
            this.valid = true;
        end%function
        function this = Create2P(x1,y1,x2,y2)
            %Create a circle from two (2) points. The center is the average
            %of the two points.
            this = circle;
            [this.xc,this.yc] = Center2P(x1,y1,x2,y2);
            this.R = sqrt((x1 - this.xc)^2 + (y1 - this.yc)^2);
        end%function          
        function this = Create3P(x1,y1,x2,y2,x3,y3)
            %Create a circle from three (3) points.
            this = circle;
            [this.xc,this.yc] = circle.Center3P(...
                x1,y1,...
                x2,y2,...
                x3,y3);
            this.R = sqrt((x1 - this.xc)^2 + (y1 - this.yc)^2);
            this.valid = true;
        end%function
        function this = CreateFilletWithRadius(x1,y1,x2,y2,x3,y3,x4,y4,R,ex,ey)
            %Create a circle as a fillet to a 2-segment corner
            if nargin < 9
                L1 = sqrt((x2-x1)^2 + (y2-y1)^2);
                L2 = sqrt((x3-x4)^2 + (y4-y3)^2);
                R = L1*(L1 <= L2) + L2*(L2 < L1);
            end
            if nargin < 10
                ex = 1;
                ey = 0;
            end%if
            this = circle;
            this.ex = ex;
            this.ex = ey;
            this.R = R;
            [this.xc,this.yc,t1,t2,th1,th2] = circle.CenterFilletSegments(...
                x1,y1,...
                x2,y2,...
                x3,y3,...
                x4,y4,...
                R,ex,ey);
            this.valid = true;
        end%function
    end%methods(Static)
    %High-level instance MODIFICATION and QUERY routines.
    methods
        %Compute the metric properties of a circle
        function Measure(this)
            this.perimeter = 2*pi*this.R;
            this.area = 0.5*this.perimeter*this.R;
            this.kappa = 1/this.R;
        end%function
        
        %Redefine from three points
        function Define3P(this,x1,y1,x2,y2,x3,y3)
            old_R = this.R;
            old_xc = this.xc;
            old_yc = this.yc;
            
            [this.xc,this.yc] = circle.Center3P(x1,y1,x2,y2,x3,y3);
            this.R = sqrt((x1 - this.xc)^2 + (y1 - this.yc)^2);
            this.valid = true;
            
            %Graphical update.
            if this.graphics_initialized 
                Show(this);
                UpdateByTranslation(this,this.xc - old_xc,this.yc - old_yc);
                UpdateByScaling(this,this.R/old_R);
            end%if
        end%function        
        function Define2P(this,x1,y1,x2,y2)
            [this.xc,this.yc] = circle.TwoPointCenter(x1,y1,x2,y2);
            this.updated = false;
        end%function
        
        %Affine transformations
        function Displace(this,dx,dy)
            this.xc = this.xc + dx;
            this.yc = this.yc + dy;
            if this.graphics_initialized
                UpdateByTranslation(this,dx,dy);
            end%if
        end%function
        
        %Draw this circle
        function Draw(this,ax)
            this.show_circle = true;
            if ~this.valid
                error('This circle does not have a valid definition!');
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch) || ~isvalid(this.sketch)
                this.sketch = custom_line(...
                    'Parent',ax,...
                    'DisplayName','Circle');
            end%if
            if this.show_center == true && (isempty(this.sketch_center) || ~isvalid(this.sketch_center))
                    this.sketch_center = scatter(ax,...
                        this.xc,...
                        this.yc,...
                        'Visible','off',...
                        'DisplayName','Circle center',...
                        'Marker','o',...
                        'MarkerEdgeColor',[0,0,0],...
                        'MarkerFaceColor',[0,0,1]);
            end%if
            UpdateGraphics(this,this.N_circle);
        end%function
        
        %Rotate the orientation vector included with this circle. Reverse
        %direction of rotation by reverting the sign of the angle.
        function RotateByAngle(this,theta)
            
        end%function
        
        %Rotate the orientation vector included with this circle by an
        %angle corresponding to a certain arc length that has been rolled.
        %Reverse direction of roll by inputting negative arc lengths.
        function RotateByArc(this,arc)
            %d_theta = 2*pi*arc/circle.Perimeter(this.R);
            d_theta = arc/this.R;
            sin_th = sin(d_theta);
            cos_th = cos(d_theta);
            
            old_ex = this.ex;
            old_ey = this.ey;
            
            this.ex = +old_ex*cos_th + old_ey*sin_th;
            this.ey = -old_ex*sin_th + old_ey*cos_th;
            
            if this.graphics_initialized
                this.sketches.Orientation.UData = this.R*this.ex;
                this.sketches.Orientation.VData = this.R*this.ey;
            end%if
        end%function
        
        %Test if this instance of the circle contains an input
        function inside = ContainsPoint(this,xp,yp)
            if ~this.valid
                warning('This circle does not have a valid definition!');
                return;
            end%if
            inside = circle.Inside(this.xc,this.yc,this.R,xp,yp);
        end%function
        
    end%methods (Ordinary)
    %Graphical setups.
    methods
        %Whether this circle object is to be drawn.
        function Show(this)
            %Make sure this object has a definition that is valid.
            if ~this.valid
                warning('The circle object does not have a valid definition!');
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
                this.sketches.Center.Parent = ax;
                this.sketches.Radius_Label.Parent = ax;
                this.sketches.Center_Label.Parent = ax;
                this.sketches.Orientation.Parent = ax;
            end%if
        end%function
        function SetName(this,string)
            this.name = string;
        end%function          
        function SetColor(this,RGB)
            this.color = RGB;
            this.sketches.Curve.Color = RGB;
            this.sketches.Center.MarkerFaceColor = RGB;
            this.sketches.Orientation.Color = RGB;
            this.sketches.Radius_Label.Color = RGB;
            this.sketches.Center_Label.Color = RGB;
        end%function
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',zeros(this.sketches.N,1),... %Do not initalize with empty ("[]" ) because...
                'YData',zeros(this.sketches.N,1),... %MATLAB won't allow ANY property access otherwise.
                'Color',this.color,...
                'LineStyle','-',...
                'DisplayName','Circle');
            theta = 0;
            d_theta = 2*pi/(this.sketches.N - 1);
            for ii = 1:this.sketches.N
                this.sketches.Curve.XData(ii) = this.xc + this.R*cos(theta);
                this.sketches.Curve.YData(ii) = this.yc + this.R*sin(theta);
                theta = theta + d_theta;
            end%ii
            
            %This will sketch the circle's center.
            this.sketches.Center = line(...
                'Parent',this.canvas,...
                'XData',0,... %Do not initalize with empty ("[]" ) because...
                'YData',0,... %MATLAB won't allow ANY property access otherwise.
                'Visible','off',... 
                'Marker','o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',this.color,...
                'DisplayName','Circle center');
           
            %Sketch the orientation.
            this.sketches.Orientation = quiver(...
               this.canvas,...
               this.xc,...
               this.yc,...
               this.ex*this.R,...
               this.ey*this.R,...
               'AutoScale','on',...
               'Color',this.color,...
               'Visible','off',...
               'DisplayName','Circle Orientation');
           
           %The label for the center.
           this.sketches.Center_Label = text(...
               'Parent',this.canvas,...
               'Interpreter','latex',...
               'Color',this.color,...
               'FontName','TimesNewRoman',...
               'Visible','off',...
               'String',['(',num2str(this.xc),';',num2str(this.yc),')'],...
               'Position',[this.xc,this.yc,0]);
           
           %The label for the radius.
           this.sketches.Radius_Label = text(...
               'Parent',this.canvas,...
               'Interpreter','latex',...
               'Color',this.color,...
               'FontName','TimesNewRoman',...
               'Visible','off',...
               'String',['R = ',num2str(this.R)],...
               'Rotation',180*atan(this.ey/this.ex)/pi,...
               'Position',[(this.xc + this.ex*this.R)*0.5,(this.yc + this.ey*this.R)*0.5,0]);
           
            this.graphics_initialized  = true;
        end%function.
                
        %Update graphics based on changes to different properties of the 
        %circle.
        function UpdateRaw(this)
            if ~this.graphics_initialized
                return;
            end%if
            d_theta = 2*pi(this.sketches.N - 1);
            theta = 0;
            for ii = 1:this.sketches.N
                this.sketches.Curve.XData(ii) = this.xc + this.R*cos(theta);
                this.sketches.Curve.YData(ii) = this.yc + this.R*sin(theta);
                theta = theta + d_theta;
            end%ii
            
        end%function
        function UpdateByTranslation(this,dx,dy)
            %THE INPUTS  "dx" AND "dy" ARE ONLY NEEDED BY THE CURVE.
            %This affects ALL the circle's graphics.
            if ~this.graphics_initialized
                return;
            end%if
            %Update the placement of the circle.
            if ~isempty(this.sketches.Curve)
                for ii = 1:length(this.sketches.Curve.XData)
                    this.sketches.Curve.XData(ii) = this.sketches.Curve.XData(ii) + dx;
                    this.sketches.Curve.YData(ii) = this.sketches.Curve.YData(ii) + dy;
                end%
            end%if
            
            %Update the Center's sketch.
            if ~isempty(this.sketches.Center)
                this.sketches.Center.XData = this.xc;
                this.sketches.Center.YData = this.yc;
            end%if
            
            %Update the Orientation's sketch.
            if ~isempty(this.sketches.Orientation)
                this.sketches.Orientation.XData = this.xc;
                this.sketches.Orientation.YData = this.yc;
            end%if
            
            %Update position of the center label.
            if ~isempty(this.sketches.Center_Label)
                this.sketches.Center_Label.Position(1) = this.xc;
                this.sketches.Center_Label.Position(2) = this.yc;
                this.sketches.Center_Label.String = ['(',num2str(this.xc),',',num2str(this.yc),')'];
            end%if
            
            %Update position of the radius label.
            if ~isempty(this.sketches.Radius_Label)
                this.sketches.Radius_Label.Position(1) = (this.xc + this.R*this.ex)*0.5;
                this.sketches.Radius_Label.Position(2) = (this.yc + this.R*this.ey)*0.5;
            end%if
            
        end%function
        function UpdateByRotation(this,d_theta)
        end%function.
        function UpdateByScaling(this,factor)
            %This affects the size of the circle and that of the radius
            %label.
            if ~this.graphics_initialized
                return;
            end%if
            
            %Update the curve.
            for ii = 1:length(this.sketches.Curve.XData)
                this.sketches.Curve.XData(ii) = (this.sketches.Curve.XData(ii) - this.xc)*factor + this.xc;
                this.sketches.Curve.YData(ii) = (this.sketches.Curve.YData(ii) - this.yc)*factor + this.yc;
            end%ii
            
            %Update the Orientation vector.
            this.sketches.Orientation.UData = this.sketches.Orientation.UData*factor;
            this.sketches.Orientation.VData = this.sketches.Orientation.VData*factor;
                        
            %Update position of the radius label.
            this.sketches.Radius_Label.Position(1) = this.xc + this.R*this.ex*0.5;
            this.sketches.Radius_Label.Position(2) = this.yc + this.R*this.ey*0.5;
            this.sketches.Radius_Label.String = ['R = ',num2str(this.R)];
            
        end%function
        function UpdateOrientation(this)
            %This only affects the orientation vector and the radius label.
            if ~this.graphics_initialized
                return;
            end%if
            
            %Update the direction the quiver points to.
            this.sketches.Orientation.XData = this.ex;
            this.sketches.Orientation.YData = this.ey;
            
            %Update the radius label.
            this.sketches.Radius_Label.Rotation = 180*atan(this.ey/this.ex)/pi;
            this.sketches.Radius_Label.Position(1) = (this.xc + this.ex*this.R)*0.5;
            this.sketches.Radius_Label.Position(1) = (this.yc + this.ey*this.R)*0.5;
            
        end%function
             
   
     
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
        
        %Demonstrate computation of a circular fillet.
        function [ax,circ] = FilletDemo
            circle.CleanSlate;

            ax = custom_axis;
            axis(ax,'equal');
            
            fillet_radius = 1;
            P1 = [0,0];
            P2 = [2,2];
            P3 = [-1,3];
            P4 = [1,-1];

            dir_x = cos(rand*2*pi);
            dir_y = sin(rand*2*pi);

            [xi,yi,t1,t2] = circle.Intersect2p2p(...
                P1(1),P1(2),...
                P2(1),P2(2),...
                P3(1),P3(2),...
                P4(1),P4(2));

            quiver(ax,...
                xi,...
                yi,...
                dir_x,...
                dir_y,...
                'Color',[0,0,0],...
                'linewidth',1);



            %Compute the normal for the first line.
            dx1 = P2(1) - P1(1);
            dy1 = P2(2) - P1(2);
            mag1 = sqrt(dx1*dx1 + dy1*dy1);
            nx1 = -dy1/mag1;
            ny1 = dx1/mag1;
            quiver(ax,...
                0.5*(P2(1) + P1(1)),...
                0.5*(P2(2) + P1(2)),...
                nx1,...
                ny1,...
                'Color',[1,0,0],...
                'linewidth',1);

            %Compute the normal for the first line.
            dx2 = P4(1) - P3(1);
            dy2 = P4(2) - P3(2);
            mag2 = sqrt(dx2*dx2 + dy2*dy2);
            nx2 = -dy2/mag2;
            ny2 = dx2/mag2;
            quiver(ax,...
                0.5*(P4(1) + P3(1)),...
                0.5*(P4(2) + P3(2)),...
                nx2,...
                ny2,...
                'Color',[0,0,1],...
                'linewidth',1);
            

            segment1 = plot(ax,[P1(1),P2(1)],[P1(2),P2(2)],'r','linewidth',1);
            segment2 = plot(ax,[P3(1),P4(1)],[P3(2),P4(2)],'b','linewidth',1);
            
            text(ax,P1(1),P1(2),'P$_{1}$','interpreter','latex','Color',[1,0,0])
            text(ax,P2(1),P2(2),'P$_{2}$','interpreter','latex','Color',[1,0,0])
            text(ax,P3(1),P3(2),'P$_{3}$','interpreter','latex','Color',[0,0,1])
            text(ax,P4(1),P4(2),'P$_{4}$','interpreter','latex','Color',[0,0,1])

            I = line('Parent',ax,...
                'XData',xi,...
                'YData',yi,...
                'Marker','s',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[0,0,0]);

            %{
            if dir_x*nx1 + dir_y*ny1 < 0
                [P1(1),P2(1)] = circle.swap(P1(1),P2(1));
                [P1(2),P2(2)] = circle.swap(P1(2),P2(2));
            end%if
            if dir_x*nx2 + dir_y*ny2 < 0
                [P3(1),P4(1)] = circle.swap(P3(1),P4(1));
                [P3(2),P4(2)] = circle.swap(P3(2),P4(2));
            end%if
                %}
            
            circ = circle.CreateFilletWithRadius(...
                P1(1),P1(2),...
                P2(1),P2(2),...
                P3(1),P3(2),...
                P4(1),P4(2),...
                fillet_radius,dir_x,dir_y);

            circ.SetCanvas(ax);
            circ.SetColor([0,0,1]);
            circ.Show;
            
            lgd = legend(ax,...
                [I,segment1,segment2],...
                {'Intersection','Segment 1','Segment 2'},...
                'location','best',...
                'interpreter','latex');
            title(lgd,'LEGEND','interpreter','latex');
            
        end%function
        
        %Demonstrate Circle Rolling
        function RollDemo
            ax = custom_axis;
            C = circle.CreateXYR(0,0,1);
        end%function
        
        %Demonstrate Shaolin Splitting
        function ShaolinDemo(pieces,ex,ey)
            if nargin < 1
                pieces = 2;
            end%if
            if nargin < 2
                ex = 1;
            end%if
            if nargin < 3
                ey = 1;
            end%if
            ax = custom_axis;
            C = circle.CreateXYR(0,0,1);
            C.SetCanvas(ax);
            C.Show;
            
            parts = circle.ShaolinSplit2(C.xc,C.yc,C.R,pieces,50,ex,ey);
            for ii = 2:(pieces-1)
                plot(ax,parts{ii}(:,1),parts{ii}(:,2),'b');
            end%ii
            for ii = [1,pieces]
                plot(ax,parts{ii}(:,1),parts{ii}(:,2),'r');
            end%ii
            scatters = cell(pieces,1);
            for ii = 1:pieces
                scatters{ii} = scatter(ax,parts{ii}(1,1),parts{ii}(1,2),'b','filled','MarkerEdgeColor',[0,0,0]);
            end
            scatters{1}.MarkerFaceColor = [1,0,0];
            scatters{pieces}.MarkerFaceColor = [1,0,0];
            for ii = 1:length(parts{1})
                for jj = 1:pieces
                    set(scatters{jj},'XData',parts{jj}(ii,1));
                    set(scatters{jj},'YData',parts{jj}(ii,2));
                    %input('Next?');
                end%jj
                pause(0.005)
                
            end%jj
            
        end%function
    end%methods
    
 
    %Low-level functions with no error checking specific to this class.
    methods (Static)
        
        %Is a point inside the circle.
        function inside = Inside(xc,yc,R,xp,yp)
            dx = xc - xp;
            dy = yc - yp;
            inside = (R*R > (dx*dx + dy*dy));
        end%function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Make a New class for this.
        %Arc length of a circular sector.
        function arc = SectorArc(R,theta)
            arc = R*theta;
        end%function
        %Perimeter of a circular sector
        function Ps = SectorPerimeter(R,theta)
            Ps = SectorArc(R,theta)+ 2*R;
        end%function
        %Area of a circular sector
        function As = SectorArea(R,theta)
            As = Area(R);
            As = As*theta/(2*pi);
        end%function
        %Chord of a circular sector
        function c = SectorChord(R,theta)
            c = 2*R*sin*(2*theta);
        end%function
        %Sagitta of a circular sector
        function h = SectorSagitta(R,theta)
            c = SectorChord(R,theta);
            h = R - sqrt(R*R - c*c*0.25);
            clear c;
        end%function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Make a New class for this.

        
        %Generate evenly spaced discrete points on the circle.
        function xy_out = points_on(xc,yc,R,N)
            if N < 1
                error('Input must be at least 1!');
            end%if
            xy_out = polar_loop(xc,yc,R,0,2*pi,N);
        end%function
        
        %Generate discrete points on a circular sector.
        function xy_out = points_on_sector(xc,yc,R,theta1,theta2,N)
            if N < 1
                error('N must be greater than 0!');
            end%if
            if theta2 < theta1
                [theta1,theta2] = circle.swap(theta1,theta2);
            end%if
            xy_out = polar_loop(xc,yc,R,theta1,theta2,N);
        end%function
        
        %Produce an array of polar coordinates. The radius is constant
        %(relative to a Cartesian center) and the angle varies evenly 
        %between two values.
        function xy_out = polar_loop(xc,yc,R,theta1,theta2,N)
            xy_out = zeros(N,2);
            theta = theta1;
            delta_theta = (theta2 - theta1)/(N - 1);
            for ii = 1:N
                xy_out(ii,1) = xc + R*cos(theta);
                xy_out(ii,2) = yc + R*sin(theta);
                theta = theta + delta_theta;
            end%ii
        end%function
        
        %Solve for circle's center from three (3) points
        function [xc,yc] = Center3P(x1,y1,x2,y2,x3,y3)
            A1 = x1 - x2;
            A2 = y1 - y2;
            A3 = x1 - x3;
            A4 = y1 - y3;
            detA = A1*A4 - A3*A2;
            B1 = x1^2 + y1^2 - x2^2 - y2^2;
            B2 = x1^2 + y1^2 - x3^2 - y3^2;
            xc = 0.5*(A4*B1 - A2*B2)/detA;
            yc = 0.5*(A1*B2 - A3*B1)/detA;
            %{
            A1 = 2*(x1 - x2); %Maybe there is a way to not have to multiply by 2?
            A2 = 2*(y1 - y2);
            A3 = 2*(x1 - x3);
            A4 = 2*(y1 - y3);
            detA = A1*A4 - A3*A2;
            B1 = x1^2 + y1^2 - x2^2 - y2^2;
            B2 = x1^2 + y1^2 - x3^2 - y3^2;
            xc = (A4*B1 - A2*B2)/detA;
            yc = (A1*B2 - A3*B1)/detA;
            %}
        end%function
        
        %Solve for a circle's center if part of fillet
        function [xc,yc,t1,t2,th1,th2] = CenterFilletSegments(x1,y1,x2,y2,x3,y3,x4,y4,R,ex,ey)
            %Computes the coordinates of the circle from which a fillet
            %between two segments may be constructed. Whether the fillet is
            %possible is determined from the result of finding "t1" and
            %"t2". If any of those two is less than 0 or greater than 1,
            %the fillet is not possible for the given radius.
            
            %If the segments cross, there are four (4) possible fillets
            %that can be computed. To specify which of the four the user
            %needs, prescribe a direction in the form of an x and a y
            %component. The direction must point towards the quadrant where
            %the fillet is needed. Assume that the direction emanates from
            %the intersection location.
            
            A1 = x2 - x1; %+Delta X (Segment #1)
            A2 = x3 - x4; %-Delta X (Segment #2)
            A3 = y2 - y1; %+Delta Y (Segment #1)
            A4 = y3 - y4; %-Delta Y (Segment #2)           
            
            L1 = sqrt(A1*A1 + A3*A3); %Length of Segment #1.
            L2 = sqrt(A2*A2 + A4*A4); %Length of Segment #2.
            
            nx1 = -A3/L1;
            ny1 = +A1/L1;
            nx2 = +A4/L2;
            ny2 = -A2/L2;
            
            if nx1*ex + ny1*ey < 0
                A1 = -A1;
                A3 = -A3;
            end%if
            if nx2*ex + ny2*ey < 0
                A2 = -A2;
                A4 = -A4;
            end%if
            
            B1 = x3 - x1 + R*(A3/L1 + A4/L2);
            B2 = y3 - y1 - R*(A2/L2 + A1/L1);
            detA = A1*A4 - A3*A2;
            t1 = (A4*B1 - A2*B2)/detA;
            t2 = (A1*B2 - A3*B1)/detA;
            
            %Compute the center of the blending circle.
            xc = x1 + t1*A1 - R*A3/L1;
            yc = y1 + t1*A3 + R*A1/L1;
            
            %Blend point on the first segment.
            xb1 = x1 + t1*(x2 - x1);
            yb1 = y1 + t1*(y2 - y1);
            
            %Blend point on the second segment.
            xb2 = x2 + t2*(x4 - x3);
            yb2 = y2 + t2*(y4 - y3); 
            
            %Angular limits for fillet generation. This a four-quadrant
            %inverse tangent operation.
            dxb1 = xb1 - xc;
            dyb1 = yb1 - yc;
            dxb2 = xb2 - xc;
            dyb2 = yb2 - yc;
            th1 = atan(dxb1/dyb1) + pi*(dxb1 < 0) + pi*(dyb1 < 0);
            th2 = atan(dxb2/dyb2) + pi*(dxb2 < 0) + pi*(dyb1 < 0);
        end%function
        function [xc,yc,t1,t2] = CenterFilletLines(x1,y1,dx1,dy1,x2,y2,dx2,dy2,R)
            A1 = x1;
        end%function
        
        
        %Centers circle from midpoint.
        function [xc,yc] = Center2P(x1,y1,x2,y2)
            xc = 0.5*(x1 + x2);
            yc = 0.5*(y1 + y2);
        end%function
        
        
        %Split a circle into evenly "area'd" pieces that resembled the
        %Shaolin symbol.
        function parts = ShaolinSplit2(xc,yc,R,pieces,N,ex,ey)
            if nargin < 6
                ex = 1;
            end%if
            if nargin < 7
                ey = 0;
            end%if
            
            %Allocate memory for the output pieces.
            parts = cell(pieces,1); %First piece has three arcs.
            
            %First split.
            parts{1} = circle.ShaolinTear(xc,yc,R,N,pieces,1,ex,ey);
            
            %Intermediate pieces if split is of three pieces or more.
            if pieces > 2
                for ii = 2:(pieces - 1)
                    parts{ii} = circle.ShaolinShard(xc,yc,R,N,pieces,ii,ex,ey);
                end%ii
            end%if
            
            %Last split
            parts{pieces} = circle.ShaolinTear(xc,yc,R,N,pieces,2,ex,ey);
        end%function
        function XY_out = ShaolinTear(xc,yc,R,N,pieces,tear,ex,ey)
            %THIS CODE COULD BE IMPROVED.
            XY_out = zeros(4*N,2);
            if nargin < 7
                ex = 1;
            end%if
            if nargin < 8
                ey = 0;
            end%if
            mag = sqrt(ex*ex + ey*ey);
            ex = ex/mag;
            ey = ey/mag;
            
            
            
            points = round(4*N/3);
            d_theta = pi/(points - 1);
            if tear == 1
                Ri = R/pieces;
                xci = xc - ex*R + ex*Ri;
                yci = yc - ey*R + ey*Ri;
                delta = + d_theta;
                theta = atan(ey/ex) + 0;
            else
                Ri = R/pieces;
                xci = xc + ex*R - ex*Ri;
                yci = yc + ey*R - ey*Ri;
                delta = +d_theta;
                theta = atan(ey/ex) + pi;
            end
            
            %Generate one arc CCW at the top..
            for ii = 1:points
                XY_out(ii,1) = xci + Ri*cos(theta);
                XY_out(ii,2) = yci + Ri*sin(theta);
                theta = theta + delta;
            end%ii
            %
            
            d_theta = pi/(points - 1);
            if tear == 1
                xci = xc;
                yci = yc;
                Ri = R;
                theta = atan(ey/ex) + pi;
                delta = +d_theta;
            else
                Ri = R;
                xci = xc;
                yci = yc;
                theta = atan(ey/ex) + 0;
                delta = +d_theta;
            end
            
            for ii = (1 + points):(2*points)
                XY_out(ii,1) = xci + Ri*cos(theta);
                XY_out(ii,2) = yci + Ri*sin(theta);
                theta = theta +delta;
            end%ii
            
            Ri = R - R/pieces;
            d_theta = pi/(4*N - 2*points - 1);
            if tear == 1
                xci = xc + ex*R - ex*Ri;
                yci = yc + ey*R - ey*Ri;
                theta = atan(ey/ex) + 0;
                delta = -d_theta;
            else
                xci = xc - ex*R + ex*Ri;
                yci = yc - ey*R + ey*Ri;
                theta = atan(ey/ex) + pi;
                delta = -d_theta;
            end
            
            for ii = (1 + 2*points):4*N
                XY_out(ii,1) = xci + Ri*cos(theta);
                XY_out(ii,2) = yci + Ri*sin(theta);
                theta = theta +delta;
            end%ii
        end%function
        function XY_out = ShaolinShard(xc,yc,R,N,pieces,shard,ex,ey)
            XY_out = zeros(4*N,2); %Allocate output.
            if pieces < 3
                warning('Cannot compute a Shaolin shard when the split is in 2 pieces.');
                return;
            end%if
            if shard == 1 || shard == pieces 
                warning('The Shard # cannot be "1" or the last piece. Default to shard #2.');
                shard = 2;
            end%if
            if nargin < 7
                ex = 1;
            end%if
            if nargin < 8
                ey = 0;
            end%if
            mag = sqrt(ex*ex + ey*ey);
            ex = ex/mag;
            ey = ey/mag;
            
            for jj = 1:4 %All four arcs that make the shard.
                d_theta = (jj == 1)*pi/(N-1) + (jj > 1)*pi/N;
                %Branchless version of
                delta = d_theta*(jj == 1 || jj == 3) - d_theta*(jj == 2 || jj == 4);
                
                %Branchless version of what is bracketed below.
                Ri = (shard*(jj <= 2) + (pieces - shard)*(jj > 2) - (jj == 2) + (jj == 3))*R/pieces;
                %{
                    Ri = ww*R/pieces; %jj= 1
                    Ri = (ww - 1)*R/pieces; %jj = 2
                    Ri = (pieces - ww + 1)*R/pieces; %jj = 3
                    Ri = (pieces - ww)*R/pieces; % jj = 4;
                %}
                
                %Branchless version of what is bracketed below.
                xci = xc + ex*((Ri - R)*(jj <= 2) + (shard*(jj >= 3) - 1*(jj == 3))*R/pieces);
                yci = yc + ey*((Ri - R)*(jj <= 2) + (shard*(jj >= 3) - 1*(jj == 3))*R/pieces);
                %{
                    xci = xc - R + Ri; %jj= 1
                    xci = xc - R + Ri; %jj= 2
                    xci = xc + (ww-1)*R/pieces; %jj= 3
                    xci = xc + (ww)*R/pieces; %jj= 4
                %}
                
                
                %Branchless version of what is bracketed below.
                theta = atan(ey/ex) + 0 + pi*(jj == 2 || jj ==3) - d_theta*(jj == 2 || jj == 4) + d_theta*(jj == 3);
                %{
                    theta = 0; %jj = 1
                    theta = pi - d_theta; %jj = 2
                    theta = pi + d_theta; %jj = 3
                    theta = 0 - d_theta; %%jj = 4
                %}
                %General updater.
                for ii = (1 + (jj - 1)*N):(jj*N)
                    XY_out(ii,1) = xci + Ri*cos(theta);
                    XY_out(ii,2) = yci + Ri*sin(theta);
                    theta = theta + delta;
                end%ii
            end%jj
            
        end%function
        
    end%methods (Static)
    %Elementary low-level functions that are not unique in application to
    %the class.
    methods (Static)
        %"Quality of Life" routines.
        function CleanSlate
            clc
            clear
            close all
        end%function
        function [A,B] = swap(A,B)
            buffer = B;
            B = A;
            A = buffer;
            clear buffer;
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
    end%methods(Public)
end%classef
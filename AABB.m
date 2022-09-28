classdef AABB < handle
    properties (SetAccess = public)
        name %Name of this AABB's graphics as they will appear on axes legends.
        color %Color of AABB's graphics as will be rendered in an axes object.
        canvas %Handle to the axes on which to render the AABB.
        sketches %Structure with the handles to the AABB's graphics.
    end%properties (Public)
    properties (SetAccess = private)
        dim %Dimension of the space this AABB exists in.
        p1 %Coordinates of the minimum point (min_x, min_y,min_z,...).
        p2 %Coordinates of the maximum point (max_x, max_y,max_z,...).
    end%properties(Private)
    properties (SetAccess = protected)
        center
        area
        sizes
        perimeter
    end%properties (protected)
    properties (Hidden = true)
        canvas_set %Whether this AABB was an axes set for rendering.
        graphics_initialized %Whether memory for the graphics has been allocated.
        
        degenerate %Flag AABB that share the same coordinate in p1 and p2.
        valid %Flag to see if the basic operations can be performed on the object.
    end%properties(Protected)
    %High-level functions that MODIFY specific instances of the class
    %object.
    methods
        %Constructor
        function this = AABB(varargin)
            this.name = 'AABB';
            this.color = [0,0,0];
            this.graphics_initialized = false;
            this.canvas_set = false;
            this.sketches = struct(...
                'Box',[],...
                'Center',[]);
            
            this.valid = false;
            this.dim = [];
            this.p1 = [];
            this.p2 = [];
            this.sizes = [];
            %this.sketch = [];
            
            %This part cannot be emulated in C because there is no way to
            %mix the datatypes of the inputs. Everything below this line is
            %a convenience available in MATLAB. Maybe replicable in Python?
            if mod(nargin,2) > 0
                error('Incomplete argument pair!');
            end%if
            ii = 0;
            valid_input = false;
            while ii < nargin
                ii = ii + 1;
                iip1 = ii + 1;
                if strcmp(varargin{ii},'dim') == 1
                    this.dim = varargin{iip1};
                    valid_input = true;
                end%if
                if strcmp(varargin{ii},'p1') == 1
                    this.p1 = varargin{iip1};
                    valid_input = true;
                end%if
                if strcmp(varargin{ii},'p2') == 1
                    this.p2 = varargin{iip1};
                    valid_input = true;
                end%if
                if ~valid_input
                    error('One of the inputs is unrecognizable!');
                end%if
                ii = iip1;
                valid_input = false;
            end%while
        end%function
        
        %Destructor
        %{
        function delete(this)
            if ~isempty(this.p1)
                this.p1 = [];
            end%if
            if ~isempty(this.p2)
                this.p2 = [];
            end%if
            if ~isempty(this.sketches.Box)
                delete(this.sketches.Box);
            end%if
            if ~isempty(this.sketches.Center)
                delete(this.sketches.Center);
            end%if
            
            clear this.dim
            clear this
        end%function
        %}
        
        %Authenticator: (Corrects as much as it can before having to
        %error).
        function authenticate(this)
            if isempty(this.dim)
                error('This AABB does not have a dimension set!');
            end%if
            if isempty(this.p1)
                error('This AABB has an empty p1 list!');
            end%if
            if isempty(this.p2)
                error('This AABB has an empty p2 list!');
            end%if
            
            %Make sure that p1 always has smaller components value than p2.
            %Also, test for degeneracy.
            count = 0;
            for ii = 1:this.dim
                if this.p2(ii) < this.p1(ii)
                    [this.p1(ii),this.p2(ii)] = AABB.swap(this.p1(ii),this.p2(ii));
                end%if
                count = count + 1*(this.p2(ii) == this.p1(ii));
            end%ii
            this.degenerate = ~(count == 0);
            if this.degenerate
                warning('AABB is degenerate! Certain operations are unavailable!');
            end%if
            this.valid = true;
        end%function
        
        %Single factor scaling about the origin.
        function Scale(this,factor)
            for ii = 1:this.dim
                C = 0.5*this.sizes(ii);
                this.p1(ii) = (this.p1(ii) - C)*factor + C;
                this.p2(ii) = (this.p2(ii) - C)*factor + C;
            end%ii
            if this.graphics_initialized
                UpdateByScaling(this,factor);
            end%if
        end%function
        
        %Move in N-dimensional space.
        function Displace(this,displacement)
            for ii = 1:this.dim
                this.p1(ii) = this.p1(ii) + displacement(ii);
                this.p2(ii) = this.p2(ii) + displacement(ii);
            end%ii
            if this.graphics_initialized
                UpdateByTranslation(this,displacement);
            end%if
        end%function
        
        %Measure properties of the AABB
        function Measure(this)
            if ~this.valid
                error('Cannot measure an invalid AABB!');
            end%if
            if isempty(this.sizes) %If memory for this was not allocated in prior call
                this.sizes = zeros(1,this.dim);
            end%if
            this.area = 1; 
            this.perimeter = 0;
            for ii = 1:this.dim
                this.center(ii) = 0.5*(this.p2(ii) + this.p1(ii));
                this.sizes(ii) = this.p2(ii) - this.p1(ii); %This generalizes to higher dimensions.
                this.area = this.area*this.sizes(ii); %This generalizes to higher dimensions.
                this.perimeter = this.perimeter + 2*this.sizes(ii); %This is hardcoded for 2D.
            end%for
        end%function
        
        %Check if an AABB is inside another AABB
        function inside = IsInside(this,AABB2)
            if ~this.valid
                error('AABB#1 is not validated!');
            end%if
            if ~AABB2.valid
                error('AABB#2 is not validated!');
            end%if
            if this.dim ~= AABB2.dim
                error('AABBs are of different sizes.');
            end%if
            inside = AABB.AABBinsideAABB(...
                this.dim,...
                this.p1,...
                this.p2,...
                AABB2.p1,...
                AABB2.p2);
        end%function
        
        %Draw the bounding box
        function Draw(this,ax)
            if ~this.valid
                error('Cannot Draw an invalid AABB!');
            end%if
            if nargin ==1
                ax = custom_axis;
            end%if
            if isempty(this.sketch) || ~isvalid(this.sketch)
                this.sketch = custom_line(...
                    'Parent',ax,...
                    'LineStyle','--',...
                    'DisplayName','AABB');
            end%if
            this.sketch.XData = [this.p1(1),this.p2(1),this.p2(1),this.p1(1),this.p1(1)];
            this.sketch.YData = [this.p1(2),this.p1(2),this.p2(2),this.p2(2),this.p1(2)];
        end%function
        
        %Draw the AABB's center
        function DrawCenter(this,ax)
            if ~this.valid
                error('Cannot Draw an invalid AABB!');
            end%if
            if nargin == 1
                ax = custom_axis;
            end%if
            if isempty(this.sketch_center) || ~isvalid(this.sketch_center)
                this.sketch_center = custom_line(...
                    'Parent',ax,...
                    'Marker','s',...
                    'MarkerFaceColor',[0,0,1],...
                    'DisplayName','AABB center');
            end%if
            this.sketch_center.XData = (this.p2(1) + this.p1(1))/2;
            if this.dim >= 2 
                this.sketch_center.YData = (this.p2(2) + this.p1(2))/2;
            end%if
            if this.dim >= 3
                this.sketch_center.YData = (this.p2(3) + this.p1(3))/2;
            end%if
            
        end%function
        
        %Redefine a box by replacing points
        function RedefinePoints(this,dim,p1,p2)
            this.dim = AABB.branchless_max(this.dim,dim);
            this.degenerate = true;
            for ii = 1:this.dim
                this.p1(ii) = AABB.branchless_min(p1(ii),p2(ii));
                this.p2(ii) = AABB.branchless_max(p1(ii),p2(ii));
                this.degenerate = this.degenerate*(this.p1(ii) == this.p2(ii));
            end%ii
        end%function
        
        %Redefine from a list of points.
        function RedefineFromList(this,points,XY)
            if ~this.valid
                error('AABB is in invalid state.');
            end%if
            [this.p2, this.p1] = AABB.GlobalExtremaOfList(this.dim,points,XY);
            
            if this.graphics_initialized
            end%if
            
        end%function
        
    end%methods (Ordinary)
    %High-level functions that CREATE instances of the polynomial objects
    %from low-level code.Error checking involved.
    
    %Graphical setups.
    methods
        %Whether this AABB object is to be drawn.
        function Show(this)
            %Make sure this object has a definition that is valid.
            if ~this.valid
                warning('The AABB object does not have a valid definition!');
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
                this.InitializeGraphics;
            end%if
        end%function
        
        %Set the axes object onto which to draw the graphics.
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
            if this.graphics_initialized
                this.sketches.Box.Parent = ax;
                this.sketches.Center.Parent = ax;
            end%if
        end%function
        
        %Set the name of this AABB as it would appear in axes legends.
        function SetName(this,string)
            this.name = string;
            if this.graphics_initialized
                this.sketches.Box.DisplayName = string;
                this.sketches.Center.DisplayName = [string, 'Center'];
            end%if
        end%function
        
        %Set the color of the AABB's graphics.
        function SetColor(this,RGB)
            this.color = RGB;
            if this.graphics_initialized
                this.sketches.Box.Color = RGB;
                this.sketches.Center.MarkerFaceColor = RGB;
            end%if
        end%function
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            %This will sketch the circle.
            this.sketches.Box = line(...
                'Parent',this.canvas,...
                'XData',[this.p1(1),this.p1(1),this.p2(1),this.p2(1),this.p1(1)],... %Do not initalize with empty ("[]" ) because...
                'YData',[this.p1(2),this.p2(2),this.p2(2),this.p1(2),this.p1(2)],... %MATLAB won't allow ANY property access otherwise.
                'Visible','on',...
                'LineStyle','--',...
                'DisplayName',this.name);
            
            %This will sketch the AABB's center.
            this.sketches.Center = line(...
                'Parent',this.canvas,...
                'XData',0.5*(this.p1(1) + this.p2(1)),... %Do not initalize with empty ("[]" ) because...
                'YData',0.5*(this.p1(2) + this.p2(2)),... %MATLAB won't allow ANY property access otherwise.
                'Visible','off',... 
                'Marker','s',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[0,0,1],...
                'DisplayName',[this.name,' Center']);
            this.graphics_initialized  = true;
        end%function.
        
        %Show or Hide the box from the canvas.
        function ToggleBox(this)
            if ~this.graphics_initialized
                return;
            end%if
            if strcmp(this.sketches.Box.Visible,'on') == 1
                this.sketches.Box.Visible = 'off';
            else
                this.sketches.Box.Visible = 'on';
            end%if
        end%function
        function ToggleCenter(this)
            if ~this.graphics_initialized
                return;
            end%if
            if strcmp(this.sketches.Center.Visible,'on') == 1
                this.sketches.Center.Visible = 'off';
            else
                this.sketches.Center.Visible = 'on';
            end%if
        end%function
        
        %Update graphics due to repositioning of the box.
        function UpdateByTranslation(this,displacement)
            %Update the graphics for the box.
            for ii = 1:length(this.sketches.Box.XData)
                this.sketches.Box.XData(ii) = this.sketches.Box.XData(ii) + displacement(1);
                this.sketches.Box.YData(ii) = this.sketches.Box.YData(ii) + displacement(2);
            end%ii
            %Update the graphics for the center.
            this.sketches.Center.XData = this.sketches.Center.XData + displacement(1);
            this.sketches.Center.YData = this.sketches.Center.YData + displacement(2);
        end%function
        function UpdateByScaling(this,factor)
            for ii = 1:length(this.sketches.Box.XData)
                this.sketches.Box.XData(ii) = (this.sketches.Box.XData(ii) - this.center(1))*factor + this.center(1);
                this.sketches.Box.YData(ii) = (this.sketches.Box.YData(ii) - this.center(2))*factor + this.center(2);
            end%ii
        end%function
        
        
    end%methods(Graphics)
    
    %Graphical demonstrations
    methods (Static)
        %Showcase the AABB vs. AABB overlap and inclusion detection on
        %randomly generated AABB's in 2D.
        function ax = TestOverlapInclusion
            clc;
            clear;
            close all;
            
            
            boxes = 7;%Number of boxes to generate.
            BOXES = cell(boxes,1);
            colors = {...
                [0,0,1],... %1
                [0,1,0],... %2
                [1,0,0],... %3
                [0,1,1],... %4
                [1,0,1],... %5
                [1,1,0],... %6
                [1,1,1]};   %7
            
            ax = custom_axis; %Create a Canvas.
            ax.Color = [0,0,0]; %Make background black for constrast.
            for ii = 1:boxes
                BOXES{ii} = AABB.CreateRandom(2);
                BOXES{ii}.SetName(['Box #',num2str(ii)]);
                BOXES{ii}.SetCanvas(ax);
                BOXES{ii}.Show;
                BOXES{ii}.SetColor(colors{ii});
            end%ii
            lgd = legend(ax,...
                ax.Children,...
                'location','northeastoutside');
            title(lgd,'LEGEND');
            lgd.Color = [1,1,1];
            
            %Do an O(n^2) search on the boxes.
            for ii = 1:boxes
                for jj = ii+1:boxes
                    overlap = AABB.OverlapsAABB(BOXES{ii},BOXES{jj});
                    jj_inside_ii = AABB.ContainsAABB(BOXES{ii},BOXES{jj});
                    ii_inside_jj = AABB.ContainsAABB(BOXES{jj},BOXES{ii});
                    if overlap
                        fprintf('Box # %i overlaps Box # %i\n',ii,jj);
                    end%if
                    if ii_inside_jj
                        fprintf('Box # %i is inside Box # %i\n',ii,jj);
                    end%if
                    if jj_inside_ii
                        fprintf('Box # %i is inside Box # %i\n',jj,ii);
                    end%if
                end%jj
            end%ii
        end%function
    end%methods (Demonstrations)
    
    methods (Static)
        %Create a random AABB from an input dimension
        function box = CreateRandom(dim,max_sizes)
            box = AABB;
            if nargin == 0
                dim = 2; %Default dimension is 2.
            end%if
            if nargin == 1
                %Deliberately avoiding the use of the "ones" function
                %because it is not available in C. "Zeros" is close enough
                %to C's malloc.
                max_sizes = zeros(1,dim);
                for ii = 1:dim
                    max_sizes(ii) = 1;
                end%end
            end%%if
            box.dim = dim;
            box.p1 = zeros(1,dim);
            box.p2 = zeros(1,dim);
            for ii = 1:dim
                box.p1(ii) = rand*max_sizes(ii);
                box.p2(ii) = rand*max_sizes(ii);
                %In the extremely unlikely event that rand returns the same
                %value
                while box.p1(ii) == box.p2(ii)
                    box.p1(ii) = rand*max_sizes(ii);
                    box.p2(ii) = rand*max_sizes(ii);
                    fprintf('Murphy''s Law! \n');
                end%while
            end %ii
            box.authenticate;
        end%function        
        
        %Create a bounding box from a list of points.
        function box = CreateFromList(dim,points,XY)
            box = AABB;
            box.dim = dim;
            [box.p2, box.p1] = AABB.GlobalExtremaOfList(dim,points,XY);
        end%function
        
    end%methods (Static)
    %High-level routines with substantial error checking to catch
    %user-input errors.
    methods (Static)
        %Create from an AABB from two points.
        function this = Create2P(dim,p1,p2)
            this = AABB;
            this.dim = dim;
            this.p1 = zeros(1,dim);
            this.p2 = zeros(1,dim);
            this.degenerate = true; %Assume worst-case.
            for ii = 1:dim
                this.p1(ii) = AABB.branchless_min(p1(ii),p2(ii));
                this.p2(ii) = AABB.branchless_max(p1(ii),p2(ii));
                this.degenerate = this.degenerate*(this.p1(ii) == this.p2(ii));
            end%ii
            this.valid = true;
        end%function
        
        %Check if two AABB's overlap
        function overlap = OverlapsAABB(AABB1,AABB2)
            if ~AABB1.valid
                error('AABB #1 is flagged as invalid!')
            end%if
            if ~AABB2.valid
                error('AABB #2 is flagged as invalid!')
            end%if
            if AABB1.dim ~= AABB2.dim
                error('AABBs are of different dimensions!');
            end%if
            overlap = AABB.BoxesOverlap(...
                AABB1.dim,...
                AABB1.p1,...
                AABB1.p2,...
                AABB2.p1,...
                AABB2.p2);
        end%function
                
        %Check if one AABB is inside another AABB
        function inside = ContainsAABB(AABB1,AABB2)
            if ~AABB1.valid
                error('AABB #1 is flagged as invalid!')
            end%if
            if ~AABB2.valid
                error('AABB #2 is flagged as invalid!')
            end%if
            if AABB1.dim ~= AABB2.dim
                error('AABBs are of different dimensions!');
            end%if
            %Check the definition of this function further below.
            inside = AABB.AABBinsideAABB(...
                AABB1.dim,...
                AABB2.p1,...
                AABB2.p2,...
                AABB1.p1,...
                AABB1.p2);
        end%function
        
        %Check if a point is inside an AABB
        function inside = ContainsPoint(AABB,point)
            inside = true;%"Innocent until proven guilty."
            for ii = 1:(AABB.dim)
                inside = inside*(AABB.p1(ii) < point(ii) & point(ii) < AABB.p2(ii));
            end%ii
        end%function
        
    end%methods(Static)
    %Low-level functions with no error checking specific to this class.
    methods(Static)
        % Check whether the bounding boxes overlap. This method does not
        % check if one AABB is inside the other.
        function overlap = BoxesOverlap(dim,p1,p2,p3,p4)
            %p1: Minimal coordinates of AABB#1.
            %p2: Maximal coordinates of AABB#1.
            %p3: Minimal coordinates of AABB#2.
            %p4: Maximal coordinates of AABB#2.
            overlap = true; %"Innocent until proven guilty"
            for ii = 1:dim %Check ALL dimensions.
                %These check for overlap. If all dimensions overlap, the 
                %AABB'soverlap.
                if p1(ii) <= p3(ii)
                    overlap = overlap*(p2(ii) > p3(ii));
                end%if
                if p1(ii) > p3(ii)
                    overlap = overlap*(p4(ii) > p1(ii));
                end%if
            end%ii
        end%function
        
        %Check if an AABB#1 is inside AABB#2.
        function inside = AABBinsideAABB(dim,p1,p2,p3,p4)
            %p1: Minimal coordinates of AABB#1.
            %p2: Maximal coordinates of AABB#1.
            %p3: Minimal coordinates of AABB#2.
            %p4: Maximal coordinates of AABB#2.
            %This checks if the AABB defined by p1 and p2 is inside the
            %AABB defined by p3 and p4.
            
            inside = true; %"Innocent until proven guilty."
            for ii = 1:dim %For all dimensions
                size1 = p2(ii) - p1(ii);
                size2 = p4(ii) - p3(ii);
                
                %If AABB#2 is smaller than AABB#1 in any dimension, there
                %is no way AABB#1 can fit inside AABB#2.
                if size2 < size1
                    inside = false;
                    return;
                end%if
               
                %If AABB#1 is not fully between AABB#2 in all coordinates,
                %then it is  not inside AABB#2.
                if ~(p3(ii) < p2(ii) && p2(ii) < p4(ii)) || ~(p3(ii) < p1(ii) && p1(ii) < p4(ii))
                    inside = false;
                    return;
                end%if
            end%ii
        end%function
        
        %Given a list of points with many components, find the maximum and
        %minimum values found in said list.
        function [p_max, p_min] = GlobalExtremaOfList(dim,points,XY)
            p_max = zeros(1,dim);
            p_min = zeros(1,dim);
            for ii = 1:dim
                p_max(ii) = XY(1,ii);
                p_min(ii) = XY(1,ii);
                for jj = 2:points
                    if p_max(ii) < XY(jj,ii)
                        p_max(ii) = XY(jj,ii);
                    end%if
                    if p_min(ii) > XY(jj,ii)
                        p_min(ii) = XY(jj,ii);
                    end%if
                end%jj
            end%ii
        end%function
        
        %Simple swapping function.
        function [A,B] = swap(A,B)
            buffer = B;
            B = A;
            A = buffer;
            clear buffer;
        end%function
        
        %Simple branchless max function.
        function max_AB = branchless_max(A,B)
            max_AB = A*(A >= B) + B*(B > A);
        end%function
        
        %Simple branchless in function.
        function min_AB = branchless_min(A,B)
            min_AB = A*(A <= B) + B*(A > B);
        end%function
        
    end%methods(Static)
end%classdef
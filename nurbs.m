%% nurbs.m (Non-Uniform Rational Basis Splines)
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% NURBS are a collection of a sequence of control points, a corresponding
% sequence of weights, and a "many-to-one" corresponding sequence of 
% "knots" that parametrically define a curve in n-dimensional space as a 
% function of a single parameter. The "knots" are the discrete observations
% of the defining parameter, whereas the control points are the discrete
% observations of quantities of interest (typically XYZ coordinates). The
% weight sequence is a localized tuning parameter which attracts or repels
% the curve (to and from respectively) specific control points. Finally,
% NURBS are characterized by a "degree," which is the same degree used to
% characterize polynomials. Parametrically, NURBS are ratios of polynomials
% and thus, the degree must be specified. The polynomials that form the
% NURBS are the so-called Basis Splines.

% NURBS are an industry standard for conveying the description of shapes
% because they are a highly general case of almost all shapes used in
% engineering. NURBS can exactly reproduce conic sections (circles,
% ellipses, parabolas, and hyperbolas), as well as Bezier curve. The
% advantage of NURBS is that with a single formula, all other shapes can be
% reproduced (as oppossed to having custom formulae for each individual
% shape).
%% Formulae
% NURBS defintion
%%
% $C(u) = \frac{\sum_{i=0}^{n}{w_{i}P_{i}N_{i,p}(u)}}{\sum_{i=0}^{n}{w_{j}N_{i,p}(u)}}$
%%
% Curvature (of any parametric curve).
%%
% N-th derivative of a rational function.
%% Class definition
classdef nurbs < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the nurbs' graphics as they appear on axes legend.
        color %Color of nurbs' graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw.
        sketches %Structured array containing handles to the graphical objects.
        refresh
        updated
        %generated
        
    end%properties(public)
    %DEFINING DATA variables
    properties (SetAccess = private)
        ID %An integer tag. %DEPRECATE?
        d %Highest derivative computed.
        p %Degree of the basis splines.
        o %Order of the NURBS
        kn %Definining knot sequence.
        k1_val %Lower limit of the knot span.
        k2_val %Upper limit of the knot span.
        knots %Number of knots in the sequence.
        intervals %Number of distinct knot intervals.
        granules %Number of discrete points to evaluate between intervals.
        closed % Whether the first CP is identical to the last CP.
        w %Defining weight vector.
        weights %Number of weights. %DEPRECATE, THIS IS REDUNDANT WITH "points".
        CP %Defining control points.
        points %Number of points (length of CP and w).
        dim %Dimension of control points.
        g %Discrete granularity parameter.
        txy
        nxy
        
        XYZ%Discretely interpolated values and derivatives.
        valid %Whether the minimum inputs are valid.
        generated %Boolean flag (whether discrete points have been generated for this curve).
    end
    properties (SetAccess = protected)        
        XYZ_plot %handle for plotting C(u). %DEPRECATE
        CP_scatter %Handle for plotting the control points. %DEPRECATE
        vel_quiver %Handle for plotting tangent vectors. %DEPRECATE
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        %State variables
        CP_defined %DEPRECATE
        kn_defined %DEPRECATE
        w_defined %DEPRECATE
        
        %Graphics state variables
        canvas_set %Whether an axis on which to plot has been set.
        graphics_initialized %Whether the data structures for visualization have been created.
        tangents_computed
        normals_computed
        
        %Computational buffers %DEPRECATE with the future polynomial.m and
        %spline.m
        R_buffer %Stores polynomial coefficients of the rational part of the NURBS.
        S_buffer %Stores polynomial coefficients of the B-splines that are active over a knot interval.        
    end%properties
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = nurbs(varargin)
            this.name = 'NURBS';
            this.color = [0,0,0];
            this.valid = false; %Inputs registered but not validated.
            this.tangents_computed = false;
            this.normals_computed = false;
            this.canvas_set = false;
            this.graphics_initialized = false;
            this.sketches = struct(...
                'Curve',[],...
                'Polygon',[],...
                'Tangents',[],...
                'Normals',[],...
                'CP_Labels',[]);
            %{
            this.generated = struct(...
                'Curve',false,...
                'Polygon',false,...
                'Tangents',false,...
                'Normals',false,...
                'CP_Labels',false);
            %}
            this.updated = struct(...
                'Curve',false,...
                'Polygon',false,...
                'Tangents',false,...
                'Normals',false,...
                'CP_Labels',false);
            this.refresh = struct(...
                'Curve',true,...
                'Polygon',false,...
                'Tangents',false,...
                'Normals',false,...
                'CP_Labels',false);
            
            if nargin == 0
                this.XYZ_plot = [];
                this.CP_scatter = [];
                this.vel_quiver = [];
                this.g = 1;
                this.o = [];
                this.d = [];
                this.XYZ = [];
                this.dim = [];
                this.kn = [];
                this.p = [];
                this.w = [];
                this.k1_val = [];
                this.k2_val = [];
                this.R_buffer = []; %WILL DEPRECATE IN THE FUTURE ONCE...
                this.S_buffer = []; %polynomial.m AND spline.m MATURE.
                return;
            end%if
            if mod(nargin,2) ~= 0
                error('Must input key-value pairs!')
            end%if
            for ii = 1:(nargin/2)
                idx  = 2*ii-1;
                switch varargin{idx}
                    case 'p'
                        this.p = varargin{idx + 1};
                    case 'kn'
                        this.knots = mustbeflat(varargin{idx + 1});
                        this.kn = varargin{idx + 1};
                    case 'w'
                        this.weights = mustbeflat(varargin{idx + 1});
                        this.w = varargin{idx + 1};
                    case 'CP'
                        [array,r,c] = stand(varargin{idx + 1},'tall');
                        this.points = r;
                        this.dim = c;
                        this.CP = array;
                        
                        %%%%% NOT SURE IF THIS IS ROBUST
                        dx = 0;
                        for jj = 1:this.dim
                            dx = dx + (this.CP(1,jj) - this.CP(this.points,jj))^2;
                        end%jj
                        if sqrt(dx) < 1e-16
                            this.closed = true;
                        else
                            this.closed = false;
                        end%if
                        %%%%%
                    case 'g'
                        this.g = varargin{idx + 1};
                    case 'd'
                        this.d = varargin{idx + 1};
                    otherwise
                        error('Unrecognizable NURBS parameter.')
                end%switch
            end%ii
        end%function
        
        %Custom creation routines.
        function nrb = CreateFromPolygon(poly,p,g)
            if ~(exist('polygon','file') == 2)
                error('This operation needs the file polygon.m')
            end%if
            if poly.open
                nrb = nurbs(...
                'CP',poly.XY,...
                'p',p,...
                'g',g);
            else %polygon is closed
                nrb = nurbs(...
                'CP',[poly.XY;poly.XY(1,:)],...
                'p',p,...
                'g',g);
            end%if
            nrb.SetColor(poly.color);
            nrb.generate(2);
        end%function.
        function nrb = CreateFromNURBSOffset(progenitor,offset)
            %This routine takes as input an already existing NURBS. It will
            %then offset said NURBS by a user-specified amount and define
            %the offset points as the control points for *this* nurbs.
            %Granularity and degree are copied.
            nrb = nurbs(...
                'CP',progenitor.Offset(offset),...
                'p',progenitor.p,...
                'g',progenitor.g);
            nrb.generate(2);
        end%function
        
    end%methods (Creation)
    %High-level instance MODIFICATION and QUERY routines.
    methods
        %Memory Allocation and reallocation.s
        function Calloc_P(this,points)
            this.CP = zeros(points,this.dim);
            this.w = zeros(points,this.dim);
            this.kn = zeros(points+this.o,this.dim); %MAKE SURE THIS IS CORRECT!!!
        end%function
        function ReCalloc_P(this,new_points)
            
        end%function
        function Calloc_G(this,g)            
            %Sizes the data buffers according to the granularity to the
            %number of points, knots, and granules.
            this.granules = this.intervals*(g + 1) + 1;
            
            %this.xyz = zeros(this.granules,this.dim);
            this.txy = zeros(this.granules,this.dim);
            this.nxy = zeros(this.granules,this.dim);

        end%function
        function ReCalloc_G(this,new_g)
            %Resizes the data buffers according to the granularity parameter.
            new_granules = this.intervals*(new_g + 1) + 1;
            if new_granules == this.granules
                return;
            elseif new_granules > this.granules
            elseif new_granules < this.granules
            end%if
        end%function
        
        %Definition subroutines.
        function SetKnots(this,knots,kn)
            %WIP
        end%function
        function SetWeights(this,weights,w)
            %WIP
            
            for ii = 1:weights
                this.w = w(ii);
            end%ii
            
        end%function
        
        function EditWeights(this,weights,idx,w_vals)
            %Replace the Old weights with the new ones.
            for ii = 1:weights
                if idx(ii) > this.points || idx(ii) < 0
                    warning([...
                        'Invalid index detected, skipping... (input was )',...
                        num2str(idx(ii)),...
                        'NURBS has only',...
                        num2str(this.points),'weights.'])
                    continue;
                end%if
                this.w(idx(ii)) = w_vals(ii);
            end%ii
            %The curve has been changed, all graphics are outdated...
            this.updated.Curve = false;
            this.updated.Polygon = false;
            this.updated.Tangents = false;
            this.updated.Normals = false;
            this.updated.CP_Labels = false;
            
            %Set a refresh routine here.            
            
        end%function
        function EditKnots(this,knots,idx,kn_vals)
            for ii = 1:knots
                if idx(ii) > this.knots || idx(ii) < 0
                    warning([...
                        'Invalid index detected, skipping... (input was )',...
                        num2str(idx(ii)),...
                        'NURBS has only',...
                        num2str(this.points),'knots.'])
                    continue;
                end%if
                this.kn(idx(ii)) = kn_vals(ii);
            end%ii
            %Run a check to make sure the newly edited sequence is ok.
            if ~nurbs.NonDecreasing(this.knots,this.kn)
                error('Knot sequence is no-longer non-decreasing after edit!');
            end%if
            
            %The curve has been changed, all graphics are outdated...
            this.updated.Curve = false;
            this.updated.Polygon = false;
            this.updated.Tangents = false;
            this.updated.Normals = false;
            this.updated.CP_Labels = false;
            
        end%function
        function EditControlPoints(this,points,idx,vals)
            for ii = 1:points
                if idx(ii) > this.points || idx(ii) < 0
                    warning([...
                        'Invalid index detected, skipping... (input was )',...
                        num2str(idx(ii)),...
                        'NURBS has only',...
                        num2str(this.points),'CP.'])
                    continue;
                end%if
                for jj = 1:this.dim
                    this.CP(idx(ii),jj) = vals(ii,jj);
                end%jj
            end%ii
            
            %The curve has been changed, all graphics are outdated...
            this.updated.Curve = false;
            this.updated.Polygon = false;
            this.updated.Tangents = false;
            this.updated.Normals = false;
            this.updated.CP_Labels = false;
            
            %Set a refresh routine here.            
            
        end%function
        function EditGranules(this,new_granules)
            
        end%function
        
        %Input validation prior to NURBS generation
        function validate(this)
            %Check the control points.
            if isempty(this.CP) == 1
                error('This NURBS does not have control points input.');
            end
            %Should also check for a minimum number of control points
            %required..
            
            %Check the weight array
            if isempty(this.w) == 1
                warning('This NURBS does not have weights set (even weighing applied by default ).');
                this.w = ones(1,this.points);
                this.weights = this.points;
            end%if
            if this.weights ~= this.points
                error('The number of weights does not match number of control points!');
            end%if
            
            %Check degree
            % Must also check for integer input.
            if isempty(this.p) == 1
                this.p = 2;
                warning('This NURBS did not have a degree set (Default of 2 applied).');
            end
            if this.p < 1
                error('NURBS must have B-splines of degree 1 or above.');
            end
            this.o = this.p + 1;
            
            %Check knot sequence.
            if isempty(this.kn) == 1
                warning('This NURBS does not have a knot sequence set (Default sequence applied).');
                this.k1_val = 1;
                this.k2_val = 2;
                this.knots = this.points + this.o; %Number of knots in the parametric sequence.
                this.kn = zeros(1,this.knots); %Allocate memory for the knot sequence.
                delta_kn = (this.k2_val - this.k1_val)/(this.knots - 2*this.o + 1); %Default knot spacing (other than at the endpoints).
                for ii = 1:this.o
                    this.kn(ii) = this.k1_val;
                end%ii
                for ii = (this.p + 2):(this.knots - this.o)
                    this.kn(ii) = this.kn(ii - 1) + delta_kn;
                end%ii
                for ii = (this.knots - this.p):this.knots
                    this.kn(ii) = this.k2_val;
                end%ii
                this.intervals = this.knots - 2*this.o + 1;
            else %Knot sequence was provided.
                this.intervals = 0;
                this.k1_val = this.kn(1);
                %Ensure that knot sequence is nondecreasing.
                for ii = 1:(this.knots - 1)
                    if this.kn(ii + 1) < this(ii)
                        error(['Non-decreasing knot pair detected at index',num2str(ii)]);
                    end%if
                    if this.kn(ii + 1) ~= this.kn(ii)
                        this.intervals = this.intervals +1;
                    end%if
                end%ii
                this.k2_val = this.kn(this.knots);
            end%if
            
            %Final checks
            if this.points + this.o ~= this.knots
                error(['"#points+p+1" = ',...
                    num2str(this.knots),...
                    '+',...
                    num2str(this.p),...
                    '+1 = ',...
                    num2str(this.knots + this.p + 1),...
                    ', does not equal #knots! (',...
                    num2str(this.points),')']);
            end%if
            if isempty(this.g) == 1
                this.g = 10;
            end%if
            if isempty(this.d) == 1
                this.d = 0;
            end%if
            
            %Allocate rational and spline buffers.
            if isempty(this.R_buffer)
                this.R_buffer = zeros(1,this.o);
            end%if
            if isempty(this.S_buffer)
                this.S_buffer = cell(1,this.o);
                for ii = 1:this.o
                    this.S_buffer{ii} = zeros(this.o,this.o);
                end%ii
            end%if
            
            this.valid = true;
        end%function
        
        function generate(this,d)
            this.d = d;
            [this.granules,this.XYZ] = this.evaluate(this.g,this.d);
            this.nxy = zeros(this.granules,2);
            this.txy = zeros(this.granules,2);
            this.ComputeTangents;
            this.ComputeNormals;
            
            this.generated = true;
        end%function
        
        %Orientation.
        function ComputeTangents(this)
            if ~this.generated
                warning('NURBS is not generated! Cannot compute unit tangent vectors.')
                return;
            end%if
            if this.d > 0 %First derivative available.
                for ii = 1:this.granules
                    
                    %Compute the magnitude of the 1st derivative.
                    mag = 0;
                    for jj = 1:this.dim
                        mag = mag + this.XYZ{jj}(ii,2)*this.XYZ{jj}(ii,2);
                    end%jj
                    mag = sqrt(mag);
                    
                    %Compute a normalized tangent vector.
                    for jj = 1:this.dim
                        this.txy(ii,jj) = this.XYZ{jj}(ii,2)/mag;
                    end%jj
                end%ii
                this.tangents_computed = true;
                return;
            end%if
        end%function
        function ComputeNormals(this)
            if ~this.generated
                warning('NURBS is not generated! Cannot compute unit normal vectors.')
                return;
            end%if
            this.dim
            this.tangents_computed
            %This depends on whether the NURBS is 2D or 3D.
            switch this.dim
                case 2 %In 2D one can use the cross product (cheaper).
                    if this.tangents_computed
                        for ii = 1:this.granules
                            [this.nxy(ii,1),this.nxy(ii,2)] = nurbs.UnitNormal2D(this.txy(ii,1),this.txy(ii,2));
                        end%ii
                    end%if
                case 3
                    fprint('WIP!\n');
                otherwise
            end%switch
            this.normals_computed = true;
        end%function
        function ReverseNormals(this)
            for ii = 1:this.granules
                for jj = 1:this.dim
                    this.nxy(ii,jj) = -this.nxy(ii,jj);
                end%jj
            end%ii
            if this.graphics_initialized
                this.sketches.Normals.UData = this.nxy(:,1);
                this.sketches.Normals.VData = this.nxy(:,2);
            end%if
            
        end%function
        
        %Evaluate the rational curve evenly between intervals.
        function [particles,Cu] = evaluate(this,g,d)
            if g < 0
                warning('Granular parameter must be greater than 0! Defaulting to g = 0.')
            end%if
            if d < 1
                error('If querying derivatives, must be one or greater.');
            end%if
            if d > this.p
                error(['Queried derivative (d = ',num2str(d),') greater than NURBS degree ( p =',num2str(this.p),').']);
            end%if
            if this.valid == false
                validate(this);%Callout to validation function above
                if this.valid == false
                    error('Atleast one input is invalid!');
                end%if
            end%if
            particles = this.intervals*(1 + g) + 1;
            this.R_buffer = zeros(1,this.o);%Buffer for the rational part.
            Cu = cell(this.dim,1);
            for ii = 1:this.dim
                Cu{ii} = zeros(particles,d+1);
            end%ii
            S = cell(this.o,1);%Allocate memory for remembering B-splines (In C, this will be a *** pointer).
            sp_idx = zeros(1,this.o); %Stores indices into S.
            kn_idx = zeros(1,this.p + 2); %Knot index buffer (Stores ID's of knots to be used during B-spline construction).
            kn_dif = zeros(1,this.o); %Buffer to store knot spaces when building splines.
            k1 = 1; %Index into knot vector for beginning of a knot span.
            k2 = 1; %Index into knot vector for ending of a knot span.
            for ii = 1:this.o
                kn_idx(ii) = ii; %In C, set to zero instead of 1.
                sp_idx(ii) = ii; %Initialize spline list indexer.
                S{ii} = zeros(this.o,this.o); %In C, this is akin to calling malloc to each ** pointer.
            end%ii
            kn_idx(this.p + 2) = this.p + 2; %In C, set to +1 instead of +2.
            % |First "p+1" B-splines|
            % Here, the minimum number of B-splines is generated to produce the
            % rational part of a NURBS. The first B-spline is constructed using the
            % "triangle fan" scheme. The "p" B-splines after the first one are
            % generated using the "sloped fan" scheme. No interpolation takes place
            % yet. The objective of this first section is to populate the "S" buffer
            % with the first "o" splines.
            sp = 1; %This will be an index into "S."
            for ii = 1:this.o %For the remaining multiplicity of the first knot.
                if ii == 1 %Triangle fan is needed for the first spline.
                    [W,~,~] = bsplgen(this.kn(kn_idx),this.p);%"W" contains lower order splines.
                else %Subsequent splines can be generated using the sloped fan.
                    W = bsplregen(W,this.kn(kn_idx)); %This reuses the lower order splines.
                end
                %Copy from work buffer "W{1}" to spline register "S{sp}".
                for jj = 1:this.o %Rows of buffer (piecewise polynomials).
                    for kk = 1:this.o %Columns of buffer (polynomial coefficients).
                        S{sp}(jj,kk) = S{sp}(jj,kk) + W{1}(jj,kk);
                    end%kk
                end%jj
                %plotspline(S{sp},this.kn(kn_idx),ax4,'b'); %REMOVE ONCE DONE.
                sp = sp + 1;%Tally another spline as being generated.
                %Update knot indexer.
                for jj = 1:(this.p + 2)
                    kn_idx(jj) = kn_idx(jj) + 1;
                end%jj
            end%ii
            
            % |Concurrent interpolation and B-spline generation.|
            % At this stage, "S" contains all B-splines needed to build the rational
            % part. Thus, interpolation can begin. However, further B-splines can be
            % generated on the go. Each loop iteration presents the opportunity for
            % interpolating over a knot span, generating the next B-spline, and
            % rebuilding the rational part over the next knot span.
            
            %Without specific knowledge of the knot sequence, one cannot know the
            %number of knot intervals.
            intvl = 0;
            k1 = 0;
            splines = this.o;
            sp = 1; %Reset the spanner variable.
            PI = 1; %Index of point knot span (PI).
            while intvl < this.intervals
                k1 = k1 + 1;
                while this.kn(k1) == this.kn(k1 + 1) && k1 < (this.knots - 1) %Find the last consecutive repetition.
                    k1 = k1 + 1;
                end%while
                k2 = k1 + 1;
                
                if this.kn(k1) ~= this.kn(k2) %A new distinct knot interval has been discovered.
                    intvl = intvl + 1; %Interpolation can happen.
                    
                    %* Select input knots and spacings. *%
                    
                    dk = (this.kn(k2) - this.kn(k1))/(g + 1); %Uniform knot-spacing within interval.
                    
                    %Assemble the rational buffer.
                    this.R_buffer = zeros(1,this.o);%Rational buffer.
                    for jj = 1:this.o %All splines that are a piece of R.
                        row = this.p + 2 - sp_idx(jj); %Relevant row of S{jj} to this interval.
                        p_idx = PI + sp_idx(jj) - 1; %Control point indexer.
                        for kk = 1:(this.p + 1)
                            this.R_buffer(kk) = this.R_buffer(kk)+ S{jj}(row,kk)*this.w(p_idx);
                        end
                    end
                    
                    for jj = 1:this.o %Indexer for S and active point.
                        row = this.p + 2 - sp_idx(jj); %Relevant row of S{jj} to this interval.
                        p_idx = PI + sp_idx(jj) - 1; %Index of weight and CP corresponding to S{jj}.
                        %if p_idx > this.points %Can this "if" be removed?
                        %    p_idx = this.points;
                        %end%if
                        if ~(splines < this.points)
                            lim = g + 2;
                        else
                            lim = g + 1;
                        end%if
                        %for kk = 1:(g + 1) %Local granule indexer.
                        for kk = 1:lim %Local granule indexer.
                            g_idx = (PI - 1)*(1 + g)+ kk; %Global granule indexer.
                            knot = this.kn(k1) + dk*(kk - 1);
                            dF = polydiffR(S{jj}(row,:),this.R_buffer,d,knot);% DISADVANTAGE: Rational part is unnecesarily derived over an over again.
                            for ww = 1:this.dim %Indexes components of interpolation buffers.
                                for xx = 1:(d + 1) %For function and all queried derivatives.
                                    Cu{ww}(g_idx,xx) = Cu{ww}(g_idx,xx) + dF(xx)*this.CP(p_idx,ww)*this.w(p_idx);
                                end%xx
                            end%ww
                        end%kk
                    end%jj
                    PI = PI + 1; %Tally another point as interpolated.
                    
                    if splines < this.points
                        %Regenerate splines and update spline indexer and spanner.
                        W = bsplregen(W,this.kn(kn_idx));
                        splines = splines + 1;
                        if splines <= this.points - this.o
                            %plotspline(W{1},this.kn(kn_idx),ax4,'b'); %REMOVE ONCE DONE.
                        else
                           % plotspline(W{1},this.kn(kn_idx),ax4,'r'); %REMOVE ONCE DONE.
                        end%if
                        
                        %Update knot indexers.
                        for jj = 1:(this.p + 2)
                            kn_idx(jj) = kn_idx(jj) + 1;
                        end%jj
                        
                        %Reset spanner every "p+2" splines.
                        if sp - this.p - 2 == 0 
                            sp = 1;
                        end%if
                        if sp - 1 == 0 %Permute the indexers.
                            sp_idx(sp) = sp_idx(this.p + 1) + 1;
                        else
                            sp_idx(sp) = sp_idx(sp-1) + 1;
                        end%if
                        for jj = 1:this.p+1 %Knock indexers down by 1 to keep them within range.
                            sp_idx(jj) = sp_idx(jj) - 1;
                        end%jj
                        
                        %Overwrite the "forgetable" spline with new one.
                        for jj = 1:this.o %All rows.
                            for kk = 1:this.o %All polynomial coefficients.
                                S{sp}(jj,kk) = W{1}(jj,kk);
                            end%kk
                        end%jj
                        sp = sp + 1; %Update the spanner when a B-spline is made.
                    end%if
                end%if
            end%while    

            %}            
            %{
            Cu{1}(end,1) = this.CP(end,1); %Fix this eventually.
            Cu{2}(end,1) = this.CP(end,2); %Fix this eventually.
            Cu{1}(end,2) = Cu{1}(end-1,2); %Fix this eventually.
            Cu{2}(end,2) = Cu{2}(end-1,2); %Fix this eventually.
            %}
        end%function
        
        %Evaluate the rational curve at discrete knot values.
        function Cu = evaluate2(this,N,kn_in,d)
            if d < 1
                error('If querying derivatives, must be one or greater.');
            end%if
            if d > this.p
                error(['Queried derivative (d = ',num2str(d),') greater than NURBS degree ( p =',num2str(this.p),').']);
            end%if
            if this.valid == false
                validate(this);%Callout to validation function above
                if this.valid == false
                    error('Atleast one input is invalid!');
                end%if
            end%if
            this.R_buffer = zeros(1,this.o);%Buffer for the rational part.
            Cu = cell(this.dim,1);
            for ii = 1:this.dim
                Cu{ii} = zeros(N,d+1);
            end%ii
            S = cell(this.o,1);%Allocate memory for remembering B-splines (In C, this will be a *** pointer).
            sp_idx = zeros(1,this.o); %Stores indices into S.
            kn_idx = zeros(1,this.p + 2); %Knot index buffer (Stores ID's of knots to be used during B-spline construction).
            kn_dif = zeros(1,this.o); %Buffer to store knot spaces when building splines.
            k1 = 1; %Index into knot vector for beginning of a knot span.
            k2 = 1; %Index into knot vector for ending of a knot span.
            for ii = 1:this.o
                kn_idx(ii) = ii; %In C, set to zero instead of 1.
                sp_idx(ii) = ii; %Initialize spline list indexer.
                S{ii} = zeros(this.o,this.o); %In C, this is akin to calling malloc to each ** pointer.
            end%ii
            kn_idx(this.p + 2) = this.p + 2; %In C, set to +1 instead of +2.
            % |First "p+1" B-splines|
            % Here, the minimum number of B-splines is generated to produce the
            % rational part of a NURBS. The first B-spline is constructed using the
            % "triangle fan" scheme. The "p" B-splines after the first one are
            % generated using the "sloped fan" scheme. No interpolation takes place
            % yet. The objective of this first section is to populate the "S" buffer
            % with the first "o" splines.
            sp = 1; %This will be an index into "S."
            for ii = 1:this.o %For the remaining multiplicity of the first knot.
                if ii == 1 %Triangle fan is needed for the first spline.
                    [W,~,~] = bsplgen(this.kn(kn_idx),this.p);%"W" contains lower order splines.
                else %Subsequent splines can be generated using the sloped fan.
                    W = bsplregen(W,this.kn(kn_idx)); %This reuses the lower order splines.
                end
                %Copy from work buffer "W{1}" to spline register "S{sp}".
                for jj = 1:this.o %Rows of buffer (piecewise polynomials).
                    for kk = 1:this.o %Columns of buffer (polynomial coefficients).
                        S{sp}(jj,kk) = S{sp}(jj,kk) + W{1}(jj,kk);
                    end%kk
                end%jj
                %plotspline(S{sp},this.kn(kn_idx),ax4,'b'); %REMOVE ONCE DONE.
                sp = sp + 1;%Tally another spline as being generated.
                %Update knot indexer.
                for jj = 1:(this.p + 2)
                    kn_idx(jj) = kn_idx(jj) + 1;
                end%jj
            end%ii
            
            % |Concurrent interpolation and B-spline generation.|
            % At this stage, "S" contains all B-splines needed to build the rational
            % part. Thus, interpolation can begin. However, further B-splines can be
            % generated on the go. Each loop iteration presents the opportunity for
            % interpolating over a knot span, generating the next B-spline, and
            % rebuilding the rational part over the next knot span.
            
            %Without specific knowledge of the knot sequence, one cannot know the
            %number of knot intervals.
            intvl = 0;
            k1 = 0;
            splines = this.o;
            sp = 1; %Reset the spanner variable.
            PI = 1; %Index of point knot span (PI).
            kn_evaluated = 1;
            while intvl < this.intervals
                k1 = k1 + 1;
                while this.kn(k1) == this.kn(k1 + 1) && k1 < (this.knots - 1) %Find the last consecutive repetition.
                    k1 = k1 + 1;
                end%while
                k2 = k1 + 1;
                
                if this.kn(k1) ~= this.kn(k2) %A new distinct knot interval has been discovered.
                    intvl = intvl + 1; %Interpolation can happen.
                                          
                    %Assemble the rational buffer.
                    this.R_buffer = zeros(1,this.o);%Rational buffer.
                    for jj = 1:this.o %All splines that are a piece of R.
                        row = this.p + 2 - sp_idx(jj); %Relevant row of S{jj} to this interval.
                        p_idx = PI + sp_idx(jj) - 1; %Control point indexer.
                        for kk = 1:this.o
                            this.R_buffer(kk) = this.R_buffer(kk)+ S{jj}(row,kk)*this.w(p_idx);
                        end%kk
                    end%jj
                    
                    
                    while this.kn(k1) <= kn_in(kn_evaluated) && kn_in(kn_evaluated) < this.kn(k2)
                        for jj = 1:this.o %Indexer for S and active point.
                            row = this.p + 2 - sp_idx(jj); %Relevant row of S{jj} to this interval.
                            p_idx = PI + sp_idx(jj) - 1; %Index of weight and CP corresponding to S{jj}.                           
                            dF = polydiffR(S{jj}(row,:),this.R_buffer,d,kn_in(kn_evaluated));% DISADVANTAGE: Rational part is unnecesarily derived over an over again.
                            for ww = 1:this.dim %Indexes components of interpolation buffers.
                                for xx = 1:(d + 1) %For function and all queried derivatives.
                                    Cu{ww}(kn_evaluated,xx) = Cu{ww}(kn_evaluated,xx) + dF(xx)*this.CP(p_idx,ww)*this.w(p_idx);
                                end%xx
                            end%ww
                        end%jj
                        kn_evaluated = kn_evaluated + 1;
                    end%if
                    PI = PI + 1; %Tally another point as activated.
                    
                    if splines < this.points
                        %Regenerate splines and update spline indexer and spanner.
                        W = bsplregen(W,this.kn(kn_idx));
                        splines = splines + 1;
                                                
                        %Update knot indexers.
                        for jj = 1:(this.p + 2)
                            kn_idx(jj) = kn_idx(jj) + 1;
                        end%jj
                        
                        %Reset spanner every "p+2" splines.
                        if sp - this.p - 2 == 0 
                            sp = 1;
                        end%if
                        if sp - 1 == 0 %Permute the indexers.
                            sp_idx(sp) = sp_idx(this.p + 1) + 1;
                        else
                            sp_idx(sp) = sp_idx(sp-1) + 1;
                        end%if
                        for jj = 1:this.p+1 %Knock indexers down by 1 to keep them within range.
                            sp_idx(jj) = sp_idx(jj) - 1;
                        end%jj
                        
                        %Overwrite the "forgetable" spline with new one.
                        for jj = 1:this.o %All rows.
                            for kk = 1:this.o %All polynomial coefficients.
                                S{sp}(jj,kk) = W{1}(jj,kk);
                            end%kk
                        end%jj
                        sp = sp + 1; %Update the spanner when a B-spline is made.
                    end%if
                end%if
                if kn_evaluated >= N
                    break;
                end%if
            end%while    
        end%function
        
        %Add thickness to the NURBS
        function XYU = thicken(this,thickness)
        %function [XYU,XYL] = thicken(this,thickness,sign)
            if this.generated == false
                error('NURBS is not generated. Cannot thicken.');
            end
            if this.d < 1
                error(['Need atleast first derivative for this operation. This NURBS has d = ',num2str(this.d)]);
            end
            %In C, the preallocations are outside the function.
            XYU = zeros(this.granules,this.dim); %Preallocate memory for the "upper".
            %XYL = zeros(this.granules,this.dim); %Preallocate memory for the "lower".
            for ii = 1:this.granules
                vel = [this.XYZ{1}(ii,2),this.XYZ{2}(ii,2)];
                [un,~] = vec2nor(vel);
                for jj = 1:this.dim
                    XYU(ii,jj) = this.XYZ{jj}(ii,1) + un(jj)*thickness;
                    %XYL(ii,jj) = this.XYZ{jj}(ii,1) - un(jj)*thickness;
                end%jj
            end%ii
        end%function
        
        %Create an orthogonal offset.
        function [XY] = Offset(this,distance)
            if this.generated == false
                error('NURBS is not generated. Cannot thicken.');
            end%if
            if this.d < 1
                error(['Need atleast first derivative for this operation. This NURBS has d = ',num2str(this.d)]);
            end%if
            XY = zeros(this.granules,this.dim); %Preallocate memory for the "upper".
            for ii = 1:this.granules
                vel = [this.XYZ{1}(ii,2),this.XYZ{2}(ii,2)];
                [un,~] = vec2nor(vel);
                for jj = 1:this.dim
                    XY(ii,jj) = this.XYZ{jj}(ii,1) + un(jj)*distance;
                end%jj
            end%ii
            clear vel
            clear un            
        end%function
        
        function force_close(this)
            %In C, we will need to use the realloc function.
            this.CP = [this.CP,this.CP(1,:)];
            this.kn = [];
            this.w = [];
            this.valid = false;
            this.validate; %Rebuil weight vectors and 
            if this.generated == true
                this.generate(this.d);
            end%
        end%function
        
        function SaveAsIGES(this)
            
        end%function
        
        % Generate control points by extruding NURBS orthogonally along the
        % curve.
        function CPx = extrapolate(this,thickness)
            if ~this.valid 
                error('NURBS has not been validated.')
            end%if
            if this.d < 1
                error(['Need atleast first derivative for this operation. This NURBS has d = ',num2str(this.d)]);
            end%if
            %
            %[particles,Cu] = this.evaluate(0,1); %Evaluate evenly spaced points on progenitor NURBS.
            
            knot_vals = linspace(this.kn(1),this.kn(this.knots),this.points);
            Cu = this.evaluate2(this.points,knot_vals,1); %Evaluate evenly spaced points on progenitor NURBS.
            
            %Copy first point to last.
            %[Cu{1}(1,1),Cu{2}(1,1)]
            %[Cu{1}(end,1),Cu{2}(end,1)]
            for ii = 1:this.dim %dimensions
                for jj = 1:2 %derivatives - 1
                        Cu{ii}(this.points,jj) = Cu{ii}(1,jj);
                end%jj
            end%ii
            
            CPx = zeros(this.points,this.dim); %Allocate memory for CP of child NURBS.
            vel = zeros(1,this.dim);
            for ii = 1:this.points
                for jj = 1:this.dim
                    vel(jj) = Cu{jj}(ii,2);
                end%for
                [un,~] = vec2nor(vel); %Gather unit normal
                for jj = 1:this.dim
                    CPx(ii,jj) = Cu{jj}(ii,1) + thickness*un(jj);
                end%jj
            end%ii
            clear vel Cu particles
            %

        end%function
        
        %Evaluate NURBS at discrete knot values.
        function XYZ_out = eval(this,kn_in)
            knots_in = mustbeflat(kn_in);
            for ii = 1:knots_in
                if kn_in(ii) < this.k1_val || kn_in(ii) > this.k2_val
                    error(['An input knot with index',...
                        num2str(ii),...
                        'and value of ',...
                        num2str(kn_in(ii)),...
                        ' is outside this NURBS''s knot range (',...
                        'num2str(',...
                        num2str(this.k1_val),...
                        ',',num2str(this.k2_val),...
                        ')']);
                end%if
            end%ii
            
            %Memory Allocation
            XYZ_out = zeros(knots_in,this.dim); %Allocate memory for output.
            knots_evaluated = 0; %Initialize a counter for # of knots evaluated.
            
            S = cell(this.o,1);%Allocate memory for remembering B-splines (In C, this will be a *** pointer).
            sp_idx = zeros(1,this.o); %Stores indices into S.
            
            
            while knots_evaluated < knots_in %Until all knots are evaluated.
                
            end%while
            
        end%function
        
        %Compute the curvature of the NURBS
        function kappa = curvature(this)
            if this.generated == false
                error('Cannot perform this operation if NURBS is not generated!');
            end%if
            if this.d < 2
                error(['Need atleast second derivative for this operation. This NURBS has d = ',num2str(this.d)]);
            end%if
            kappa = zeros(this.granules,1);
            
            %THIS IS HARDCODED FOR 2D, FOR THIS TO GENERALIZE TO 3D, NEED A
            %DEDICATED CROSS PRODUCT FUNCTION.
            for ii = 1:this.granules
                kappa(ii) = this.XYZ{1}(ii,2)*this.XYZ{2}(ii,3) - this.XYZ{2}(ii,2)*this.XYZ{1}(ii,3);
                kappa(ii) = kappa(ii)/power(this.XYZ{1}(ii,2)^2 + this.XYZ{2}(ii,2)^2   ,1.5);
            end%ii
        end%ii
        
        %DEPRECATE
        %Draw the NURBS curve
        function draw(this,ax)
            if nargin == 1 %No axes specified.
                ax = custom_axis;
            end%if
            if this.generated == false
                error('Cannot produce plots if the NURBS is not generated!');
            end%if
            if isempty(this.XYZ_plot) == 1
                this.XYZ_plot = line; %Instantiate the line.
                this.XYZ_plot.DisplayName = 'NURBS';
                this.XYZ_plot.LineWidth = 1; %
                this.XYZ_plot.Parent = ax; %Link this line plot to the axes.
                this.XYZ_plot.Color = [0,0,0]; %Default Color is black.
            end%if
            this.XYZ_plot.XData = this.XYZ{1}(:,1);
            if this.dim > 1
                this.XYZ_plot.YData = this.XYZ{2}(:,1);
            end%if
            if this.dim > 2
                this.XYZ_plot.ZData = this.XYZ{3}(:,1);
            end%if
        end%function
        %DEPRECATE
        
        %DEPRECATE
        %Draw Tangent vectors to NURBS
        function draw_velocity(this,ax)
            if nargin == 1 %No axes specified.
                ax = custom_axes;
            end%if
            if this.d < 1
                error('NURBS needs atleast its first derivative computed for a "velocity" plot!');
            end%if
            if this.generated == false
                error('Cannot produce plots if the NURBS is not generated!');
            end
            this.vel_quiver = quiver(ax,...
                this.XYZ{1}(:,1),... %XData
                this.XYZ{2}(:,1),... %YData
                this.XYZ{1}(:,2),... %UData
                this.XYZ{2}(:,2),... %VData
                'DisplayName','NURBS "Velocity"',...
                'Visible','on',...
                'Color',[1,0,0]...
                );            
        end
        %DEPRECATE
        
        %DEPRECATE
        %Draw the control points that define this NURBS.
        function draw_CP(this,ax)
            if nargin == 1 %No axes specified.
                ax = custom_axes;
            end%if
            if this.generated == false
                error('Cannot produce plots if the NURBS is not generated!');
            end
            this.CP_scatter = scatter(ax,... %Axis handle
                this.CP(:,1),...%XData
                this.CP(:,2),...%YData;
                'Marker','o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[1,0,0]...
                );
            if this.dim > 2
                this.CP_scatter.ZData = this.CP(:,3);
            end
            this.CP_scatter.DisplayName = 'NURBS CP';
        end
        %DEPRECATE
        %DEPRECATE
        %Label controls points
        function txt = label_CP(this,ax)
            %Create the numbered labels.
            labels = cell(1,this.points);
            for ii = 1:this.points
                labels{ii} = num2str(ii);
            end%ii
            txt = text(ax,...
                this.CP(:,1),... %XData
                this.CP(:,2),... %YData
                labels,...
                'Visible','on',... %Visibility toggle.
                'interpreter','latex',...
                'FontName','Helvetica');
        end%function.
        %DEPRECATE
        
        %Plot the Basis Splines
        function draw_bsplines(this,ax)
            if nargin == 1 %No axes specified.
                ax = custom_axes;
            end%if
            if this.generated == false
                error('Cannot produce plots if the NURBS is not generated!');
            end
            
           %Some preallocation
            S = cell(this.o,1);%Allocate memory for remembering B-splines (In C, this will be a *** pointer).
            kn_idx = zeros(1,this.p + 2); %Knot index buffer (Stores ID's of knots to be used during B-spline construction).
            sp_idx = zeros(1,this.o);
            for ii = 1:(this.o)
                kn_idx(ii) = ii;
                sp_idx(ii) = ii; %Initialize spline list indexer.
                S{ii} = zeros(this.o,this.o); %In C, this is akin to calling malloc to each ** pointer.
            end%ii
            kn_idx(this.o + 1) = ii + 1;
            
            sp = 1;
            for ii = 1:this.o %For the remaining multiplicity of the first knot.
                if ii == 1 %Triangle fan is needed for the first spline.
                    [W,~,~] = bsplgen(this.kn(kn_idx),this.p);%"W" contains lower order splines.
                else %Subsequent splines can be generated using the sloped fan.
                    W = bsplregen(W,this.kn(kn_idx)); %This reuses the lower order splines.
                end
                %Copy from work buffer "W{1}" to spline register "S{sp}".
                for jj = 1:this.o %Rows of buffer (piecewise polynomials).
                    for kk = 1:this.o %Columns of buffer (polynomial coefficients).
                        S{sp}(jj,kk) = S{sp}(jj,kk) + W{1}(jj,kk);
                    end
                end
                plotspline(S{sp},this.kn(kn_idx),ax,'b'); %REMOVE ONCE DONE.
                sp = sp + 1;%Tally another spline as being generated.
                %Update knot indexer.
                for jj = 1:(this.p + 2)
                    kn_idx(jj) = kn_idx(jj) + 1;
                end%jj
            end%ii
            PI = 1; %Index of point knot span (PI).
            intvl = 0;
            k1 = 0 ;
            splines = this.o;
            while intvl < this.intervals %Plot splines, until it is done.
                k1 = k1 + 1;
                while this.kn(k1) == this.kn(k1 + 1) && k1 < (this.knots - 1) %Find the last consecutive repetition.
                    k1 = k1 + 1;
                end%while
                k2 = k1 + 1;
                if this.kn(k1) ~= this.kn(k2) %A new distinct knot interval has been discovered.
                    intvl = intvl + 1; %Interpolation can happen.
                    
                    %Assemble the rational buffer.
                    this.R_buffer = zeros(1,this.o);%Rational buffer.
                    for jj = 1:this.o %All splines that are a piece of R.
                        row = this.p + 2 - sp_idx(jj); %Relevant row of S{jj} to this interval.
                        p_idx = PI + sp_idx(jj) - 1; %Control point indexer.
                        for kk = 1:(this.p + 1)
                            this.R_buffer(kk) = this.R_buffer(kk)+ S{jj}(row,kk)*this.w(p_idx);
                        end%kk
                    end%jj
                    PI = PI + 1;
                    
                    plotspline(this.R_buffer,[this.kn(k1),this.kn(k2)],ax,'k'); %REMOVE ONCE DONE.
                    if splines < this.points
                        %Regenerate splines and update spline indexer and spanner.
                        W = bsplregen(W,this.kn(kn_idx));
                        splines = splines + 1;
                        if splines <= this.points - this.o
                            plotspline(W{1},this.kn(kn_idx),ax,'b'); %REMOVE ONCE DONE.
                        else
                            plotspline(W{1},this.kn(kn_idx),ax,'r'); %REMOVE ONCE DONE.
                        end%if
                        %Update knot indexers.
                        for jj = 1:(this.p + 2)
                            kn_idx(jj) = kn_idx(jj) + 1;
                        end%jj
                    end%if
                    if sp - this.p - 2 == 0 %Reset spanner every "p+2" splines.
                        sp = 1;
                    end%if
                    if sp - 1 == 0 %Permute the indexers.
                        sp_idx(sp) = sp_idx(this.p+1) + 1;
                    else
                        sp_idx(sp) = sp_idx(sp-1) + 1;
                    end%if
                    for jj = 1:this.o %Knock indexers down by 1 to keep them within range.
                        sp_idx(jj) = sp_idx(jj) - 1;
                    end%jj
                    
                    %Overwrite the "forgetable" spline with new one.
                    for jj = 1:this.o %All rows.
                        for kk = 1:this.o %All polynomial coefficients.
                            S{sp}(jj,kk) = W{1}(jj,kk);
                        end%kk
                    end%jj
                    sp = sp + 1; %Update the spanner when a B-spline is made.
                end%if
            end%while
            
            
            function plotspline(B,kn,ax,color)
                [rB,cB] = size(B);
                for pp = 1:rB
                    u = linspace(kn(pp),kn(pp+1),10);
                    val = horner(B(pp,:),u);
                    plot(ax,u,val,color,'linewidth',1);
                end%pp
            end%function (nested)            
        end%function
        
        %Replace the NURBS CP
        function CP_new = cusp_sniper(this,kappa_tol,ax)
            if this.d < 2
                error('Need at least second derivative for this operation!');
            end%if
            kappa = curvature(this); %Compute the curvature of the NURBS.
            [L,idx,xtr] = xtrloc(kappa); %Detect the maxima/minima in the curvature.
            [f,ax1] = xtr_viewer2(kappa,idx,xtr,{'$\kappa$'},[]);
            %Make a copy of the discrete points on the NURBS.
            CP_new = [
                this.XYZ{1}(:,1),...
                this.XYZ{2}(:,1)...
                ];
            
            idx1 = zeros(L,1);
            idx2 = zeros(L,1);
            e1 = ellipse; %The blending is accomplished using an ellipse.
            for ii = 1:L
                %Need only look at the maxima
                if xtr(ii) == 1 %If the extremum is a local maximum
                    if abs(kappa(idx(ii))) > kappa_tol %Cruvature too high ?
                        kappa = abs(kappa(idx(ii)));
                        idx1(ii) = idx(ii) - 1;%Index of point before cusp.
                        idx2(ii) = idx(ii) + 1;%Index of point after cusp.
                        %Enforce the curvature constraint but don't go out
                        %of bounds.
                        while kappa > kappa_tol && idx1(ii) > 1 && idx2(ii) < length(this.XYZ{1}) %&& kk <= 3
                            %Fit ellipses
                            p1 = [this.XYZ{1}(idx1(ii),1),this.XYZ{2}(idx1(ii),1)];
                            d1 = [this.XYZ{1}(idx1(ii),2),this.XYZ{2}(idx1(ii),2)];
                            p2 = [this.XYZ{1}(idx2(ii),1),this.XYZ{2}(idx2(ii),1)];
                            d2 = [this.XYZ{1}(idx2(ii),2),this.XYZ{2}(idx2(ii),2)];
                            e1.genfrom2dir2p(p1,d1,p2,d2);
                            kappa = (e1.a^2)/e1.b;%Curvature at the major radius.
                            e1.plot(ax);
                            idx1(ii) = idx1(ii) - 1;
                            idx2(ii) = idx2(ii) + 1;
                            input('Next Blend?');
                        end%while
                        
                        %If while loop terminates with indices within
                        %bounds.
                        if idx1(ii) >= 1 && idx2(ii) <= length(this.XYZ{1})
                            %Implode points onto the ellipse. And overwrite
                            %the copy.
                            CP_new(idx1(ii):idx2(ii),:) = e1.attract(...
                                [this.XYZ{1}(idx1(ii):idx2(ii),1),...
                                 this.XYZ{2}(idx1(ii):idx2(ii),1)]);
                        end %if
                        
                    end%if
                end%if
            end%ii
            plot(ax,CP_new(:,1),CP_new(:,2),'Color',[1,0,0],'linewidth',1);
            clear kappa
        end%function
        
    end%methods (Ordinary)
    %Graphical setup
    methods 
        %Signal that graphics are to be generated.
        function Show(this)
            %Make sure this object has a definition that is valid.
            if ~this.valid
                warning('The NURBS object does not have a valid definition!');
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
                this.sketches.Polygon.Parent = ax;
                this.sketches.Tangents.Parent = ax;
                this.sketches.Normals.Parent = ax;
            end%if
        end%function
        function SetColor(this,RGB)
            this.color = RGB;
            if this.graphics_initialized
                this.sketches.Curve.Color = RGB;
                this.sketches.Polygon.Color = RGB;
                this.sketches.Polygon.MarkerFaceColor = RGB;
                this.sketches.Tangents.Color = RGB;
                this.sketches.Normals.Color = RGB;
                for ii = 1:this.points
                    this.sketches.CP_Labels(ii).Color = RGB;
                end%ii
            end%if
        end%function
        function SetName(this,name)
            this.name = name;
            if this.graphics_initialized
                this.sketches.Curve.DisplayName = name;
                this.sketches.Polygon.DisplayName = [name,'Control Curve'];
                this.sketches.Tangents.DisplayName = [name, 'Unit Tangents'];
                this.sketches.Normals.DisplayName = [name,'Unit Normals'];
            end%if
        end%function
        
        %Create the Graphical objects.
        function InitializeGraphics(this)
            %This function assumes that a canvas has been set already.
            if ~this.generated
                warning('NURBS is not generated, cannot display.');
                return;
            end%if
            
            %This will sketch the rational curve C(u).
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',this.XYZ{1}(:,1),... %Do not initalize with empty ("[]" ) because...
                'YData',this.XYZ{2}(:,1),... %MATLAB won't allow ANY property access otherwise.
                'Color',this.color,...
                'LineStyle','-',...
                'DisplayName',this.name);
            
            %This will sketch NURBS's control polygon.
            this.sketches.Polygon = line(...
                'Parent',this.canvas,...
                'XData',this.CP(:,1),... %Do not initalize with empty ("[]" ) because...
                'YData',this.CP(:,2),... %MATLAB won't allow ANY property access otherwise.
                'Visible','off',... 
                'LineStyle','--',...
                'Marker','o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',this.color,...
                'DisplayName',[this.name,'Control Curve']);
           
            %This will sketch the NURBS's unit tangent vectors.
            this.sketches.Tangents = quiver(...
                this.canvas,...
                this.XYZ{1}(:,1),... %Do not initalize with empty ("[]" ) because...
                this.XYZ{2}(:,1),... %MATLAB won't allow ANY property access otherwise.
                this.txy(:,1),...
                this.txy(:,2),...
                'Color',this.color,...
                'DisplayName',[this.name,'Unit Tangents'],...
                'Visible','off');
            
            %This will sketch the NURBS's unit normal vectors.
            this.sketches.Normals = quiver(...
                this.canvas,...
                this.XYZ{1}(:,1),... %Do not initalize with empty ("[]" ) because...
                this.XYZ{2}(:,1),... %MATLAB won't allow ANY property access otherwise.
                this.nxy(:,1),...
                this.nxy(:,2),...
                'DisplayName',[this.name,'Unit Normals'],...
                'Color',this.color,...
                'Visible','off');
            
            this.graphics_initialized  = true;
        end%function.
        function GenerateCPLabels(this)
            %This will render labels for the control vertices.
            %Create the labels for Vertices
            this.sketches.CP_Labels = text(...
                this.CP(:,1),...
                this.CP(:,2),...
                '');
            for ii = 1:this.points
                this.sketches.CP_Labels(ii).Color = this.color;
                this.sketches.CP_Labels(ii).Visible = 'off';
                this.sketches.CP_Labels(ii).Interpreter = 'latex';
                this.sketches.CP_Labels(ii).String = {...
                    ['P$_{',num2str(ii),'}$'],...
                    ['w = ',num2str(this.w(ii))]...
                    };
            end%ii
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
                elseif strcmpi(varargin{ii},'Polygon') == 1
                    VisibilityToggle(this,this.sketches.Polygon);
                elseif strcmpi(varargin{ii},'Normals') == 1
                    VisibilityToggle(this,this.sketches.Normals);
                elseif strcmpi(varargin{ii},'Tangents') == 1
                    VisibilityToggle(this,this.sketches.Tangents);
                elseif strcmpi(varargin{ii},'CP_Labels') == 1
                    for jj = 1:this.points
                        VisibilityToggle(this,this.sketches.CP_Labels(jj));
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
    end%methods (Graphics)
    %Graphical demonstrations
    methods (Static)
        
        %Draw a NURBS circle.
        function [ax,nrb] = CircleDemo
            clc;
            clear;
            close all;
            
            ax = custom_axis;
            axis(ax,'equal');
            R = 1; %Radius of the circle.
            xc = 0.5;
            yc = 0.5;
            
            %Define the control points.
            control_points = [...
                xc + R, yc;...     %Point 0
                xc + R, yc + R;... %Point 1
                xc - R, yc + R;... %Point 2
                xc - R, yc;...     %Point 3
                xc - R, yc - R;... %Point 4
                xc + R, yc - R;... %Point 5
                xc + R, yc];       %Point 6 (Repeat of Point 0).
            
            %Define the weights.
            control_weights = [...
                1.0,... %Point 0
                0.5,... %Point 1
                0.5,... %Point 2
                1.0,... %Point 3
                0.5,... %Point 4
                0.5,... %Point 5
                1.0];   %Point 6 (Repeat of Point 0).
            
            %NOTE: The knot vector can be arbitrary. Let the nurbs routine
            %define a default.
            nrb = nurbs(...
                'CP',control_points,... %The XY points defined above.
                'w',control_weights,... %The weights defined above.
                'p',2,... %Degree of the piecewise polynomials.
                'g',9); %Granularity parameter for plotting.
            nrb.generate(2);
            nrb.SetCanvas(ax);
            nrb.Show;
        end%function
        
    end%methods (Graphical demonstrations)
    methods (Static)
        function [nx,ny] = UnitNormal2D(dx,dy)
            %Compute a unit normal.
            mag = sqrt(dx*dx + dy*dy); %Magnitude of the input direction.
            A = dx/mag; %Normalized X-component of the direction.
            B = dy/mag; %Normalized Y-component of the direction.
            detA = A*A + B*B; %Determinant of the matrix.
            nx = -B/detA;
            ny = A/detA;
        end%function
        
        function status = NonDecreasing(entries,sequence)
            %Check that a sequence is non-decreasing.
            status = true; %"Innocent until proven guilty."
            for ii = 1:(entries - 1)
                status = status*(sequence(ii) <= sequence(ii + 1));
            end%ii
        end%function
    end%methods (Static)
end
classdef spline  < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        name %Name of the spline's graphics as they appear on axes legend.
        color %Color of spline's graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the spline.
        sketches %Structure containing handles to the spline's graphics.
        refresh_rate %A parameter that determines how long (ideally) calls to Refresh shoul take.
    end%properties (Public)
    %DEFINING DATA variables
    properties (SetAccess = protected)
        kn_seq %Knot sequence.
        knots %Number of knots in the sequence.
        coeff %Coefficients of the polynomial pieces.
        degree %Degree of each polynomial.
        pieces %Number of polynomial pieces.
        dim %Dimension of the interpolated space.
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        %Graphics related
        refresh %Structure holding boolean values that control which graphics refresh.
        updated %Structure holding boolean values that track graphics state of date.
        generated %Structure holding boolean values that track whether certain graphics are generated.
        canvas_set
        
        %Object related.
        valid
    end%properties (Protected)
    %High-level instance CREATION routines.
    methods (Static)
        
        %Constructor
        function this = spline
            %Graphics-related
            this.canvas_set = false;
            this.valid = false;
            this.sketches = struct(...
                'Curve',[],...
                'BreakPoints',[]);
            this.generated = struct(...
                'Curve',false,...
                'BreakPoints',false);
            this.updated = struct(...
                'Curve',false,...
                'BreakPoints',false);
            this.refresh = struct(...
                'Curve',false,...
                'BreakPoints',false);
            
            this.degree = [];
            this.pieces = [];
            this.kn_seq = [];
            this.coeff = [];
        end
        
        %Custom creation routines.
        function spl = CreateN3FromList(knots,KY)
            %Create Natural cubic splines from an XY list.
            spl = spline;
            spl.knots = knots;
            spl.kn_seq = KY;
            
        end%function
        function spl = CreateBasisFromList(knots,kn_seq)
            %Create A Basis spline from scratch
            spl = spline;
            spl.knots = knots;
            spl.kn_seq = kn_seq; %This needs to be sorted.
            spl.degree = knots - 1;
            spl.pieces = spl.degree;
            spl.coeff = zeros(spl.pieces,knots);
            order = knots;
            buffer = BasisTriangleFan(order,kn_seq);
            spl.coeff = buffer{order};
        end%function
        function spl = CreateBasis
        end%function
        
    end%methods (Static 1)
    %High-level instance MODIFICATION and QUERY routines.
    methods
        
    end%methods (Ordinary)
    %Graphical setups
    methods
    end%methods
    %GRAPHICAL SETUP routines.
    methods
        
        %Graphics Refresh Routines
        function Refresh(this)
        end%function
        function RefreshCurve(this)
            
            this.updated.Curve = true;
        end%function
        function RefreshBreakPoints(this)
            
            this.updated.BreakPoints = true;
        end%function
        
        %Graphics generation functions.
        function GenerateCurve(this)
            GenerateDefaultCanvas(this);
            this.sketches.Curve = line(...
                'Parent',this.canvas,...
                'XData',0,...
                'YData',0,...
                'UserData',100,...
                'Color',this.color,...
                'LineWidth',1,...
                'DisplayName',this.name);
            this.generated.Curve = true;
        end%function
        function GenerateBreakPoints(this)
            GenerateDefaultCanvas(this);
            this.sketches.BreakPoints = line(...
                'Parent',this.canvas,...
                'XData',0,...
                'YData',0,...
                'Color',this.color,...
                'LineStyle','none',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',this.color,...
                'Marker','o',...
                'DisplayName',[this.name,'Break Points']);
            this.generated.BreakPoints = true;
        end%function
        function GenerateDefaultCanvas(this)
            if ~this.canvas_set
                warning('Default canvas created!');
                this.canvas = custom_axis;
                this.canvas_set = true;
            end%if
        end%function
    end%methods (Graphics)
    %GRAPHICAL DEMONSTRATION routines.
    methods
        
        function PieglEX2p2
            %Replicates example 2.2 from Les Piegl's "The NURBS Book"
            n_kn = 10;
            sequence = [0,0,0,1,2,3,4,5,5,5];
            p = 2;
            splines = 7;
            spl = cell(splines,1);
            for ii = 1:splines
                iip1 = ii + p + 1;
                iip2 = iip1 + p + 1;
                spl{ii} = CreateBasisFromList(4,sequence(iip1:iip2));
                spl{ii}.GenerateCurve;
                spl{ii}.sketches.Curve.Visible = 'on';
            end%ii
            
        end%function
        
    end%methods (Demonstrations)
    %Low-level SPECIALIZED routines
    methods (Static)
        
        %Specialized Basis Spline routines.
        function buffer = BasisBuffer(order)
            %Creates a cell array that contains arrays for storing Basis
            %Splines. The last array in the cell is the desired Basis
            %Spline by the input order. All other arrays store intermediate
            %Basis Splines used to create the final Basis Spline. This
            %routine only allocates memory. Does not create the splines.
            
            %order = degree + 1; ""
            buffer = cell(order,1); %Need pointers to buffer arrays.
            for ii = 1:order
                buffer{ii} = zeros(order + 1 - ii,order + 1 -ii);
                buffer{ii}(1,1) = 1; %This one is a B^0 spline.
            end%ii
        end%function
        function buffer = CoxDeBoor(buffer,order,kn_seq)
            %Recursive formula for generating higher-order Basis splines
            %from lower-order ones.
            %buffer: memory where polynomial coefficients are stored.
            
            %order: ... of the input basis spline. Degree goes up by one
            %after calling this routine.
            %kn_seq: Sequence of knots used to upgrade the order of the
            %Basis spline.
            for ii = 1:order 
                idx = 1:kk;
                for jj = 1:order
                    den = kn_seq(jj + ii) - kn_seq(jj); %Denominator
                    if den ~= 0 %"Bootstrap" component.
                        buffer{jj}(idx  ,idx  ) = -buffer{jj}(idx  ,idx)*kn_seq(jj)/den;
                        buffer{jj}(idx  ,idx+1) = +buffer{jj}(idx  ,idx+1) - buffer{jj  }(idx,idx)/kn_seq(jj);
                    end%if
                    den = kn_seq(jj+ii+1) - kn_seq(jj+1); %Denominator
                    if den ~= 0 %"Next" spline component.
                        buffer{jj}(idx+1,idx  ) = +buffer{jj}(idx+1,idx  ) + buffer{jj+1}(idx,idx)*kn_seq(jj+1+kk)/den;
                        buffer{jj}(idx+1,idx+1) = +buffer{jj}(idx+1,idx+1) - buffer{jj+1}(idx,idx)/den;
                    end%if
                end%jj
            end%ii
        end%function
        function buffer = BasisTriangleFan(order,kn_seq)
            %Allocates a specialized memory buffer, and builds a Basis
            %Spline of prescribed order using the CoxDeBoor formula. It 
            %does so by building all intermediate lower-order splines
            %according to the "triangular" pattern. The knot sequence is
            %not error checked and must be of a minimum size. Excess knots
            %are ignored.
            
            %order:Desired order of the Basis Spline.
            %kn_seq: Defining knot sequence.
            buffer = BasisBuffer(order); %Allocate space for the Splines
            for ii = 1:(order - 1) %ii = intermediate order of splines.
                buffer = spline.CoxDeBoor(buffer,ii,kn_seq);
            end%ii
        end%function
        function buffer = BasisSlopedFan(buffer,degree)
        end%function
        
    end%methods (Static)
    %Low-level "QUALITY OF LIFE" routines
    methods (Static)
    end
end%classdef
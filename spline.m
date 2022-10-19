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
        graphics_initialized
        
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
                'Curve',[]);
            
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
        function spl = CreateBFromList(knots,KY)
            %Create A Basis spline from scratch
            spl = spline;
            spl.knots = knots;
            spl.kn_seq = KY;
            
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
    end%methods (Graphics)
    %GRAPHICAL DEMONSTRATION routines.
    methods
    end%methods (Demonstrations)
    %Low-level SPECIALIZED routines
    methods (Static)
    end%methods (Static)
    %Low-level "QUALITY OF LIFE" routines
    methods (Static)
    end
end%classdef
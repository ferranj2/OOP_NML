classdef spline  < handle
    properties (SetAccess = public)
        canvas
        sketches
    end%properties (Public)
    properties (SetAccess = protected)
        kn_seq %Knot sequence.
        knots %Number of knots in the sequence.
        coeff %Coefficients of the polynomial pieces.
        degree %Degree of each polynomial.
        pieces %Number of polynomial pieces.
    end%properties (Protected)
    properties (Hidden = true)
        %Graphics related
        canvas_set
        graphics_initialized
        
        %Object related.
        valid
    end%properties (Protected)
    
    %High-level functions that MODIFY specific instances of the class
    %object.
    methods
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
    end%methods (Ordinary)
    
    %Graphical setups
    methods
    end%methods
    
    %High-level functions that CREATE instances of the polynomial objects
    %from low-level code.Error checking involved.
    methods (Static)
        
        %Create Natural cubic splines from an XY list.
        function spl = CreateN3FromList(knots,KY)
            spl = spline;
            spl.knots = knots;
            spl.kn_seq = KY;
            
        end%function
        
        %Create A Basis spline from scratch
        function spl = CreateBFromList(knots,KY)
            spl = spline;
            spl.knots = knots;
            spl.kn_seq = KY;
            
        end%function
        
        
    end%methods (Static 1)
    %Low-level functions with no error checking specific to this class.
    methods (Static)
    
    end%methods (Static 2)
    %Elementary low-level functions that are not unique in application to
    %the class.
    methods (Access = public)
    
    end%methods (Static 2)
end%classdef
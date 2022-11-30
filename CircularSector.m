classdef  CircularSector < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        N %Number of vertices to place on the polyline.
        Polyline %Discrete polyline.
        name %Name of the OBJECT'S graphics as they appear on axes legend.
        color %Color of OBJECT'S graphics as will be rendered in an axes object.
        canvas %Axes object on which to draw the OBJECT.
        sketches %Structure containing handles to the OBJECT'S graphics.
        refresh_rate %A parameter that determines how long (ideally) calls to Refresh shoul take.
    end%properties (Public)
    %DEFINING DATA variables
    properties (SetAccess = protected)
        th1 %First angular limit
        th2 %second angular limit
        C% The defining circle.
    end%properties (Protected)
    %METRIC variables
    properties (SetAccess = protected)
        theta %Angular span of the sector.
        arc %Length of the arc.
        perimeter %Length around the whole sector.
        area %Area of the sector.
        chord %Length between the two points on the defining circle.
        sagitta %Radial length between the circle and the chord.
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        %Graphics related
        refresh %Structure holding boolean values that control which graphics refresh.
        updated %Structure holding boolean values that track graphics state of date.
        generated %Structure holding boolean values that track whether certain graphics are generated.
        canvas_set %Whether the sector has an axes associated with it.
        
        %Object related.
        valid
    end%properties (Protected)
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = CircularSector
            if exist('circle.m','file') ~= 2
                error('The "circle.m" class is a prerequisite for usage of this class.');
            end%if
            
            this.name = 'Circular Sector';
            this.color = [0,0,0];
            this.canvas = [];
            this.refresh_rate = 60;
            this.th1 = [];
            this.th2 = [];
            this.Polyline = polygon.empty;
            this.N = 30;
            
            this.canvas_set = false;
            this.valid = false;
            
            this.sketches = struct(...
                'Curve',[],...
                'Label',[],...
                'LabelArc',[]);
            this.generated = struct(...
                'Label',false,...
                'LabelArc',false);
            this.refresh = struct(...
                'Label',false,...
                'LabelArc',false);
            this.updated = struct(...
                'Label',false,...
                'LabelArc',false);
        end%function
        
        %Custom creation routines.
        function sector = CreateByCuttingCircle(circ,th1,th2)
            %Given an instance of circle.m, define the sector by
            %prescribing angular bounds.
            if th1 == th2
                error('The bounding angles for the circular sector cannot be the same!');
            end%if
            sector = CircularSector;
            sector.C = circ;
            if th1 > th2
                [th1,th2] = CircularSector.swap(th1,th2);
            end%if
            sector.th1 = th1;
            sector.th2 = th2;
            sector.Measure;
        end%function
    end%methods (Ordinary)
    %High-level instance MODIFICATION and QUERY routines.
    methods 
        function Measure(this)
            %Obtain ALL metric properties of the sector.
            this.theta = this.th2 - this.th1;
            this.chord = 2*this.C.R*sin(0.5*this.theta);
            this.sagitta = 0.5*this.chord*tan(this.theta*0.25);
            this.arc = this.C.R*this.theta;
            this.area = this.C.area*this.theta*0.5/pi;
            this.perimeter = this.arc + 2*this.C.R;
        end%function

    end%methods (Static)
    %GRAPHICAL SETUP routines.
    methods
        %"Set" functions
        function InheritFromCircle(this)
            %Have the circular arc inherit as many graphical settings from
            %the progenitor circle as applicable.
            if this.C.canvas_set
                this.SetCanvas(this.C.canvas);
            end%if
            this.color = this.C.color;
            this.name = ['Arc of',this.C.name];
        end%function
        function SetCanvas(this,ax)
            this.canvas = ax;
            this.canvas_set = true;
        end%function
        function SetColor(this,RGB)
            this.color = RGB;
        end%function
        function SetName(this,string)
            this.name = string;
        end%function
        
        function GenerateArcPolyline(this)
            this.Polyline = polygon.CreateRegularByRadius(...
                this.C.R,... %Progenitor circle's radius.
                this.N,...
                this.C.xc,... %Progenitor circle's x-centroid.
                this.C.yc,... %Progenitor circle's y-centroid.
                this.th1,...
                this.th2);
            this.Polyline.SetLineWidth(0.5);
            if this.canvas_set
                this.Polyline.SetCanvas(this.canvas);
                this.Polyline.Toggle('Curve');
            end%if
        end%function
        

    end%methods (Graphics)
    %GRAPHICAL DEMONSTRATION routines.
    methods (Static)
        function [C,S] = Test0
            C = circle.CreateXYR(2,2,4);
            ax = custom_axis;
            axis(ax,'equal');
            C.SetCanvas(ax)
            S = CircularSector.CreateByCuttingCircle(C,pi/4,3*pi/4);
            S.InheritFromCircle;
            S.GenerateArcPolyline;
        end%function
    end%methods (Demonstrations)
    %Low-level SPECIALIZED routines
    methods (Static)
    end%methods
    %Low-level "QUALITY OF LIFE" routines
    methods (Static)
        function [A,B] = swap(A,B)
            buffer = B;
            B = A;
            A = buffer;
            clear buffer;
        end%function
    end%methods
end%classdef
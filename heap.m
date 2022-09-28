classdef heap < handle
    properties (SetAccess = protected)
        length %Elements in the heap array.
        keys %The actual elements.
        base %Numeric integer base 2,3,4... (Cannot be 1 or less.)
    end%properties (protected)
    properties (Hidden = true)
        order %"+1" if max-heap, "-1" if min-heap or "0" if neither.
        valid %Sanity check flag.
    end%properties (Hidden)
    methods
        %Constructor
        function this = heap
            this.length = [];
            this.base = [];
            this.order = [];
            this.keys = [];
        end%function
    end %methods (Ordinary)
    methods (Static)
    end%methods (Static)
    
end%classdef
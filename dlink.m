classdef dlink < handle
    properties (SetAccess = public)
        dimlock
    end
    properties (SetAccess = protected)
        ID
        dim %Number of components.
        components %The components.
    end
    properties (SetAccess = private)
        chain = dchain.empty %Pointer to a chain.
        next = dlink.empty %Pointer to a dlink ahead of this one.
        prev = dlink.empty %Pointer to a dlink behind this one.
    end
    methods
        %% Constructor
        function this = dlink(components)
            if nargin == 0
                this.dim = [];
                this.components = [];
            else
                [rc,cc] = size(components);
                if rc ~= 1 && cc ~= 1
                    error('Input array for dlink components is not a flat array!');
                end
                if rc > cc
                    this.dim = rc;
                else
                    this.dim = cc;
                end
                this.components = components;
           end
        end
        %% Edit con
        function edit(this,varargin)
            %WARNING: "this" is included in nargin.
            if nargin == 1
                error('No dlink properties input!');
            end
            if mod(nargin - 1,2) ~= 0
                error('Must input key-value pairs!');
            end
            for ii = 1:(nargin - 1)/2
                switch varargin{2*ii-1} %"this" is included in the nargin.
                    case 'components'
                        if isnumeric(varargin{2*ii}) ~= 1
                            error('Components must be numeric.');
                        end
                        this.components = varargin{2*ii-1}; 
                    case 'ID'
                        if mod(varargin{2*ii},1) ~= 0
                            error('Input must be an integer.')
                        end
                        this.ID = varargin{2*ii};
                    otherwise
                        error('Unrecognizable dlink property.');
                end
            end
        end
        %% Link to the "next" dlink.
        function Lleftof(this,next_dlink)
            this.next = next_dlink;
            next_dlink.prev = this;
        end
        %% Link to the previous dlink.
        function Lrightof(this,prev_dlink)
            this.prev = prev_dlink;
            prev_dlink.next = this;
        end
        %% Link between two vertices.
        function Lbetween(this,prev_dlink,next_dlink)
            this.prev = prev_dlink;
            prev_dlink.next = this;
            this.next = next_dlink;
            next_dlink.prev = this;
        end

        %% Viewer
        function view(this)
            fprintf('%d-dimensional dlink.\n',this.dim);
            fprintf('ID: %d \n',this.ID);
            fprintf('Components: \n');
            for ii = 1:this.dim
                fprintf('%f\n',this.components(ii));
            end
        end
        %% Destructor
        function delete(this)
            if ~isempty(this.next)
                this.next.prev = dlink.empty;
            end
            if ~isempty(this.prev)
                this.prev.next = dlink.empty;
            end
            this.components = [];
            this.dim = [];
            this.ID = [];
            clear this;
        end
    end
end
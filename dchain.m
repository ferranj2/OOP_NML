classdef dchain < handle
    properties (SetAccess = protected)
        all_vertices_same_dim %BOOLEAN FLAG: Ensure that all vertices in this list are the same dimension.
        all_vertices_dim %INTEGER: If the above boolean flag is used, what dimension should be enforced?
        links
        closed %BOOLEAN FLAG: That indicates whether the dchain loops back into itself.
        needle %Where along the dchain is the pointer.
        length %Length of the dchain.
        first %Pointer to the first member of the dchain.
        last %Pointer to the last member of the dchain.
        current %Vertex corresponding to the needle.
    end
    methods
        %% Constructor
        function this = dchain(varargin)
            if nargin == 0
                this.all_vertices_same_dim = false;
                this.needle = 0;
                this.length = 0;
                this.first = dlink.empty;
                this.last = dlink.empty;
                return;
            end
            %{
            if mod(nargin,2) ~=0
                error('Must input key-value pairs.')
            end
            for ii = 1:nargin/2
                switch varargin{2*ii-1}
                    case 'first'
                    case 'last'
                end
            end
            %}
        end
        %% Link addition
        function append(this,dlink,side)
            if nargin == 2
                side = 'right'; 
            end
            %In case the dchain is uninitialized.
            if isempty(this.last) == 1 && isempty(this.first) == 1
                if isempty(dlink.ID) == 1
                    dlink.edit('ID',1);
                end
                this.last = dlink;
                this.first = dlink;
            else %Chain is initialized.
                if strcmp(side,'right') == 1 
                    this.last.Lleftof(dlink); %Current "last" is left-coupled to the newly added dlink.
                    this.last = dlink; %"Last" pointer is now towards newly added dlink.
                    dlink.edit('ID',dlink.prev.ID + 1);
                elseif strcmp(side,'left') == 1
                    this.first.Lrightof(dlink);
                    this.first = dlink;
                    dlink.edit('ID',dlink.next.ID - 1);
                end
            end
            this.current = dlink; %"Current" Vertex is the newly added one.
            this.length = this.length + 1; %Update the dchain's length.
            this.needle = this.needle + 1; %Set needle to the newly added dlink.
        end
        %% 
        % |Insert between existing entries|
        function displace(this,dlink,side)
        end
        %%
        function unveil(this)
            this.current = this.first; %Point to the first link.
            for ii = 1:this.length 
                fprintf('ID: %-4d \t',this.current.ID);
                for jj = 1:this.current.dim
                    fprintf('%5f\t',this.current.components(jj));
                end
                fprintf('\n')
                this.current = this.current.next;
            end
        end
        %%
        function array = chain2array(this)
            if this.length == 0 %If the dchain is uninitialized.
                array = [];
                return;
            end
            array = zeros(this.length,this.first.dim); %Allocate contiguous array to write dchain to.
            this.current = this.first; %Set dlink pointer to the first dlink.
            for ii = 1:this.length %Row counter.
                for jj = 1:this.current.dim % Column counter.
                    array(ii,jj) = this.current.components(jj);
                end
                this.current = this.current.next;
            end
        end
        %%
        function closechain(this)
            this.first.Lrightof(this.last);
            this.closed = true; %Set the flag to closed.
        end
        %%
        function breakchain(this)
            
        end
        %% Destructor
        function delete(this)
            this.current = this.first;
            this.needle = 1;
            %Destroy the links.
            for ii = 1:this.length
                next = this.current.next;
                this.current.delete;
                this.current = next;                
            end
        end
    end
end
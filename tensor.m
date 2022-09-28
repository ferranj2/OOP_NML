classdef tensor < handle
        properties (SetAccess = private)
            dim
            rank
            data
        end%properties
        properties (SetAccess = protected)
            components
        end%properties (Protected)
        properties (Hidden = true)
            generated
        end%properties (Hidden)
        methods
            %Constructor
            function this = tensor(varargin)
                if nargin == 0
                    this.dim = [];
                    this.rank = [];
                    this.data = [];
                    this.components = [];
                    return;
                end%if
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
                    if strcmp(varargin{ii},'rank') == 1
                        this.rank = varargin{iip1};
                        valid_input = true;
                    end%if
                    if strcmp(varargin{ii},'data') == 1
                        this.data = varargin{iip1};
                        valid_input = true;
                    end%if
                    if ~valid_input
                        error('One of the inputs is unrecognizable!');
                    end%if
                    ii = iip1;
                    valid_input = false;
                end%while
                %this.generated = IsGenerated(this);
            end%function
            
            %Access tensor element
            function el = index(this,idx)
                if length(idx) ~= this.rank
                    error('More indices than the rank of the tensor!')
                end%if
                el_idx = 0;
                for ii = 1:this.rank
                    el_idx = el_idx + idx(ii)*(this.dim^(ii-1));
                end%ii
                el = this.
            end%function
            
            
            
        end%methods
        methods (Static)
            %Create a tensor from dimension,rank, and data.
            function this = Create(data,dim,rank)               
                if isempty(data)
                    error('Cannot input empty data!');
                end%if
                this = tensor;
                if nargin == 1
                    this.rank = 1;
                    this.dim = length(data);
                    this.data = data;
                    this.components = this.dim;
                    return;
                end%if
                this.dim = dim;
                this.rank = rank;
                this.components = tensor.tensor_components(dim,rank);
                %Use of the length function in C can be emulated with
                %pointer arithmetic.
                if length(data) ~= this.components
                    error('Length of data arrays does not match D^R!');
                end%if
                this.data = data; 
                this.generated = true;
            end%function
            
            %Create a random tensor of given dimension and rank.
            function W = CreateRandom(dim,rank)
                if nargin == 1
                    W = tensor.Create(rand(1,dim));
                    return;
                end%if
                n_comp = tensor.tensor_components(dim,rank);
                W = tensor.Create(rand(1,n_comp),dim,rank);
            end%function
            
            %Add two tensors together inplace. (PETSc nomenclature)
            function AXPY(Y,a,X)
                %Add Alpha*X Plus Y (AXPY) inplace
                if nargin < 3
                    for ii = 1:Y.components
                        Y.data(ii) = Y.data(ii) + a;
                    end%ii
                    return;
                end%if
                tensor.samerank(X,B);
                tensor.samedim(X,B);
                for ii = 1:X.components
                    Y.data(ii) = Y.data(ii) + a*X.data(ii);
                end%ii
            end%function
            
            %Create a new tensor by adding two already existing ones.
            function W = WXPY(X,Y)
                tensor.samerank(X,B);
                tensor.samedim(X,B);
                W = Copy(X);
                for ii = 1:X.components
                    W.data(ii) = X.data(ii) + Y.data(ii);
                end%for
            end%function
            
            %Create a new tensor by copying an alrady existing one.
            function W = copy(A)
                W = tensor;
                W.dim = A.dim;
                W.rank = A.rank;
                W.components = A.components;
                W.data = A.data;
            end%function
            
            function W = WAXPBY
            end%function
            
            %Inner product between two tensors
            function dot = inner_product(A,B)
                samerank(A,B);
                samedim(A,B);
                dot = 0;
                for ii = 1:A.components
                    dot = dot + A.data(ii)*B.data(ii);
                end%ii
            end%function
            
            
            %# of components of a tensor.
            function n_comp = tensor_components(dim,rank)
                n_comp = power(dim,rank);
            end%function
            
            
        end%methods (Static)
        
        methods (Access = public)
            
            
            %Ensure tensors are of same dimension.
            function samedim(A,B)
                if A.dim ~= B.dim
                    error('Input tensors must be of same dimension!');
                end%if
            end%function
            
            %Ensure tensors are of same rank
            function samerank(A,B)
                if A.rank ~= B.rank
                    error('Input tensors must be of same rank!');
                end
            end%function
            
        end%methods (public)
        
end%classdef
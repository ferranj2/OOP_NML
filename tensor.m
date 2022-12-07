%% tensor.m
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  Department of Aerospace Engineering (AE)
%% Description
% A data structure used to represent tensors of arbitrary rank and
% dimension. This is meant as a convenience feature for computing and
% storing tensor fields during pre-processing or post-processing. This is
% not meant as a competitor to Linear Algebra abstractions like the good
% old fashioned 'N by N' matrix. Most of the provisions in this class hold
% for rank-2 tensors in 3D space. However, some effort has been put to
% generalize to rank-M and ND space.
%% Class definition
classdef tensor < handle
    %CUSTOMIZATION variables
    properties (SetAccess = public)
        basis_names %Name of the tensor field components.
        name %Name of the tensor field.
    end%properties
    %DEFINING DATA variables
    properties (SetAccess = private)
        comp_arrays %Cell array (**double in C) that stores the tensor components.
        dim %Dimension of the space that the tensor exists in.
        rank %The number of free indices in the tensor.
        entries %Number of copies to store (to "arrayify" the tensor).
        type %Storage format
    end%properties
    %METRIC variables
    properties (SetAccess = protected)
        components %Number of components per tensor.
    end%properties (Protected)
    %FLAG and STATE variables.
    properties (Hidden = true)
        rank_spanner
        rank_sizes %[D^1,D^2,..., D^R]
        first_offset %The integer ID of the first component of a queried tensor.
        count_offset %A counter variable that spans between the first and last integer ID's of the queried tensor.
        final_offset %The integer ID of the last component of a queried tensor.
        query %The # of the tensor (of potentially many arrays) to query.
    end%properties (Hidden)  
    %High-level instance CREATION routines.
    methods (Static)
        %Constructor
        function this = tensor(varargin)
            this.comp_arrays = [];
            this.dim = [];
            this.rank = [];
            this.entries = [];
            this.rank_sizes = [];
            this.first_offset = 0;
            this.count_offset = 1;
            this.final_offset = [];
            this.rank_spanner = 0;
            this.query = 1; %By default, the queries are on the first tensor.
            this.components = [];
            this.name = 'Unnamed Tensor';
        end%function
        
        %Custom Creation routines
        function T = Create(dim,rank,n_entries,data)
            %Create a tensor manually by defining the dimension of the
            %space, the rank, the number of entries if creating a tensor
            %field, and the data. The size of the data is implied by the
            %choice of dimension and rank.
            T = tensor;
            T.SetDimension(dim);
            T.SetRank(rank);
            T.SetEntries(n_entries);
            T.Measure;
            T.SetComponents(data);
            T.SetQueryTensor(1);
        end%function
        function T = CreateRandom(dim,rank,n_entries)
            %Create a dense tensor of given dimension and rank with
            %randomized components.
            if nargin < 3
                n_entries = 1; %At least one tensor must be stored.
            end%if
            T = tensor;
            T.SetRank(rank);
            T.SetDimension(dim);
            T.SetEntries(n_entries);
            T.SetName('Random Tensor');
            T.CallocDense;
            for ii = 1:T.components
                for jj = 1:T.entries
                    T.comp_arrays{ii}(jj) = rand;
                end%jj
            end%ii
            T.SetQueryTensor(1);
        end%function
        function S = CreateScaling(dim,factors)
            %Creates a rank-2 tensor with scalar factors on its diagona
        end%function
        function R = Create2DRotation(theta)
            %Creates a rank-2 direction cosines tensor that encodes a
            %rotation about the origin in the 2D Euclidean plane. The
            %convention is that positive angles produce counter-clockwise
            %rotations. Negative angles produce clockwise rotations.
            R = tensor;
            R.SetRank(2);
            R.SetDimension(2);
            R.SetEntries(1);
            R.SetName('Rotation Direction Cosines');
            R.CallocDense;
            
            c = cos(theta);
            s = sin(theta);
            R.comp_arrays{1}(1) = +c;
            R.comp_arrays{2}(1) = +s;
            R.comp_arrays{3}(1) = -s;
            R.comp_arrays{4}(1) = +c;
            R.SetQueryTensor(1);
            clear c;
            clear s;
        end%function
        function R = Create3DRotationX(theta)
            %Create a rank-2 tensor in R^3 that encodes a single axis
            %rotation about the X-axis (1st basis direction). Units for the
            %angle are radians. Positive angles imply CCW rotations.
            R = tensor;
            R.SetRank(2);
            R.SetDimension(3);
            R.SetEntries(1);
            R.SetName('X-Rotation Direction Cosines');
            R.SetBasisNames({'X','Y','Z'});
            R.CallocDense;
            
            c = cos(theta);
            s = sin(theta);
            R.comp_arrays{1}(1) = 1;
            R.comp_arrays{2}(1) = 0;
            R.comp_arrays{3}(1) = 0;
            R.comp_arrays{4}(1) = 0;
            R.comp_arrays{5}(1) = +c;
            R.comp_arrays{6}(1) = +s;
            R.comp_arrays{7}(1) = 0;
            R.comp_arrays{8}(1) = -s;
            R.comp_arrays{9}(1) = +c;
            clear c;
            clear s;
            R.SetQueryTensor(1);
        end%function
        function R = Create3DRotationY(theta)
            %Create a rank-2 tensor in R^3 that encodes a single axis
            %rotation about the Y-axis (2nd basis direction). Units for the
            %angle are radians. Positive angles imply CCW rotations.
            R = tensor;
            R.SetRank(2);
            R.SetDimension(3);
            R.SetEntries(1);
            R.SetName('Y-Rotation Direction Cosines');
            R.SetBasisNames({'X','Y','Z'});
            R.CallocDense;
                        c = cos(theta);
            s = sin(theta);
            R.comp_arrays{1}(1) = +c;
            R.comp_arrays{2}(1) = 0;
            R.comp_arrays{3}(1) = -s;
            R.comp_arrays{4}(1) = 0;
            R.comp_arrays{5}(1) = 1;
            R.comp_arrays{6}(1) = 0;
            R.comp_arrays{7}(1) = +s;
            R.comp_arrays{8}(1) = 0;
            R.comp_arrays{9}(1) = +c;
            clear c;
            clear s;
            R.SetQueryTensor(1);
        end%function
        function R = Create3DRotationZ(theta)
            %Create a rank-2 tensor in R^3 that encodes a single axis
            %rotation about the Z-axis (3nd basis direction). Units for the
            %angle are radians. Positive angles imply CCW rotations.
            R = tensor;
            R.SetRank(2);
            R.SetDimension(3);
            R.SetEntries(1);
            R.SetName('Z-Rotation Direction Cosines');
            R.SetBasisNames({'X','Y','Z'});
            R.CallocDense;
            
            c = cos(theta);
            s = sin(theta);
            R.comp_arrays{1}(1) = +c;
            R.comp_arrays{2}(1) = +s;
            R.comp_arrays{3}(1) = 0;
            R.comp_arrays{4}(1) = -s;
            R.comp_arrays{5}(1) = +c;
            R.comp_arrays{6}(1) = 0;
            R.comp_arrays{7}(1) = 0;
            R.comp_arrays{8}(1) = 0;
            R.comp_arrays{9}(1) = 1;
            clear c;
            clear s;
            R.SetQueryTensor(1);
        end%function
        function T = Copy(A)
            %Create a new tensor by copying an alrady existing one.
            T = tensor;
            T.SetDimension(A.dim);
            T.SetRank(A.rank);
            T.SetName(['Copy of ',A.name]);
            T.CallocDense;
            for ii = 1:T.components
                for jj = 1:T.entries
                    T.comp_arrays{ii}(jj) = A.comp_arrays{ii}(jj);
                end%jj
            end%ii
            T.SetQueryTensor(1);
        end%function
        
        function C = Contraction(A,idx1,idx2)
            %Given a tensor "A" of rank-2 or above, and two indices upon
            %which to perform a contraction, computes a tensor "C" of
            %rank(C) = rank(A) - 2.
            if A.rank < 2
                error('Contraction is defined for rank-2 tensors or above!');
            end%if
            if nargin == 1 && A.rank == 2
                idx1 = 1;
                idx2 = 2;
            end%if
            if idx1 > A.rank || idx1 < 0
                error('Index 1 is out of bounds!');
            end%if
            if idx2 > A.rank || idx2 < 0
                error('Index 2 is out of bounds!');
            end%if
            C = tensor;
            C.SetRank(A.rank - 2);
            C.SetDimension(A.dim);
            
            
        end%function
        function C = OuterProduct(A,B)
            %Computes the outer tensor product of "A" with "B' to produce a
            %third tensor "C" of rank(C) = rank(A) + rank (B). "A" exists
            %in R^"N" and "B" exists in R^"M" such that N = M.
            tensor.EnforceSameDimension(A,B);
            C = tensor;
            C.SetRank(A.rank + B.rank);
            C.SetDimension(A.dim);
            C.SetName(['Outer Product of ',A.name,' with ',B.name]);
            C.SetEntries(A.entries);
            C.CallocDense;
            
            idx_A = zeros(1,A.rank);
            idx_B = zeros(1,B.rank);
            %The use of the "ones" MATLAB function is avoided because it
            %does not exist in the C standard library.
            for ii = 1:A.rank
                idx_A(ii) = 1;
            end%ii
            for ii = 1:B.rank
                idx_B(ii) = 1;
            end%ii
            
            %Make a function later
            for ii = 1:C.entries %All field entries.
                for jj = 1:A.rank
                    for kk = 1:A.dim
                        %Get the element from the "A" tensor.
                        el_A = A.comp_arrays{tensor.Indicial2Offset(A.dim,A.rank,idx_A)}(ii);
                        for ww = 1:B.rank
                            for zz = 1:B.dim
                                %Get the element from the "B" tensor.
                                offB = tensor.Indicial2Offset(B.dim,B.rank,idx_B);
                                el_B = B.comp_arrays{offB}(ii);
                                
                                %Compute the array index into "C" were to place the
                                %product.
                                C.comp_arrays{tensor.Indicial2Offset(C.dim,C.rank,[idx_A,idx_B])}(ii) = el_A*el_B;
                                idx_B(1) = idx_B(1) + 1;
                            end%zz
                            for zz = 1:(B.rank - 1)
                                if idx_B(zz) > B.dim
                                    idx_B(zz) = 1;
                                    idx_B(zz + 1) = idx_B(zz + 1) + 1;
                                end%if
                            end%zz
                            if idx_B(B.rank) == (B.dim + 1)
                                idx_B(B.rank) = 1;
                            end%if
                        end%ww
                        idx_A(1) = idx_A(1) + 1;
                    end%kk
                    for kk = 1:(A.rank - 1)
                        if idx_A(kk) == (A.dim + 1)
                            idx_A(kk) = 1;
                            idx_A(kk + 1) = idx_A(kk + 1) + 1;
                        end%if
                    end%kk
                end%jj
                clear idx_A;
                clear idx_B;
                
              
            end%ii
            %Make a function later
            
            C.SetQueryTensor(1);
        end%function
        function C = InnerProduct(A,B,idxA,idxB)
            %Computes the inner tensor product of "A" with "B' to produce a
            %third tensor "C". The rank of C is rank(A) + rank(B) - 2. The
            %input for "A" can be a tensor field or an individual tensor.
            %"B" must be a single tensor or a tensor field of equal entries
            %to "A". The operation involves a contraction that involves one
            %indes from "A" and another from "B". These indices are
            %specified via "idxA" and "idxB".
            
            tensor.EnforceSameDimension(A,B);
            C = tensor;
            C.SetRank(A.rank + B.rank - 2);
            C.SedDimension(A.dim);
            C.CallocDense;
            
            %Make a function later
            for ii = 1:C.entries %All field entries.
                for jj = 1:C.components %All components
                    C.comp_arrays{jj}(ii) = 0;
                end%jj
            end%ii
            %Make a function later
            
            C.SetQueryTensor(1);
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
        
    end%methods
    %High-level instance MODIFICATION and QUERY routines.
    methods
        %Memory allocation.
        function Measure(this)
            %Measure basic integer characteristics of the tensor based on
            %its dimension and rank.
            this.components = power(this.dim,this.rank);
            this.rank_sizes = zeros(this.rank+1,1);
            this.rank_sizes(1) = 1;
            for ii = 2:this.rank+1
                this.rank_sizes(ii) = this.dim*this.rank_sizes(ii-1);
            end%ii
        end%function
        function CallocDense(this)
            %Allocates memory for a dense tensor instance.
            if ~this.Authentic
                error('The tensor''s definition is innapropriate.')
            end%if
            this.Measure;
            this.comp_arrays = cell(this.components,1);
            for ii = 1:this.components
                this.comp_arrays{ii} = zeros(1,this.entries);
            end%ii
            this.type = 'dense';
        end %function
        
        %"Set" functions
        function SetRank(this,rank)
            %Sets the rank of the tensor.
            if rank < 0
                error('Input rank cannot be negative!');
            end%if
            this.rank = rank;
        end%function
        function SetName(this,string)
            this.name = string;
        end%function
        function SetBasisNames(this,strings)
           this.basis_names = strings; 
        end%function
        function SetDimension(this,dim)
            %Sets the dimension of the space that the tensor lives in.
            if dim < 1
                error('Input dimension cannot be less than 1!')
            end%if
            this.dim = dim;
        end%function
        function SetArray(this,in_array)
            %Associates data with the tensor.
            this.array = in_array;
        end%function
        function SetEntries(this,n_entries)
            if n_entries < 1
                warning('Cannot set a negative number of entries. Defaulting to 1.')
                n_entries = 1;
            end%
            this.entries = n_entries;
        end%function
        function SetComponents(this,data)
            %In C, the following check is implicit because data will be of
            %the "**double" type.
            if ~iscell(data)
                error('Input data must be a flat cell array!');
            end%if
            this.comp_arrays = data;
        end%function
        function SetQueryTensor(this,n_query)
            %Sets the active entry if the tensor is instantiated as an
            %array of tensors for individual operations.
            if n_query > this.entries
                error('Query number exceeds number of entries.');
            end%if
            if n_query < 0
                error('Cannot set a negative query number.');
            end%if
            this.query = n_query;
        end%function
        function SetStorageFormat(this,format)
            %Sets the storage format for keeping the tensor in memory.
            this.type = format;
        end%function
        
        %"Get" functions
        function el = GetElement(this,varargin)
            %Get an element from the tensor's data array correspoding
            %to the input indices in "varargin."
            indices = length(varargin);
            if  indices ~= this.rank
                error('Number of query indices does not match tensor rank.')
            end%if
            for ii = 1:indices
                if varargin{ii} > this.dim
                    error(['Index #',num2str(ii),' exceeds tensor''s dimenion.'])
                end%if
            end%ii
            mult = 1;
            idx = 1;
            R = length(varargin); %rank
            for ii = 1:R
                idx = idx + mult*(varargin{R - ii + 1}-1);
                mult = mult*this.dim;
            end%%ii
            el = this.comp_arrays{idx}(this.query);            
        end%function
        
        %Visualization functions
        function PrintDense(this)
            if this.count_offset == 1
                fprintf('Entry %i out of %i of rank-%i tensor in R^%i "%s":\n',this.query,this.entries,this.rank, this.dim,this.name);
            end%if
            if this.count_offset <= this.components
                if this.rank > 2
                    if mod(this.count_offset - 1,this.rank_sizes(3)) == 0
                        this.rank_spanner = this.rank_spanner + 1;
                        fprintf(['(',num2str(this.rank_spanner)]);
                        fprintf(',j,k)\n');
                    end%if
                end%if
                fprintf('[');
                for jj = 1:this.dim-1
                    fprintf('%5.4f,',this.comp_arrays{jj+this.count_offset-1}(this.query));
                end%jj
                fprintf('%5.4f]\n',this.comp_arrays{jj+this.count_offset}(this.query));
                this.count_offset = this.count_offset + this.dim;
                this.PrintDense;%Recursive call;
            else
                this.count_offset = 1;
                this.rank_spanner = 0;
            end%if
        end%function
        
        %Validation
        function valid = Authentic(this)
            %Make sure that the defining settings for the tensor are
            %valid.
            infractions = 0;
            
            %Validate the rank setting.
            if isempty(this.rank)
                warning('Tensor does not have a defined rank.');
                infractions = infractions + 1;
            elseif this.rank < 0
                warning('Tensor cannot have a negative rank');
                infractions = infractions + 1;
            end%if
            
            %Validate the dimenion setting.
            if isempty(this.dim)
                warning('Tensor does not have a defined dimension.');
                infractions = infractions + 1;
            elseif this.dim < 1
                warning('Dimension of the tensor''s space cannot be less than 1!')
                infractions = infractions + 1;
            end%if
            
            %Validate the "copies" setting.
            if isempty(this.entries)
                warning('Number of tensor copies was not specified, defaulting to 1');
                this.entries = 1;
            elseif this.entries < 1
                warning('Number of tensor copies cannot be negative!');
                infractions = infractions + 1;
            end%if
            
            valid = infractions == 0;
        end%function
        
        %Modification functions
        function AXPYOneToMany(Y,a,X)
            %Add two tensors together inplace. (PETSc nomenclature)
            %Add Alpha*X Plus to Y (AXPY) inplace
            if nargin < 3
                for ii = 1:Y.components
                    for jj = 1:Y.entries
                        Y.comp_arrays{ii}(jj) = Y.comp_arrays{ii}(jj) + a;
                    end%jj
                end%ii
                return;
            end%if
            tensor.EnforceSameRank(X,B);
            tensor.EnforceSameDimension(X,B);
            for ii = 1:Y.components
                for jj = 1:Y.entries
                    Y.comp_arrays{ii}(jj) = Y.comp_arrays{ii}(jj) + a*X.comp_arrays{ii}(jj);
                end%jj
            end%ii
        end%function
        
    end%methods (Static)
    %Low-level ELEMENTARY routines
    methods (Static)
        
        function idx = TensorIndices3ArrayIndex(dim,varargin)
            %Convert the indices of a tensor to an equivalent array
            %index should the components of said tensor be stored in a
            %contiguous memory chunk.
            mult = 1;
            idx = 1;
            R = length(varargin); %rank
            for ii = 1:R
                idx = idx + mult*(varargin{R - ii + 1}{1}-1);
                mult = mult*dim;
            end%%ii
        end%function
        function idx = ArrayIndex2TensorIndices(index,rank,rank_sizes)
            idx = zeros(1,rank); %Allocate memory for output.
            for ii = 1:rank
                rem = mod(index,rank_sizes(rank  + 1 - ii));
                index = index - rem;
                idx(ii) = (index)/rank_sizes(rank  + 1 - ii);
            end%ii
        end%function
        function components = RankAndDim2Comp(rank,dim)
            components = power(dim,rank);
        end%function
        function e_ijk = LeviCivita(i,j,k)
            %Evaluates the rank-3 LeviCivita tensor.
            e_ijk = 0.5*(i - j)*(j - k)*(k - i);
        end%function
        function d_ij = KroneckerDelta(i,j)
            %Evaluates the kronecker delta symbol.
            d_ij = (i == j);
        end%function
                
        
        
        function idx = Indicial2Offset(dim,rank,indices)
            %Given the dimension of the space, the rank of the tensor, and
            %the indices of the query element returns the index into the
            %array storing the elements.
            mult = 1;
            idx = 1;
            for ii = 1:rank
                idx = idx + mult*(indices(rank - ii + 1)-1);
                mult = mult*dim;
            end%%ii
        end%function
        
        function EnforceSameDimension(A,B)
            %Ensure tensors are of same dimension.
            if A.dim ~= B.dim
                error('Input tensors must be of same dimension!');
            end%if
        end%function
        function EnforceSameRank(A,B)
            %Ensure tensors are of same rank
            if A.rank ~= B.rank
                error('Input tensors must be of same rank!');
            end%if
        end%function
        function EnforceSameEntries(A,B)
            %Ensure tensors are of same rank
            if A.entries ~= B.entries
                error('Input tensors must have same number of entries!');
            end%if
        end%function
    end%methods (public)
    %Unit Tests
    methods (Static)
        function [Q1,Q2,Q12] = T00
            %Creates the rank-2 direction cosines tensor in R^2 for some
            %angle. Makes a copy of said tensor, and computes the Outer
            %Product of the original tensor with its copy. 
            Q1 = tensor.Create2DRotation(pi/3);
            fprintf('Rank of Q1 %i\n',Q1.rank);
            Q1.PrintDense;
            Q2 = tensor.Copy(Q1);
            fprintf('Rank of Q2 %i\n',Q2.rank);
            Q2.PrintDense;
            Q12 = tensor.OuterProduct(Q1,Q2);
            Q12.PrintDense
            
        end%functions
        function [A,B,AB] = T01
            %Tests the Outer Product of two rank-1 tensors in R^3. The
            %operation is akin to the following matrix multiplication:
            %[1][4,5,6]   [ 4, 5, 6]
            %[2]        = [ 8,10,12]
            %[3]          [12,15,18]
            dim = 3;
            
            %Some rank-1 tensor "A".
            rank_A = 1;
            components_A = {1,2,3};
            A = tensor.Create(...
                dim,...
                rank_A,...
                1,...
                components_A);
            A.SetName('A');
            
            %Some rank-1 tensor "B".
            rank_B = 1;
            components_B = {4,5,6};
            B = tensor.Create(...
                dim,...
                rank_B,...
                1,...
                components_B);
            B.SetName('C');
            
            AB = tensor.OuterProduct(A,B);
            
        end%function
        function [sigma,normal] = T02
            %Manually create a rank-2 tensor in R^2 and compute its inner
            %product with a manually created rank-1 tensor.
            dim = 3;
            basis_names = {'X','Y','Z'};
            
            %Some stress tensor.
            rank_sigma = 2;
            components_sigma = {...
                500;500;800;...
                500;1000;-750;...
                800;-750;-300};
            sigma = tensor.Create(...
                dim,...
                rank_sigma,...
                1,...
                components_sigma);
            sigma.SetName('Cauchy Stress Tensor');
            
            %Some normal direction.
            rank_normal = 1;
            components_normal = {0.5;0.5;1/sqrt(2)};
            normal = tensor.Create(...
                dim,...
                rank_normal,...
                1,...
                components_normal);
            normal.SetName('Normal');
            %Compute the traction in the direction of the normal.
            
            
        end%function
    end%methods
end%classdef
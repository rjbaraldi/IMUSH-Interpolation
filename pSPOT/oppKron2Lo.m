classdef oppKron2Lo < oppSpot
    %oppKron2Lo  kronecker tensor product to act on a distributed vector.
    %
    %   oppKron2Lo(A,B, gather)A and B are Spot operators or numeric 
    %   matrices. Optional param gather specifies whether the output vector
    %   should be gathered to local lab.
    %   GATHER = 0 will not gather
    %   GATHER = 1 will gather the results of forwards or adjoint
    %   multiplication.
    %   GATHER = 2 will gather only in forward mode.
    %   GATHER = 3 will gather only in backwards (adjoint) mode.
    %
    %       %% Example: Defining seperable sparsity transforms over 
    %       different axes. Here we define a sparsity transform S that 
    %       performs Wavelet analysis on the first dimension
    %       and a 2D Curvelet analysis on the second & third dimension
    %             dim=[64,32,32];
    %             C = opCurvelet(dim(2),dim(3));
    %             W = opWavelet(dim(1),1);
    %             S = oppKron2Lo(C,W,1);
    %
    %       % Make a random 3d data-array
    %       D = distributed.randn(dim(1),prod(dim(2:end)));
    %
    %       % Check to see if the analysis followed by synthesis returns 
    %       the original signal
    %       norm(D(:)-S'*S*D(:))
    %
    %   note that the second operator is applied to x and then the first
    %   operator to the transpose of the result, and so x should be of
    %   dimemsions [cols(op2),cols(op1)], and vectorized after distribution
    %   so it is distributed along the columns evenly.
    %
    %   *Now oppKron2Lo also supports local x vectors and distributes them
    %   before calculation( this distribution will be faster than using 
    %   matlabs (:) function).
    %
    %   See also: oppKron2Lo.drandn, oppKron2Lo.rrandn, oppKron2Lo.dzeros,
    %   oppKron2Lo.rzeros
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        tflag = 0;
        permutation; %Permutation vector of intergers defining the order to
        %use when the operators (children) of the Kronecker product are
        %applied to a data vector.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = oppKron2Lo(varargin)
            
            % Setup and extract gather
            gat = 0;
            if isscalar(varargin{end})
                assert(any(varargin{end} == [0 1 2 3]),...
                    'Gather must be 0 1 2 or 3')
                gat = varargin{end};
                varargin(end) = [];
            end
            
            % Check for number of operators
            if length(varargin) ~= 2
                error('Must specify only 2 operators!')
            end
            
            % Standard pSpot Check
            [opList,m,n,cflag,linear] = ...
                pSPOT.utils.stdpspotchk(varargin{:});
            m = prod(m);
            n = prod(n);
            
            % Construct operator
            op = op@oppSpot('pKron', m, n);
            op.cflag       = cflag;
            op.linear      = linear;
            op.sweepflag   = true;
            op.children    = opList;
            op.gather      = gat;
            op.permutation = (1:length(opList));
            
            %Evaluate the best permutation to use when a multiplication is
            %applied
            if ~ (m == 0 || n == 0)
                op.permutation = op.best_permutation();
            end            
            
            % Setting up implicit dimensions of output vector
            op.ms = fliplr(cellfun(@(x) size(x,1),varargin)); % Flipped
            op.ns = fliplr(cellfun(@(x) size(x,2),varargin));
            clear varargin;
            
        end % Constructor
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str=['pKron(',char(op.children{1})];
            
            % Get operators
            for i=2:length(op.children)
                str=strcat(str,[', ',char(op.children{i})]);
            end
            str=strcat(str,')');
            if op.tflag
                str = strcat(str, '''');
            end
        end % Char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mtimes is overloaded so as to call multiplication on a
        % distributed array. This multiplication will do the expected 2D
        % transform on 'x'.
        % For the moment mtimes is only implemented for right
        % multiplication
        function y=mtimes(op,x)
            if isa(x,'dataContainer')
                y = mtimes(x,op,'swap');
            else
                
                if ~isa(op,'oppKron2Lo')
                    error('Left multiplication not taken in account')
                elseif isa(x,'opSpot')    
                    y = opFoG(op,x);
                elseif ~isa(x,'oppKron2Lo')
                    assert( isvector(x) , 'Please use vectorized matrix')
                    op.counter.plus1(op.tflag + 1 );
                    y=op.multiply(x, 1 ); % use tflag to determine mode
                                          % within multiply
                else
                    error(['unsupported data type: ' class(x)]);
                end
            end % catch
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transpose is overloaded to avoid wrapping the operator in an
        % opTranspose.
        function y = transpose(op)
            [m,n] = size(op);
            y = op;
            y.m = n;
            y.n = m;
            y.tflag =  ~op.tflag;
            y.permutation = op.permutation(end:-1:1);
        end
        function y = ctranspose(op)
            y = transpose(op);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double is overloaded to use a vectorized identity matrix
        function y = double(op, enable)
            if nargin < 2 || ~enable
                error(['oppKron2Lo is intended for large applications.' ...
                    ' The explicit representation will likely be very '...
                    ' large, \nuse double(x,1)  to proceed anyway']);
            end
            y = double(kron(op.children{1},op.children{2}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn/rrandn/zeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn is overloaded to create a distributed random vector
        function y = drandn(op)
        %DRANDN Random vector in operator domain
            if ~op.tflag
                dims = [op.children{2}.n, op.children{1}.n];
            else
                dims = [op.children{2}.m, op.children{1}.m];
            end
            y = distrandnvec( dims );
        end
        function y = rrandn(op)
            %RRANDN Random vector in operator range
            if ~op.tflag
                dims = [op.children{2}.m, op.children{1}.m];
            else
                dims = [op.children{2}.n, op.children{1}.n];
            end
            y = distrandnvec( dims );
        end
        function y = dzeros(op)
            %DZEROS Zero vector in operator domain
            if ~op.tflag
                dims = [op.children{2}.n, op.children{1}.n];
            else
                dims = [op.children{2}.m, op.children{1}.m];
            end
            y = distzeros( dims );
        end
        function y = rzeros(op)
            %RZEROS Zero vector in operator range
            if ~op.tflag
                dims = [op.children{2}.m, op.children{1}.m];
            else
                dims = [op.children{2}.n, op.children{1}.n];
            end
            y = distzeros( dims );
        end
        
    end% Methods
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            %The Kronecker product (KP) is applied to the right-hand matrix
            %taking in account the best order to apply the operators A and
            %B.
            
            %Operators
            A=op.children{1};
            B=op.children{2};
                        
            %we could have been called through opSpot.applyMultiply, then
            %we need to pay attention to mode.
            if mode == 2 && ~op.tflag || mode ==1 && op.tflag
                A = A';
                B = B';
            end
            
            %Size of the operators
            [rA,cA]=size(A);
            [rB,cB]=size(B);
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            perm = op.permutation; %Permutation to take in account.
            skipA = A.isDirac;
            skipB = B.isDirac;
            
            if perm(1)==2 %Classic multiplication order
                spmd
                    % Setup partition of columns
                    partition = codistributed.zeros(1,numlabs);
                    
                    if iscodistributed(x)
                        y=getLocalPart(x);
                        local_width=length(y)/cB;
                        assert( mod(local_width,1) == 0, ...
                            ' x must be distributed along cols before vec')
                        y = reshape(y,cB,local_width); % reshape to local 
                        partition(labindex) = local_width; % matrices
                    else
                        y = reshape(x,cB,cA);
                        y = codistributed(y);
                        y = getLocalPart(y);
                        local_width = size(y,2);
                        partition(labindex) = local_width;
                    end
                    if ~skipB
                        y=B*y;% apply B to local matrices
                    end
                    if ~skipA
                        y=y.'; % Tranpose
                        % Build y distributed across rows
                        y = codistributed.build(y, codistributor1d...
                            (1,partition,[cA,rB]));
                        y = redistribute(y,codistributor1d(2));%distributed
                        %along cols
                        y = getLocalPart(y);
                        y=A*y;%apply A to local matrices, then transpose
                        y=y.'; % Now y is distributed across rows again
                        % Rebuild y as distributed across rows
                        y=codistributed.build(y,codistributor1d(1,...
                            codistributor1d.unsetPartition,[rB,rA]));
                        % Redistribute y across columns
                        y=redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    % now vectorize y
                    local_size = numel(y);
                    partition(labindex) = local_size;
                    y = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        partition, [rA*rB,1]));
                end
            else  %Inverted multiplication order
                spmd
                    % Setup partition of columns
                    partition = codistributed.zeros(1,numlabs);
                    
                    if iscodistributed(x)
                        y=getLocalPart(x);
                        local_width=length(y)/cB;
                        assert( mod(local_width,1) == 0, ...
                            'x must be distributed along cols before vec')
                        y = reshape(y,cB,local_width); % reshape to local 
                        partition(labindex) = local_width; % matrices
                        
                        if ~skipA
                            y=y.'; % transpose since A is applied first
                            % Build y distributed across rows
                            y = codistributed.build(y, codistributor1d...
                                (1,partition,[cA,cB]));
                            % Redistribute y across cols
                            y = redistribute(y,codistributor1d(2));
                        end
                    else
                        y = reshape(x,cB,cA);
                        if ~skipA
                            y = y.';
                            y = codistributed(y);
                        end
                    end
                    
                    if ~skipA
                        y = getLocalPart(y);
                        y = A*y;%apply A to local matrices, then transpose
                        y = y.';
                        
                        % Rebuild y distributed across rows
                        y = codistributed.build(y,codistributor1d(1,...
                            codistributor1d.unsetPartition,[cB,rA]));
                        % Redistribute y across cols
                        y=redistribute(y,codistributor1d(2));
                        y = getLocalPart(y);
                    end
                    
                    if ~skipB
                        y = B*y;% apply B to local matrices, no need to 
                    end         % transpose
                    
                    %now vectorize y
                    local_size = numel(y);
                    partition(labindex) = local_size;
                    y = y(:);
                    y = codistributed.build(y, codistributor1d(1,...
                        partition, [rA*rB,1]));
                end
            end
            % if op.gather, y = gather(y); end %#ok<PROP,CPROP>
            if mode == 2 && ~op.tflag || mode ==1 && op.tflag % this is the
                if op.gather == 1 || op.gather == 3           % correct
                    y = gather(y);                            % adjoint
                end                                           % case
            else % this is the forward case
                if op.gather == 1 || op.gather == 2
                    y = gather(y);
                end
            end % gather
        end % Multiply
    end %Protected methods
    
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % best_permutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Returns the best permutation associated to this Kronecker product
        function perm=best_permutation(op)
            list=op.children; %List of 'op''s children
            cost=zeros(1,length(list)); %Computational costs of the
            %operators (children of 'op'). This is simply a numeric
            %representation of theirs shapes, which will affect computation
            %time. Operators with low computational costs should be applied
            %first.
            for i=1:length(list)
                %Cost = (nbr_rows-nbr_columns) / (size of the operator)
                cost(1,i)=(size(list{i},1)-size(list{i},2))/...
                    (size(list{i},1)*size(list{i},2));
            end
            
            perm=op.quicksort(cost,1,length(cost),op.permutation);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % quick_sort
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Function doing a quick sort on the vector containing the
        %computational costs associated to the operators of the Kronecker
        %product. The corresponding permutation 'perm' is returned as an
        %output. It contains the indices of the operators which have to be
        %successively applied to the data vector. These laters are
        %bracketed from left to right.
        
        %n: permutation enabling to follow the transpositions during the
        %recursive application of the quick sort function.
        %n initialy rates [1,2,..,n] where n is the number of operators in
        %the Kronecker product.
        
        %start and stop: indices of the sort area in the cost vector.
        
        function perm=quicksort(op,cost,start,stop,n)
            
            if start<stop
                left=start;
                right=stop;
                pivot=cost(start);
                
                while 1
                    while cost(right)>pivot,right=right-1;
                    end
                    if cost(right)==pivot && right>start
                        right=right-1;
                    end
                    while cost(left)<pivot,left=left+1;
                    end
                    
                    if(left<right)
                        temp=cost(left);
                        cost(left)=cost(right);
                        cost(right)=temp;
                        
                        temp=n(left);
                        n(left)=n(right);
                        n(right)=temp;
                    else break
                    end
                end
                n=op.quicksort(cost, start, right,n);
                n=op.quicksort(cost, right+1, stop,n);
            end
            perm=n;
        end
        
    end % Private methods
end % Classdef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distrandnvec helper funtion to create a vectorized distributed
% normal random matrix
function y = distrandnvec( sz )
spmd
    y = codistributed.randn(sz);
    codist = getCodistributor(y);
    part = codist.Partition;
    part = part.*sz(1);
    y = getLocalPart(y);
    y = y(:);
    y = codistributed.build(y, codistributor1d(1,part,...
        [sz(1)*sz(2),1]));
end
end

% distzeros helper funtion to create a vectorized distributed
% zeros matrix
function y = distzeros( sz )
spmd
    y = codistributed.zeros(sz);
    codist = getCodistributor(y);
    part = codist.Partition;
    part = part.*sz(1);
    y = getLocalPart(y);
    y = y(:);
    y = codistributed.build(y, codistributor1d(1,part,...
        [sz(1)*sz(2),1]));
end
end


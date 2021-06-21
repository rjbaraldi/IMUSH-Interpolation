classdef oplWindow1Davg < opSpot
%oplWindow1Davg tapered windowing for CARP-SG method
%
%   oplWindow1Davg(N,P,H)
%
%   ARGUMENTS:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%
%   Notes:
%       1. This is not a parallel/distributed operator
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        p = 0;
	h = 0;
	oshape = 0;
	yshape = 0;
	xshape = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oplWindow1Davg(varargin)
	  assert(nargin==3,'lWindow1Davg: wrong # of arguments')
	  n = varargin{1};
	  p = varargin{2};
	  h = varargin{3};
          [ m os ys xs ] = pSPOT.pWindow.funWindowShape1D( n, p, h );
	  op = op@opSpot('lWindow1Davg',m,n);
	  op.p = p;
	  op.h = h;
	  op.oshape = os;
	  op.yshape = ys;
	  op.xshape = xs;
       end % function oplWindow1Davg
       
       % xtratests
       function result = xtratests(op)
           T = 14;
	   x0=rand(op.n,1);
	   y=op*x0;
	   x1=op'*y;
	   check=norm(x1-x0);
	   if check < op.n*10^-T
               result = true;
	   else
               result = false;
	   end
       end % xtratests

       % utest ( skip dottest )
       function output = utest(op,k,verbose)
           try
               addpath(fullfile(spot.path,'tests','xunit'))
           catch ME
               error('Can''t find xunit toolbox.')
           end
           if nargin < 3, verbose = 0; end
           if nargin < 2, k = 5; end
           assertTrue(op.xtratests,k);
           output = 'PASSED!';
       end % utest

    end % Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if (mode == 1)
	      [ A ] = pSPOT.pWindow.funWindow1DavgFor(op.n,op.p,op.h);
              y = A * x;
           else
	      [ B ] = pSPOT.pWindow.funWindow1DavgBck(op.n,op.p,op.h);
              y = B * x;
           end
        end % Multipy
      
    end % Methods
        
end % Classdef

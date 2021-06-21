classdef opLBFGS < opSpot
%OPLBFGS   Maintain a limited-memory BFGS approximation.
%
%   opLBFGS(n, mem) creates an n-by-n operator that performs
%   matrix-vector multiplication with a limited-memory BFGS
%   approximation with memory m >= 1.
%
%   By default, the operator acts as an inverse L-BFGS approximation,
%   i.e., its inverse is an approximation of the Hessian. It is used
%   as follows:
%
%   B = opLBFGS(n, mem);
%   B = update(B, s, y);
%   d = - B \ g;          % Apply inverse L-BFGS.
%
%   The operator may also be used in forward mode, i.e., as an
%   approximation to of the Hessian. In this case, the attribute
%   update_forward should be set to true, as forward mode incurs
%   additional computational cost. It is used as follows:
%
%   B = opLBFGS(n, mem);
%   B.update_forward = true;
%   B = update(B, s, y);
%   d = - B \ g;          % Apply inverse L-BFGS.
%   Bx = B * x;           % Apply forward L-BFGS.

%   D. Orban, 2014.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        mem;
        s;       % Array of s vectors.
        y;       % Array of y vectors.
        ys;      % Array of s'y products.

        alpha;   % Multipliers (for inverse L-BFGS)

        a;       % Negative curvature components of forward L-BFGS
        b;       % Positive curvature components of forward L-BFGS

        insert;  % Current insertion point.
        gamma;   % Scaling factor.
    end

    properties (SetAccess = public)
        update_forward;  % Whether or not to update forward L-BFGS.
        scaling;
        updates;         % number of update attempts
        rejects;         % number of rejected updates
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opLBFGS(n, mem)
       %opLBFGS  Constructor.
          if nargin == 1
             mem = 1;
          end
          if nargin > 2
             error('At most one argument can be specified.')
          end

          % Check if input is an integer
          if ~(isnumeric(mem) || mem ~= round(mem))
             error('Memory parameter must be an integer.');
          end

          % Create object
          op = op@opSpot('L-BFGS', n, n);
          op.cflag  = false;
          op.sweepflag  = true;
          op.mem = min(max(mem, 1), n);
          op.s = zeros(n, op.mem);
          op.y = zeros(n, op.mem);
          op.ys = zeros(op.mem, 1);
          op.alpha = zeros(op.mem, 1);
          op.a = sparse(n, mem);
          op.b = sparse(n, mem);
          op.update_forward = false;
          op.insert = 1;
          op.scaling = false;
          op.gamma = 1;
          op.updates = 0;
          op.rejects = 0;
       end % function opLBFGS
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       function op = set.update_forward(op, val)
          if val && ~op.update_forward
            op.a = zeros(size(op.s));
            op.b = zeros(size(op.s));
          end
          op.update_forward = val;
       end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Must use H = update(H, s, y)...
       % How do you get the syntax H.update(s,y) to work???

       function op = update(op, s, y)
       %store  Store the new pair {s,y} into the L-BFGS approximation.
       %       Discard oldest pair if memory has been exceeded.

         op.updates = op.updates + 1;
         ys = dot(s, y);
         if ys <= 1.0e-20
           warning('L-BFGS: Rejecting {s,y} pair')
           op.rejects = op.rejects + 1;
         else

           op.s(:, op.insert) = s;
           op.y(:, op.insert) = y;
           op.ys(op.insert) = ys;

           if op.scaling
             op.gamma = ys / (y' * y);
           end

           % Update arrays a and b used in forward products.
           if op.update_forward
             op.b(:, op.insert) = y / sqrt(ys);

             for i = 1 : op.mem
               k = mod(op.insert + i - 1, op.mem) + 1;
               if op.ys(k) ~= 0
                 op.a(:, k) = op.s(:, k) / op.gamma;
                 for j = 1 : i - 1
                   l = mod(op.insert + j - 1, op.mem) + 1;
                   if op.ys(l) ~= 0
                     op.a(:, k) = op.a(:, k) + (op.b(:, l)' * op.s(:, k)) * op.b(:, l);
                     op.a(:, k) = op.a(:, k) - (op.a(:, l)' * op.s(:, k)) * op.a(:, l);
                   end
                 end
                 op.a(:, k) = op.a(:, k) / sqrt(op.s(:, k)' * op.a(:, k));
               end
             end
           end

           % Update next insertion position.
           op.insert = mod(op.insert, op.mem) + 1;
         end
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function x = double(op)
       %double  Convert operator to a double.

          % Can't do op * eye(n), but can do op \ eye(n).
          e = zeros(op.n, 1);
          x = zeros(op.n);
          for i = 1 : op.n
            e(i) = 1;
            x(:, i) = op * e;
            e(i) = 0;
          end
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function d = diagonal(op)
       %diag  Extract the diagonal of the limited-memory approximation

         if ~op.update_forward
           error('L-BFGS: not using forward mode. Set update_forward = true.');
         end

         d = ones(op.n, 1);
         if op.scaling
           d = d / op.gamma;
         end

         for i = 1 : op.mem
           k = mod(op.insert + i - 2, op.mem) + 1;
           if op.ys(k) ~= 0
             for j = 1 : op.n
               d(j) = d(j) + op.b(j, k)^2 - op.a(j, k)^2;
             end
           end
         end
       end % function diagonal
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % Methods


    methods ( Access = protected )

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function q = multiply(op, x, ~)
       %multiply  Multiply operator with a vector.
       % See, e.g., Nocedal & Wright, 2nd ed., Procedure 7.6, p. 184.

         if ~op.update_forward
           error('L-BFGS: not using forward mode. Set update_forward = true.');
         end

         q = x;
         if op.scaling
           q = q / op.gamma;
         end

         % B = B0 + ∑ (bb' - aa').
         for i = 1 : op.mem
           k = mod(op.insert + i - 2, op.mem) + 1;
           if op.ys(k) ~= 0
             q = q + (op.b(:, k)' * x) * op.b(:, k)- (op.a(:, k)' * x) * op.a(:, k);
           end
         end
       end % function multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function q = divide(op, b, ~)
       %divide  Solve a linear system with the operator.
       % See, e.g., Nocedal & Wright, 2nd ed., Algorithm 7.4, p. 178.

         q = b;

         for i = 1 : op.mem
           k = mod(op.insert - i - 1, op.mem) + 1;
           if op.ys(k) ~= 0
             op.alpha(k) = (op.s(:, k)' * q) / op.ys(k);
             q = q - op.alpha(k) * op.y(:, k);
           end
         end

         if op.scaling
           q = q * op.gamma;
         end

         for i = 1 : op.mem
           k = mod(op.insert + i - 2, op.mem) + 1;
           if op.ys(k) ~= 0
             beta = (op.y(:, k)' * q) / op.ys(k);
             q = q + (op.alpha(k) - beta) * op.s(:, k);
           end
         end
       end % function divide
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

end % Classdef

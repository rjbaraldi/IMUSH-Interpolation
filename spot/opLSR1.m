classdef opLSR1 < opSpot
  %OPLSR1   Maintain a limited-memory SR1 approximation.
  %
  %   opLSR1(n, mem) creates an n-by-n operator that performs
  %   matrix-vector multiplication with a limited-memory SR1
  %   approximation with memory m >= 1.
  %
  %   By default, the operator acts as an inverse L-SR1 approximation,
  %   i.e., its inverse is an approximation of the Hessian. It is used
  %   as follows:
  %
  %   B = opLSR1(n, mem);
  %   B = update(B, s, y);
  %   d = - B \ g;          % Apply inverse L-SR1.
  %
  %   The operator may also be used in forward mode, i.e., as an
  %   approximation to of the Hessian. In this case, the attribute
  %   update_forward should be set to true, as forward mode incurs
  %   additional computational cost. It is used as follows:
  %
  %   B = opLSR1(n, mem);
  %   B.update_forward = true;
  %   B = update(B, s, y);
  %   d = - B \ g;          % Apply inverse L-SR1.
  %   Bx = B * x;           % Apply forward L-SR1.

  %   D. Orban, 2016.

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
    insert;  % Current insertion point.
    a;       % Storage of limited-memory terms.
    as;
    gamma;   % Scaling factor.
    need_refresh;  % Indicate when to update scaling and rank-1 terms.
  end

  properties (SetAccess = public)
    update_forward;  % Whether or not to update forward L-SR1.
    scaling;
    updates;         % number of update attempts
    rejects;         % number of rejected updates
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function op = opLSR1(n, mem)
      %opLSR1  Constructor.
      if nargin == 1
        mem = 1;
      end
      if nargin > 2
        error('At most two arguments can be specified.')
      end

      % Check if input is an integer
      if ~(isnumeric(mem) || mem ~= round(mem))
        error('Memory parameter must be an integer.');
      end

      % Create object
      op = op@opSpot('L-SR1', n, n);
      op.cflag  = false;
      op.sweepflag  = true;
      op.mem = max(mem, 1);
      op.s = zeros(n, op.mem);
      op.y = zeros(n, op.mem);
      op.ys = zeros(op.mem, 1);
      op.a = zeros(op.n, op.mem);
      op.as = zeros(op.mem, 1);
      op.update_forward = true;
      op.insert = 1;
      op.scaling = false;
      op.gamma = 1;
      op.need_refresh = false;
      op.updates = 0;
      op.rejects = 0;
    end % function opLSR1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function op = set.update_forward(op, val)
      op.update_forward = val;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Must use H = update(H, s, y)...
    % How do you get the syntax H.update(s,y) to work???

    function op = update(op, s, y)
      %store  Store the new pair {s,y} into the L-SR1 approximation.
      %       Discard oldest pair if memory has been exceeded.

      op.updates = op.updates + 1;
      Bs = op * s;
      ymBs = y - Bs;
      ys = dot(y, s);
      upTol = 1.0e-8;
      well_defined = abs(dot(ymBs, s)) >= upTol + upTol * norm(s) * norm(ymBs);

      sufficient_curvature = true;
      scaling_condition = true;
      y_neq_s = true;
      if op.scaling
        sufficient_curvature = abs(ys) >= upTol;
        if sufficient_curvature
          scaling_factor = ys / dot(y, y);
          scaling_condition = norm(y - s / scaling_factor) >= upTol;
        end
      else
        y_neq_s = norm(y - s) >= 1.0e-8;
      end

      if ~(well_defined && sufficient_curvature && scaling_condition && y_neq_s)
        msg = 'L-SR1: Rejecting {s,y} pair';
        if ~well_defined
          msg = strcat(msg, ' (not well defined)');
        elseif ~sufficient_curvature
          msg = strcat(msg, ' (curvature)');
        elseif ~scaling_condition
          msg = strcat(msg, ' (scaling)');
        else
          msg = strcat(msg, ' (y=s)');
        end
        warning(msg);
        op.rejects = op.rejects + 1;
      else

        op.s(:, op.insert) = s;
        op.y(:, op.insert) = y;
        op.ys(op.insert) = ys;

        % Update next insertion position.
        op.insert = mod(op.insert, op.mem) + 1;

        % Need to refresh scaling factor and rank-1 terms before next product
        op.need_refresh = true;
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function x = double(op)
      %double  Convert operator to a double.

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
        error('L-SR1: not using forward mode. Set update_forward = true.');
      end

      if op.need_refresh
        op = op.refresh();
      end

      d = ones(op.n, 1);
      if op.scaling
        d = d / op.gamma;
      end

      for i = 1 : op.mem
        k = mod(op.insert + i - 2, op.mem) + 1;
        if op.ys(k) ~= 0
          for j = 1 : op.n
            d(j) = d(j) + op.a(j, k)^2 / op.as(k);
          end
        end
      end
    end % function diag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end % Methods

  methods ( Access = protected )

    function op = refresh(op)
      %refresh  Refresh internal data before computing an operator-vector product

      % Update scaling factor.
      if op.scaling
        last = mod(op.insert - 2, op.mem) + 1;
        if op.ys(last) ~= 0
          op.gamma = op.ys(last) / dot(op.y(:, last), op.y(:, last));
        end
      end

      % Update rank-1 terms.
      for i = 1 : op.mem
        k = mod(op.insert + i - 2, op.mem) + 1;
        if op.ys(k) ~= 0
          op.a(:, k) = op.y(:, k) - op.s(:, k) / op.gamma;
          for j = 1 : i-1
            l = mod(op.insert + j - 2, op.mem) + 1;
            if op.ys(l) ~= 0
              op.a(:, k) = op.a(:, k) - dot(op.a(:, l), op.s(:, k)) / op.as(l) * op.a(:, l);
            end
          end
          op.as(k) = dot(op.a(:, k), op.s(:, k));
        end
      end

      op.need_refresh = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function q = multiply(op, x, ~)
      %multiply  Multiply operator with a vector.

      if ~op.update_forward
        error('L-SR1: not using forward mode. Set update_forward = true.');
      end

      if op.need_refresh
        op = op.refresh();
      end

      q = x;

      if op.scaling
        q = q / op.gamma;
      end

      for i = 1 : op.mem
        k = mod(op.insert + i - 2, op.mem) + 1;
        if op.ys(k) ~= 0
          q = q + dot(op.a(:, k), x) / op.as(k) * op.a(:, k);
        end
      end

    end % function multiply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function r = divide(op, b, mode)
    %   %divide  Solve a linear system with the operator.
    %
    %
    % end % function divide
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end % methods

end % Classdef

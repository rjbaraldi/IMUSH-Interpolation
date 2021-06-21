function d = diag(A)
%DIAG  Diagonal operator and diagonals of an operator.
%
%   diag(OP) is the main diagonal of the Spot operator OP.
%
%   See also opDiag.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

   if ismethod(A, 'diagonal')
     % If an operator type implements its own efficient diag(), use it.
     d = A.diagonal();

   else

      [m,n] = size(A);
      k = min(m,n);
      d = zeros(k,1);
      v = zeros(n,1);
      for i=1:k
         v(i) = 1;
         w = A*v;
         d(i) = w(i);
         v(i) = 0;
      end
   end

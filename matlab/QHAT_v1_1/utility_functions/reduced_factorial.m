function [E] = reduced_factorial(A,B)
% reduced_factorial(A,B) - calculates the ratio of products, cancelling all
%   available integer factors beforehand to reduce overflow errors. A and B
%   are one-dimensional matrices containing the factors in the numerator
%   and denominator, respectively. Useful when calculating Clebsch-Gordan
%   coefficients or Wigner 6j symbols.
% 
% 
% Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
% at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
% reachable at mason47@llnl.gov.
% 
% CODE-609912. All rights reserved.
% 
% This file is part of the Quaternionic Harmonic Analysis of Texture.  
% Please read LICENSE.txt for Our Notice and GNU General Public License
% information.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License (as published by
% the Free Software Foundation) version 2, dated June 1991.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
% conditions of the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% Separate nonintegers
C = A(mod(A,1)>1e-10);
D = B(mod(B,1)>1e-10);

% Sort integer factors in increasing order
A = int32(A(mod(A,1)<1e-10));
B = int32(B(mod(B,1)<1e-10));
A = sort(A);
B = sort(B);

% Storage for integer factors after cancellation
e = zeros(1,length(A));
f = zeros(1,length(B));

c = 1;
d = 1;
g = 1;
h = 1;
while c <= length(A)
    if A(c) == B(d)
        c = c+1;
        d = d+1;
    else
        if A(c) < B(d)
            e(g) = A(c);
            g = g+1;
            c = c+1;
        else
            if A(c) > B(d)
                f(h) = B(d);
                h = h+1;
                d = d+1;
            end
        end
    end
    if d > length(B)
        while c <= length(A)
            e(g) = A(c);
            g = g+1;
            c = c+1;
        end
    end
    if c > length(A)
        while d <= length(B)
            f(h) = B(d);
            h = h+1;
            d = d+1;
        end
    end
end

% Sort remaining factors in increasing order
A = [e(1:g-1),C];
B = [f(1:h-1),D];
A = sort(A);
B = sort(B);

% Multiply out factors starting with smallest
E = 1;
if length(A) > length(B)
    for p = 1:length(B)
        E = E*A(length(A)-p+1)/B(length(B)-p+1);
    end
    E = E*prod(A(1:length(A)-length(B)));
else
    for p = 1:length(A)
        E = E*A(length(A)-p+1)/B(length(B)-p+1);
    end
    E = E/prod(B(1:length(B)-length(A)));
end

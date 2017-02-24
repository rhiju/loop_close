function [U] = complexreal(N)
% complexreal(N) - constructs a sparse matrix that performs the similarity
%   transformation to convert from the complex harmonics to the real
%   harmonics. Rows ordered by increasing l, and by decreasing m for a 
%   given l. Columns ordered by increasing l, by decreasing m for a given 
%   l, and by c before s for a given m. N is the maximum value of l.
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

U = sparse((N+1)^2,(N+1)^2);
for L = 0:N
    U((L+1)^2-L,(L+1)^2) = (1i)^L;
    for M = 1:L
        U((L+1)^2-L-M,(L+1)^2-2*M) = (-1)^M*(1i)^L/sqrt(2);
        U((L+1)^2-L+M,(L+1)^2-2*M) = (1i)^L/sqrt(2);
        U((L+1)^2-L-M,(L+1)^2-2*M+1) = (-1)^M*(1i)^(L-1)/sqrt(2);
        U((L+1)^2-L+M,(L+1)^2-2*M+1) = -(1i)^(L-1)/sqrt(2);
    end
end

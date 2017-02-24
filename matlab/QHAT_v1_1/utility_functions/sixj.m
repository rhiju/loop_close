function [A] = sixj(a,b,c,d,e,f)
% sixj(a,b,c,d,e,f) - calculates the Wigner 6j symbol with the triad a b c
%   above d e f, in the standard notation. Does not check for vanishing of
%   symbol.
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

zmax = min([2*b,b+c-e+f,b+c+e+f+1,-a+b+c,b-d+f,a+b+c+1,b+d+f+1]);
zmin = 0;
kmax = zmax-zmin;

A = 1;
for k = kmax-(0:(kmax-1))
    A = 1-((-a+b+c-k+1)*(b-d+f-k+1)*(a+b+c+1-k+1)*(b+d+f+1-k+1))/(k*(2*b-k+1)*(b+c-e+f-k+1)*(b+c+e+f+1-k+1))*A;
end

A = A*reduced_factorial([1:(a-b+c),1:(-c+d+e),1:(a+e-f),1:(-b+d+f),1:(2*b),1:(2*b),1:(b+c-e+f),1:(b+c-e+f),1:(b+c+e+f+1),1:(b+c+e+f+1)],[1:(a+b-c),1:(-a+b+c),1:(c-d+e),1:(c+d-e),1:(a-e+f),1:(-a+e+f),1:(b+d-f),1:(b-d+f),1:(a+b+c+1),1:(c+d+e+1),1:(a+e+f+1),1:(b+d+f+1)]);
A = A*(-1)^(b+c+e+f);

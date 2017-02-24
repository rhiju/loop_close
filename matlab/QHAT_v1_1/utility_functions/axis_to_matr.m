function [R] = axis_to_matr(x,a,j)
% axis_to_matr(x,a,j) - returns the rotation matrix of dimension (2j+1) for
%   the rotation of angle a around axis x. The rotation axis is given in
%   Cartesian coordinates. The rows and colums are ordered by decreasing m.
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

x = x/sqrt(x*x');
q = cos(a/2);
Q = sin(a/2).*x;

R = zeros(1,2*j+1);

a = q-1i*Q(3);
b = -Q(2)-1i*Q(1);

aa = q^2+Q(3)^2;
bb = Q(2)^2+Q(1)^2;

R(1,1) = a^(2*j);
R(2*j+1,2*j+1) = b^(2*j);
for mup = j-(0:2*j-1)
    rmin = max([0,mup]);
    rmax = min([j,j+mup]);
    c = j+mup-rmin;
    d = j-rmin;
    e = -mup+rmin;
    f = rmin;
    tau = 1;
    for n = rmax-(rmin:(rmax-1))
        tau = 1-((c-n+1)*(d-n+1)*bb)/((e+n)*(f+n)*aa)*tau;
    end
    R(j+1,j-mup+1) = (-1)^rmin*(a^c*a'^d*b^e*b'^f)*factalt({1:j,sqrt(1:(j+mup)),sqrt(1:(j-mup))},{1:c,1:d,1:e,1:f})*tau;
end

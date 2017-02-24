function [R] = rotation(ax,an,N)
% rotation(ax,an,N) - constructs irreducible representation of SO(3) of 
%   dimension (N+1) for the rotation specified by axis an and angle an. 
%   Rows and columns ordered by increasing l, and by decreasing m for a 
%   given value of l.
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

r = 5;

j = N/2;
R = zeros(N+1,N+1);

% Construct quaternions for left and right isoclinic rotations
ax = ax/sqrt(ax*ax');
q = cos(an/(2*r));
Q = sin(an/(2*r)).*ax;

a = q-1i*Q(3);
lb = -Q(2)-1i*Q(1);

% Rotation about z axis
if abs(lb) < 1e-10
    du = (q+1i*Q(3))/(q-1i*Q(3));
    R(1,1) = a^N;
    for d = 2:N+1
        R(d,d) = R(d-1,d-1)*du;
    end
    R = R^r;
    return;
end

% Rotation about axis in x-y plane
if abs(a) < 1e-10
    sd = (Q(2)-1i*Q(1))/(-Q(2)-1i*Q(1));
    R(1,N+1) = lb^N;
    for d = 2:N+1
        R(d,N-d+2) = R(d-1,N-d+3)*sd;
    end
    R = R^r;
    return;
end

% General rotation
if abs(a) > 1e-10 && abs(lb) > 1e-10
    aa = q^2+Q(3)^2;
    bb = Q(2)^2+Q(1)^2;

    R(1,1) = a^N;
    R(N+1,N+1) = lb^N;
    for mup = j-(1:2*j)
        R(1,j-mup+1) = realsqrt(prod(1:2*j)/(prod(1:(j+mup))*prod(1:(j-mup))))*a^(j+mup)*lb^(j-mup);
    end
    for mpr = j-(1:2*j)
        for mup = j-(0:2*j-1)
            rmin = max([0,mup-mpr]);
            rmax = min([j-mpr,j+mup]);
            c = j+mup-rmin;
            d = j-mpr-rmin;
            e = mpr-mup+rmin;
            f = rmin;
            ltau = 1;
            for n = rmax-(rmin:(rmax-1))
                ltau = 1-((c-n+1)*(d-n+1)*bb)/((e+n)*(f+n)*aa)*ltau;
            end
            R(j-mpr+1,j-mup+1) = realsqrt(prod(1:(j+mpr))*prod(1:(j-mpr))*prod(1:(j+mup))*prod(1:(j-mup)))*(-1)^rmin*(a^c*a'^d*lb^e*lb'^f)/(prod(1:c)*prod(1:d)*prod(1:e)*prod(1:f))*ltau;
        end
        R(j-mpr+1,N+1) = realsqrt(prod(1:2*j)/(prod(1:(j+mpr))*prod(1:(j-mpr))))*a'^(j-mpr)*lb^(j+mpr);
    end
    
    % Force matrix to be orthogonal
    R = R^r;
    [U,~,V] = svd(R);
    R = U*V';
    return;
end

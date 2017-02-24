% Construction of the symmetrized harmonics corresponding to the
%   eigenvectors calculated by make_symm.m or make_symm_loop.m. The result
%   is the function Z, stored in the file x_y.mat, where x is the
%   principal order N and y is the number of the harmonic for a given N. Z
%   is a three dimensional array where the first index corresponds to c, 
%   the second index corresponds to b, and the third index corresponds to
%   a. Notice that the density of points is different for the three
%   dimensions.
% 
% N_max - maximum principal order of the harmonics to be constructed
% pts - resolution of the numerical representation
% 
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% P, p - associated Legendre polynomials
% G, g - normalized Gegenbauer polynomials
% weight - nonzero coefficients of the real harmonics to be included in the
%   expansion of the symmetrized harmonic
% lindex, mindex, iindex - l, m and i of the real harmonics with nonzero
%   coefficients to be included in the expansion of the symmetrized 
%   harmonic. i = 0 corresponds to c, i = 1 corresponds to s.
% Z - explicit evaluation of the symmetrized harmonicm defined by the
%   expansion coefficients in eigenvectors
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

N_max = 32;

pts = 200;

c = zeros(pts,pts/2,pts/2);
p = zeros(pts,pts/2,pts/2);
g = zeros(pts,pts/2,pts/2);

A = linspace(0,pi,pts/2);
B = linspace(0,pi,pts/2);
C = linspace(0,2*pi,pts);
for j = 1:pts
    c(j,:,:) = C(j);
end

for N = 0:2:N_max
    try
        load(['../eigenvectors/',num2str(N),'.mat']);
    catch
        continue;
    end
    disp(num2str(N));
    for q = 1:size(sym,2)-3
        indices = (abs(sym(:,q+3)) > 1e-8);

        weight = sym(indices,q+3)';
        lindex = sym(indices,1)';
        mindex = sym(indices,2)';
        iindex = sym(indices,3)';

        Z = zeros(pts,pts/2,pts/2);
        for r = 1:size(weight,2)
            L = lindex(r);
            M = mindex(r);

            G = zeros(1,pts/2);
            for s = 0:N-L
                G = G+((prod((L+1):(L+s))*prod((N-L-s+1):(N-s)))/prod(1:s))*cos((2*s-N+L)*A);
            end
            G = 2^(L+1/2)*realsqrt((N+1)/(pi*prod((N-L+1):(N+L+1))))*sin(A).^L.*G;
            P = legendre(L,cos(B));
            for s = 1:pts/2;
                g(:,:,s) = G(s);
                p(:,s,:) = P(M+1,s);
            end

            if M > 0.5
                if iindex(r) < 0.5
                    Z = Z+weight(r)*(-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*g.*p.*cos(M*c);
                else
                    Z = Z+weight(r)*(-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*g.*p.*sin(M*c);
                end
            else
                Z = Z+weight(r)*(-1)^L*realsqrt((2*L+1)/(4*pi))*g.*p;
            end
        end
        save([num2str(N),'_',num2str(q)],'Z');
    end
end

function [rhsh] = real_coeffs(T,N_max)
% real_coeffs(T,N_max) - calculates the coefficients in the expansion of
%   the discrete texture represented by quaternions in the struct T over
%   the real harmonics. Returns a cell array rhsh of coefficients, where 
%   the index of the cell corresponds to N/2+1 and the coefficients within
%   the cell are ordered by increasing L, decreasing M for a given L, and 
%   with C before S for a given M.
%
% N_max - maximum principal order of the real harmonics used
% a_q, b_q, c_q - triplets of spherical angles corresponding to the 
%   quaternions in the struct T
% G - normalized Gegenbauer polynomials
% P - associated Legendre polynomials
% itgrl - the integral of the product of a real harmonic with the
%   summation of delta functions corresponding to the orientations in T
% rhsh - a cell array storing the coefficients of the expansion. The index
%   of the cell corresponds to N/2+1 and the coefficients within the cell
%   are ordered by increasing L, decreasing M for a given L, and with C
%   before S for a given M.
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

assert(mod(N_max,2)==0,'N_max must be an even integer');

T = cell2mat(T);

a_q = real(acos(T(:,1)));
b_q = real(acos(T(:,4)./sin(a_q)));
b_q(abs(sin(a_q))<1e-8) = 0;
c_q = real(acos(T(:,2)./(sin(a_q).*sin(b_q))));
c_q(T(:,3)./(sin(a_q).*sin(b_q))<0) = 2*pi-c_q(T(:,3)./(sin(a_q).*sin(b_q))<0);
c_q(abs(sin(a_q).*sin(b_q))<1e-8) = 0;

rhsh = cell(N_max/2+1,1);

for N = 0:2:N_max
    disp(num2str(N));
    rhsh{N/2+1} = zeros((N+1)^2,1);
    for L = 0:N
        for M = 0:L
            G = zeros(size(T,1),1);
            for s = 0:N-L
                G = G+((prod((L+1):(L+s))*prod((N-L-s+1):(N-s)))/prod(1:s))*cos((2*s-N+L)*a_q);
            end
            G = 2^(L+1/2)*realsqrt((N+1)/(pi*prod((N-L+1):(N+L+1))))*sin(a_q).^L.*G;

            P = legendre(L,cos(b_q));
            P = P(M+1,:)';
            
            if M > 0.5
                rhsh{N/2+1}((L+1)^2-2*M) = sum((-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*G.*P.*cos(M*c_q))/size(T,1);
                rhsh{N/2+1}((L+1)^2-2*M+1) = sum((-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*G.*P.*sin(M*c_q))/size(T,1);
            else
                rhsh{N/2+1}((L+1)^2) = sum((-1)^L*realsqrt((2*L+1)/(4*pi))*G.*P)/size(T,1);
            end
        end
    end
end

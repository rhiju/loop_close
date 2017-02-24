function [CG] = clebschgordan(N)
% clebschgordan(N) - constructs a sparse form of the Clebsch-Gordan
%   coefficient matrix.  Let j_1 and j_2 be the uncoupled total angular
%   momenta, m_1 and m_2 be the uncoupled z components, j be the coupled 
%   total angular momentum and m be the coupled z component. Notice that
%   the properties of the expansion constrain the relevant Clebsch-Gordan
%   coefficients to satisfy j_1 = j_2. The columns of the matrix are
%   ordered by increasing j, and decreasing m for a given j.  The rows of
%   the matrix are ordered by decreasing m_1, and decreasing m_2 for a
%   given m_1. N is the maximum value of j.
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

j = N/2;

CG = sparse((N+1)^2,(N+1)^2);
rma = (j-(1:N+1)'+1);
lma = (j-(1:N+1)+1);
mar = rma*ones(1,N+1)+ones(N+1,1)*lma;

CG(1,N^2+1) = 1;
CG((N+1)^2,(N+1)^2) = 1;
for M = N-(1:2*N-1)
    cm = prod(1:N)*realsqrt((prod(1:(N+M))*prod(1:(N-M)))/prod(1:2*N));
    [rd,ld] = find(abs(mar-M)<1e-10);
    for ar = 1:size(rd,1)
        CG((N+1)*(j-lma(ld(ar)))+(j-rma(rd(ar))+1),(N+1)^2-N-M) = 1/(prod(1:(N-j+lma(ld(ar))))*prod(1:(N-j-rma(rd(ar)))))*realsqrt((prod(1:(j+lma(ld(ar))))*prod(1:(j-rma(rd(ar)))))/(prod(1:(j-lma(ld(ar))))*prod(1:(j+rma(rd(ar))))))*cm;
    end
end
lma = [nan,j-(1:N)];
mar = rma*ones(1,N+1)+ones(N+1,1)*lma;
for L = 0:N-1
    CG(N-L+1,L^2+1) = realsqrt(((2*L+1)*prod((L+1):(2*L)))/prod((N+1):(N+L+1)));
    for rm = L-j+1:j
        CG((N+1)*(j-L+rm)+j-rm+1,L^2+1) = CG((N+1)*(j-L+rm-1)+j-rm+2,L^2+1)*(-1)*realsqrt(((j+rm)*(j-rm+1))/((j+L-rm+1)*(j-L+rm)));
    end
    for M = L-(1:L)
        CG(N-M+1,(L+1)^2-L-M) = realsqrt(((2*L+1)*prod((N-L+1):(N-M))*prod((L-M+1):(L+M)))/(prod((N+1):(N+L+1))*prod(1:M)));
    end
    for M = L-(1:2*L)
        cm = prod(1:L)*realsqrt(((2*L+1)*prod(1:(L+M))*prod(1:(L-M)))/prod((N-L+1):(N+L+1)));
        [rd,ld] = find(abs(mar-M)<1e-10);
        for ar = 1:size(rd,1)
            rmin = max([0,-L+j-lma(ld(ar)),-L+j+rma(rd(ar))]);
            rmax = min([N-L,j-lma(ld(ar)),j+rma(rd(ar))]);
            a = rmin;
            b = N-L-rmin;
            c = j-lma(ld(ar))-rmin;
            d = j+rma(rd(ar))-rmin;
            e = L-j+lma(ld(ar))+rmin;
            f = L-j-rma(rd(ar))+rmin;
            tau = 1;
            for n = rmax-(rmin:(rmax-1))
                tau = 1-((b-n+1)*(c-n+1)*(d-n+1))/((a+n)*(e+n)*(f+n))*tau;
            end
            CG((N+1)*(j-lma(ld(ar)))+(j-rma(rd(ar))+1),(L+1)^2-L-M) = tau*realsqrt(prod(1:(j+lma(ld(ar))))*prod(1:(j-lma(ld(ar))))*prod(1:(j+rma(rd(ar))))*prod(1:(j-rma(rd(ar)))))*(-1)^rmin/(prod(1:a)*prod(1:b)*prod(1:c)*prod(1:d)*prod(1:e)*prod(1:f))*cm;
        end
    end
end

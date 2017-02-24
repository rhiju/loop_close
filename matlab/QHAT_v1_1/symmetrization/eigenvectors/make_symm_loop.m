% Construction of the coefficients in the real expansion of the symmetrized
%   harmonics for a set of values of N. Allows the use of four symmetry
%   group generators, represented as paired left and right isoclinic 
%   rotations. Set the left isoclinic rotation to a symmetry operation and
%   the right isoclinic rotation to the identity for a sample symmetry, or
%   the right isoclinic rotation to a symmetry operation and the left
%   isoclinic rotation to the identity for a crystal symmetry. Coefficients
%   are stored in sym with rows ordered by increasing l, decreasing m for a 
%   given l, and with c before s for a given m. sym is saved in the file
%   num2str(N).mat
%   
% N_max - maximum principal order of the harmonics to be symmetrized
% max_error - magnitude of acceptable rounding error
% exchange_symm - when true, enforces the grain exchange symmetry operation
%   required for misorientations of homophase boundaries
% lxa, lxb, lxc, lxd - axes of the left isoclinic rotations of symmetry
%   group generators a, b, c, and d
% lna, lnb, lnc, lnd - rotation angles of the left isoclinic rotations of 
%   symmetry group generators a, b, c, and d
% rxa, rxb, rxc, rxd - axes of the right isoclinic rotations of symmetry
%   group generators a, b, c, and d
% rna, rnb, rnc, rnd - rotation angles of the right isoclinic rotations of 
%   symmetry group generators a, b, c, and d
% 
% Ra, Rb, Rc, Rd - (N+1)^2 dimensional irreducible representations of SO(4)
%   corresponding to the symmetry group generators a, b, c, and d. Given in
%   the uncoupled basis by the Kronecker product of irreducible 
%   representations of SO(3) corresponding to left and right isoclinic
%   rotations. Rows and colums ordered by decreasing m_1, and decreasing
%   m_2 for a given value of m_1. This basis reduces numerical overhead and
%   rounding error.
% S - simultaneous eigenvectors of eigenvalue one for all of the symmetry
%   group generators considered. Rows initially ordered by deacreasing m_1,
%   and decreasing m_2 for a given value of m_1. Eventually transformed to
%   the real coupled basis with rows ordered by increasing l, decreasing m
%   for a given l, and with c before s for a given m.
% sym - matrix containing coefficients of the symmetrized harmonics and the
%   corresponding indices. Columns one and two give values of l and m.
%   Column three contains 0 for c and 1 for s.
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

max_error = 1e-10;

exchange_symm = false;

% Symmetry generators
lxa = [0 0 1];
lna = 0;
rxa = [0 0 1];
rna = pi/2;

lxb = [1 0 0];
lnb = 0;
rxb = [1 0 0];
rnb = pi/2;

lxc = [0 0 1];
lnc = pi/180;
rxc = [0 0 1];
rnc = 0;

lxd = [1 0 0];
lnd = 0;
rxd = [1 0 0];
rnd = 0;

for N = 0:2:N_max
    j = N/2;
    
    transform = full(clebschgordan(N)*complexreal(N));

    % Construct simultaneous eigenvectors
    Ra = kron(rotation(rxa,rna,N)',rotation(lxa,lna,N));
    [va,d] = eig(Ra);
    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
    if any(col)
        S = orth(va(:,col));

        Rb = kron(rotation(rxb,rnb,N)',rotation(lxb,lnb,N));
        [vb,d] = eig(S'*Rb*S);
        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
        if any(col)
            S = S*orth(vb(:,col));

            Rc = kron(rotation(rxc,rnc,N)',rotation(lxc,lnc,N));
            [vc,d] = eig(S'*Rc*S);
            col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            if any(col)
                S = S*orth(vc(:,col));

                Rd = kron(rotation(rxd,rnd,N)',rotation(lxd,lnd,N));
                [vd,d] = eig(S'*Rd*S);
                col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                if any(col)
                    S = S*orth(vd(:,col));

                    % Enforces the grain exchange symmetry
                    if exchange_symm
                        Ri = ones((N+1)^2,1);
                        for L = 0:N
                            Ri((L+1)^2-2*L:(L+1)^2) = Ri((L+1)^2-2*L:(L+1)^2)*(-1)^L;
                        end
                        Ri = transform*diag(Ri)*transform';
                        [vi,d] = eig(S'*Ri*S);
                        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                        if any(col)
                            S = S*orth(vi(:,col));
                        end
                    end
                end
            end
        end
    end

    if any(col)
        error = max([max(max(abs(Ra*S-S))),max(max(abs(Rb*S-S))),max(max(abs(Rc*S-S))),max(max(abs(Rd*S-S)))]);
        disp(['Rounding error: ',num2str(error)]);
        if error > max_error
            disp('Error threshold exceeded. Exiting.');
            break;
        end
        S = clean(transform'*S,N);
    else
        S = [];
    end

    % Assign indices
    sym = [zeros((N+1)^2,3) S];
    for L = 0:N
        sym((L+1)^2,1:3) = [L 0 0];
        for M = 1:L
            sym((L+1)^2-2*M,1:3) = [L M 0];
            sym((L+1)^2-2*M+1,1:3) = [L M 1];
        end
    end

    % Orders eigenvectors by indices of first nonzero coefficient
    for L = 0:N
        for M = 1:L
            sym((L+1)^2-2*M+(0:1),:) = flipud(sym((L+1)^2-2*M+(0:1),:));
        end
        sym((L+1)^2-2*L:(L+1)^2,:) = flipud(sym((L+1)^2-2*L:(L+1)^2,:));
    end
    order = zeros(1,size(sym,2)-3);
    for m = 1:size(sym,2)-3
        order(m) = find(abs(sym(:,m+3))>0,1);
    end
    [~,IX] = sort(order);
    sym(:,4:end) = sym(:,IX+3);
    for L = 0:N
        sym((L+1)^2-2*L:(L+1)^2,:) = flipud(sym((L+1)^2-2*L:(L+1)^2,:));
        for M = 1:L
            sym((L+1)^2-2*M+(0:1),:) = flipud(sym((L+1)^2-2*M+(0:1),:));
        end
    end
    
    if any(col)
        disp([num2str(size(sym,2)-3),' symmetrized harmonics for N = ',num2str(N)]);
        save(num2str(N),'sym');
    end
end

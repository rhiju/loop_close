function [S] = clean(S,N)
% clean(S,N) - attempts to minimize the number of nonzero coefficients of
%   the eigenvectors in S while maintaining orthonormality.
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

% Reorder rows for increasing l, and increasing m for given value of l
for L = 0:N
    S((L+1)^2-2*L:(L+1)^2,:) = flipud(S((L+1)^2-2*L:(L+1)^2,:));
end

% Forward sweep
for m = 1:size(S,2)
    mr = find(abs(S(:,m))>1e-10,1);
    S(:,m) = S(:,m)*conj(S(mr,m));
    S(:,m) = S(:,m)/sqrt(S(:,m)'*S(:,m));
    for n = m+1:size(S,2)
        nr = find(abs(S(:,n))>1e-10,1);
        if nr < mr
            swap = S(:,m);
            S(:,m) = S(:,n);
            S(:,n) = swap;
            
            mr = find(abs(S(:,m))>1e-10,1);
            S(:,m) = S(:,m)*conj(S(mr,m));
            S(:,m) = S(:,m)/sqrt(S(:,m)'*S(:,m));
        elseif nr == mr
            S(:,n) = S(:,n)-(S(mr,n)/S(mr,m))*S(:,m);
            S(:,n) = S(:,n)/sqrt(S(:,n)'*S(:,n));
        end
    end
end

% Backward sweep
for m = fliplr(1:size(S,2))
    S(:,m) = S(:,m)*conj(S(find(abs(S(:,m))>1e-10,1),m));
    S(:,m) = S(:,m)/sqrt(S(:,m)'*S(:,m));
    for n = fliplr(1:m-1)
        S(:,n) = S(:,n)-(S(:,m)'*S(:,n))*S(:,m);
        S(:,n) = S(:,n)/sqrt(S(:,n)'*S(:,n));
    end
end

% Attempt to remove remainders
Re = real(S);
Re(abs(Re)<1e-10) = 0;
Im = imag(S);
Im(abs(Im)<1e-10) = 0;
S = Re+1i*Im;

% First nonzero coefficient of the eigenvector is postive.
for m = 1:size(S,2)
    S(:,m) = sign(S(find(abs(S(:,m))>1e-10,1),m))*S(:,m);
end

% Reorder rows for increasing l, and decreasing m for given value of l
for L = 0:N
    S((L+1)^2-2*L:(L+1)^2,:) = flipud(S((L+1)^2-2*L:(L+1)^2,:));
end

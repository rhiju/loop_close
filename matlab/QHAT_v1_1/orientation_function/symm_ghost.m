function [shsh] = symm_ghost(shsh)
% symm_ghost(shsh) - modifies the coefficients of the expansion over the
%   symmetrized harmonics to remove regions of negative probability
%   density. Returns the modified coefficients in the cell array shsh,
%   where the index of the cell corresponds to N/2+1 and the coefficients
%   within the cell correspond to the symmetrized harmonics for the given
%   value of N.
%
% shsh - a cell array storing the coefficients of the expansion. The indiex
%   of the cell corresponds to N/2+1, while the index within the cell
%   corresponds to the symmetrized harmonic for the given value of N.
% N_max - maximum principal order of the symmetrized harmonics used
% pts - resolution of the symmetrized harmonics as generated by
%   symm_construct
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% dp - the invariant Haar measure for the space of normalized quaternions
% output - the values of the function expressed by the coefficients in shsh
% error - the portion of output for negative values, inverted in sign to be
%   positive
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

N_max = 2*(size(shsh,1)-1);

pts = 200;

a = zeros(pts,pts/2,pts/2);
b = zeros(pts,pts/2,pts/2);
c = zeros(pts,pts/2,pts/2);

A = linspace(0,pi,pts/2);
B = linspace(0,pi,pts/2);
C = linspace(0,2*pi,pts);
for j = 1:pts/2
    a(:,:,j) = A(j);
    b(:,j,:) = B(j);
end
for j = 1:pts
    c(j,:,:) = C(j);
end
dp = sin(a(1:end-1,1:end-1,1:end-1)+diff(A(1:2))/2).^2.*sin(b(1:end-1,1:end-1,1:end-1)+diff(B(1:2))/2)*diff(A(1:2))*diff(B(1:2))*diff(C(1:2));

output = zeros(pts,pts/2,pts/2);
for N = 0:2:N_max
    for q = 1:size(shsh{N/2+1},1)
        load(['../symmetrization/eigenfunctions/',num2str(N),'_',num2str(q),'.mat'],'Z');
        output = output+shsh{N/2+1}(q)*Z;
    end
end

out_min = min(output(output<0));
out_max = max(output(output>0));
disp(abs(out_max/out_min));
while ~isempty(out_min) && abs(out_max/out_min) < 10
    error = -output;
    error(error<0) = 0;
    
    scale = shsh{1}(1);
    for N = 0:2:N_max
        for q = 1:size(shsh{N/2+1},1)
            load(['../symmetrization/eigenfunctions/',num2str(N),'_',num2str(q),'.mat'],'Z');
            Y = error.*Z;
            Y = (Y(1:end-1,1:end-1,1:end-1)+Y(2:end,1:end-1,1:end-1)+Y(1:end-1,2:end,1:end-1)+Y(1:end-1,1:end-1,2:end)+Y(2:end,2:end,1:end-1)+Y(2:end,1:end-1,2:end)+Y(1:end-1,2:end,2:end)+Y(2:end,2:end,2:end))/8;
            shsh{N/2+1}(q) = shsh{N/2+1}(q)+sum(sum(sum(Y.*dp)));
        end
    end
    scale = scale/shsh{1}(1);
    for N = 0:2:N_max
        shsh{N/2+1} = shsh{N/2+1}*scale;
    end
    
    output = zeros(pts,pts/2,pts/2);
    for N = 0:2:N_max
        for q = 1:size(shsh{N/2+1},1)
            load(['../symmetrization/eigenfunctions/',num2str(N),'_',num2str(q),'.mat'],'Z');
            output = output+shsh{N/2+1}(q)*Z;
        end
    end
    
    out_min = min(output(output<0));
    out_max = max(output(output>0));
    disp(abs(out_max/out_min));
end

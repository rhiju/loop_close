% Displays the coefficients of the rth symmetrized harmonic of principal
%   index n, as calculated by make_symm.m or make_symm_loop.m
%
% N - order of the representation, or equivalently index N
% r - number of the symmetrized harmonic for the given value of N
%
% out - human readable format of the expansion of the given symmetrized
%   harmonic. Values within the parentheses correspond to l, m, and i,
%   respectively.
% weight - nonzero coefficients of the real harmonics included in the
%   expansion of the symmetrized harmonic
% lindex, mindex, iindex - l, m and i of the real harmonics with nonzero
%   coefficients included in the expansion of the symmetrized harmonic. 
%   i = 0 corresponds to c, i = 1 corresponds to s.
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

N = input('Principal index: ');
r = input('Harmonic number: ');

load(['eigenvectors/',num2str(N),'.mat']);

out = [];
weight = [];
lindex = [];
mindex = [];
iindex = [];

for L = 0:N
    sym((L+1)^2-2*L:(L+1)^2,:) = flipud(sym((L+1)^2-2*L:(L+1)^2,:));
end
for a = 1:size(sym,1)
    if abs(sym(a,r+3)) > 1e-8
        sn = [];
        wsn = [];
        if sym(a,r+3) < 0
            sn = '-';
            wsn = '-';
        elseif size(out,2) > 0
            sn = '+';
        end
        digits = num2str(abs(sym(a,r+3)),'%0.8f');
        if sym(a,3) < 0.5
            par = 'c';
        else
            par = 's';
        end
        out = [out,sn,digits,'(',num2str(sym(a,1)),',',num2str(sym(a,2)),',',par,')'];
        weight = [weight,wsn,digits,' '];
        lindex = [lindex,num2str(sym(a,1)),' '];
        mindex = [mindex,num2str(sym(a,2)),' '];
        iindex = [iindex,num2str(sym(a,3)),' '];
    end
end
for L = 0:N
    sym((L+1)^2-2*L:(L+1)^2,:) = flipud(sym((L+1)^2-2*L:(L+1)^2,:));
end

disp(' ');
disp('human readable output:')
disp(out);
disp('weight:');
disp(weight);
disp('lindex:');
disp(lindex);
disp('mindex:');
disp(mindex);
disp('iindex:');
disp(iindex);
disp(' ');

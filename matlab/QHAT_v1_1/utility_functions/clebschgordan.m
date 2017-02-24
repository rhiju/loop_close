function [C] = clebschgordan(lb,mb,la,ma,l,m)
% clebschgordan(lb,mb,la,ma,l,m) - calculates Clebsch-Gordan coefficient
%   for l m as the quantum numbers of the coupled momentum. Does not check
%   for vanishing of coefficient.
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

zmax = min([2*l,l+la+mb,l-lb+la,l+m,lb+la+l+1]);
zmin = 0;
kmax = zmax-zmin;

C = 1;
for k = kmax-(0:(kmax-1))
    C = 1-((l-lb+la-k+1)*(l+m-k+1)*(lb+la+l+1-k+1))/(k*(2*l-k+1)*(l+la+mb-k+1))*C;
end

C = C*realsqrt(reduced_factorial([1:(lb+la-l),1:(lb-mb),1:(2*l),1:(2*l),1:(l+la+mb),1:(l+la+mb)],[1:(lb+la+l+1),1:(lb-la+l),1:(l-lb+la),1:(lb+mb),1:(la+ma),1:(la-ma),1:(l+m),1:(l-m)]));
C = C*(-1)^(la+ma)*realsqrt(2*l+1);

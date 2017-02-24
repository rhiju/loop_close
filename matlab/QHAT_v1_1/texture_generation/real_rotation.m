function [R] = real_rotation(ax,an)
% real_rotation(ax,an) - constructs the canonical representation of SO(3) 
%   for the rotation specified by axis ax and angle an.
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

ax = ax/sqrt(ax*ax');
R = [1-2*(ax(2)^2+ax(3)^2)*sin(an/2)^2, -ax(3)*sin(an)+2*ax(1)*ax(2)*sin(an/2)^2, ax(2)*sin(an)+2*ax(3)*ax(1)*sin(an/2)^2; ax(3)*sin(an)+2*ax(1)*ax(2)*sin(an/2)^2, 1-2*(ax(3)^2+ax(1)^2)*sin(an/2)^2, -ax(1)*sin(an)+2*ax(2)*ax(3)*sin(an/2)^2; -ax(2)*sin(an)+2*ax(3)*ax(1)*sin(an/2)^2, ax(1)*sin(an)+2*ax(2)*ax(3)*sin(an/2)^2, 1-2*(ax(1)^2+ax(2)^2)*sin(an/2)^2];

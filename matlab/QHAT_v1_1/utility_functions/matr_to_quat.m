function [q,ang] = matr_to_quat(o)
% matr_to_quat(o) - returns the quaternion q and the angle and axis of
%   rotation ang for the rotation matrix o. The rotation axis is given in
%   spherical coordinates.
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

e = 1e-10;

% Quaternion construction
q = [0,0,0,0];
if sum(diag(o)) > 0
    q(1) = realsqrt(max([sum(diag(o))+1,0]))/2;
    q(2) = (o(3,2)-o(2,3))/(4*q(1));
    q(3) = (o(1,3)-o(3,1))/(4*q(1));
    q(4) = (o(2,1)-o(1,2))/(4*q(1));
else
    q(2) = realsqrt(max([1+o(1,1)-o(2,2)-o(3,3),0]))/2;
    q(3) = realsqrt(max([1-o(1,1)+o(2,2)-o(3,3),0]))/2;
    q(4) = realsqrt(max([1-o(1,1)-o(2,2)+o(3,3),0]))/2;
    [~,I] = max(q(2:4));
    switch I
        case 1
            q(1) = (o(3,2)-o(2,3))/(4*q(2));
            q(3) = (o(2,1)+o(1,2))/(4*q(2));
            q(4) = (o(1,3)+o(3,1))/(4*q(2));
        case 2
            q(1) = (o(1,3)-o(3,1))/(4*q(3));
            q(2) = (o(2,1)+o(1,2))/(4*q(3));
            q(4) = (o(3,2)+o(2,3))/(4*q(3));
        case 3
            q(1) = (o(2,1)-o(1,2))/(4*q(4));
            q(2) = (o(1,3)+o(3,1))/(4*q(4));
            q(3) = (o(3,2)+o(2,3))/(4*q(4));
    end
    if q(1) < 0
        q = -q;
    end
end

% Axis-angle construction
a = 2*acos(q(1));
if abs(sin(a/2)) > e
    b = acos(q(4)/sin(a/2));
    if abs(sin(a/2)*sin(b)) > e
        c = acos(q(2)/(sin(a/2)*sin(b)));
        if q(3)/(sin(a/2)*sin(b)) < 0
            c = 2*pi-c;
        end
    else
        c = 0;
    end
else
    b = 0;
    c = 0;
end
ang = [a,b,c];

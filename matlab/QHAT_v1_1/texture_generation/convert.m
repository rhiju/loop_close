function [T] = convert(C)
% convert(C) - converts the crystal orientations stored in C as rotation
%   matrices to crystal orientations stored in T as normalized quaternions
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

T = cell(size(C,1),1);

for r = 1:size(C,1)
    m = C{r};
    an = acos(1/2*(trace(m)-1));
    if abs(sin(an))>1e-10
        ax(1) = (m(3,2)-m(2,3))/(2*sin(an));
        ax(2) = (m(1,3)-m(3,1))/(2*sin(an));
        ax(3) = (m(2,1)-m(1,2))/(2*sin(an));
    else
        if abs(an)>1e-10
            q(1) = sqrt(1+m(1,1)-m(2,2)-m(3,3))/2;
            q(2) = sqrt(1+m(2,2)-m(3,3)-m(1,1))/2;
            q(3) = sqrt(1+m(3,3)-m(1,1)-m(2,2))/2;
            qnt = find(abs(q-max(q))<1e-10, 1);
            if abs(qnt-1)<1e-10
                ax(1) = q(1)/sin(an/2);
                ax(2) = (m(2,1)+m(1,2))/(4*q(1)*sin(an/2));
                ax(3) = (m(1,3)+m(3,1))/(4*q(1)*sin(an/2));
            end
            if abs(qnt-2)<1e-10
                ax(1) = (m(1,2)+m(2,1))/(4*q(2)*sin(an/2));
                ax(2) = q(2)/sin(an/2);
                ax(3) = (m(2,3)+m(3,2))/(4*q(2)*sin(an/2));
            end
            if abs(qnt-3)<1e-10
                ax(1) = (m(3,1)+m(1,3))/(4*q(3)*sin(an/2));
                ax(2) = (m(2,3)+m(3,2))/(4*q(3)*sin(an/2));
                ax(3) = q(3)/sin(an/2);
            end
        else
            ax = [1,1,1];
        end
    end
    ax = ax/sqrt(ax*ax');
    T{r} = [cos(an/2),sin(an/2).*ax];
end

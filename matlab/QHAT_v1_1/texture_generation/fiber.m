function [C] = fiber(OR,ax,noise,mult)
% fiber(OL,OR,ax,noise,mult) - constructs crystal orientations for a fiber 
%   textured material crystal symmetry OR, fiber axis ax, an angle 
%   deviation given by noise and a number of crystals given by mult. 
%   Returns crystal orientations as rotation matrices in the struct C.
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

C = {};
for q = 1:mult
	an = rand*2*pi;
    a = real_rotation(ax,an);
    c = {};
    for r = 1:length(OR)
        o = a*OR{r};
        add_o = true;
        for t = 1:size(c,1)
            if add_o && max(max(abs(o\c{t}-eye(3)))) < 1e-5
                add_o = false;
            end
        end
        if add_o
            c{end+1,1} = o;
        end
    end
    rax = randn(1,3);
    ran = randn*noise;
    nse = real_rotation(rax,ran);
    for r = 1:size(c,1)
        C{end+1,1} = nse*c{r};
    end
end

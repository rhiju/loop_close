function [C] = rolling(OL,OR,type,noise,mult)
% rolling(OL,OR,type,noise,mult) - constructs crystal orientations for a
%   rolled material with sample symmetry OL, crystal symmetry OR, an angle
%   deviation given by noise and a number of crystals given by mult. 
%   Texture is copper when type is true and brass when type is false.
%   Rolling direction corresponds to x and normal direction corresponds to
%   z. Returns crystal orientations as rotation matrices in the struct C.
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

permute = perms([1 2 3]);
sgn = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; -1 1 -1; 1 -1 -1; -1 -1 -1];

if type
    A = [1 1 2; 8 12 23; 18 24 51; 20 35 64; 1 4 6; 1 1 0];
    B = [1 1 1; 7 3 4; 3 2 2; 9 4 5; 2 1 1; 2 1 1];
    ratio = [15; 3; 3; 10; 3; 21];
else
    A = [1 1 0; 1 3 2; 1 4 6; 1 1 0];
    B = [1 1 2; 6 4 3; 2 1 1; 0 0 1];
    ratio = [15; 2; 2; 2];
end

mult = round(mult*ratio/sum(ratio));

C = {};

for q = 1:size(A,1)
    d = A(q,:);
    d = d/sqrt(d*d');
    a = [];
    for r = 1:size(permute,1)
        for s = 1:size(sgn,1)
            o = d(permute(r,:)).*sgn(s,:);
            add_o = true;
            for t = 1:size(a,1)
                if add_o && max(abs(a(t,:)-o)) < 1e-5
                    add_o = false;
                end
            end
            if add_o
                a(end+1,:) = o;
            end
        end
    end
    d = B(q,:);
    d = d/sqrt(d*d');
    b = [];
    for r = 1:size(permute,1)
        for s = 1:size(sgn,1)
            o = d(permute(r,:)).*sgn(s,:);
            add_o = true;
            for t = 1:size(b,1)
                if add_o && max(abs(b(t,:)-o)) < 1e-5
                    add_o = false;
                end
            end
            if add_o
                b(end+1,:) = o;
            end
        end
    end
    
    c = {};
    for r = 1:size(a,1)
        for s = 1:size(b,1)
            if abs(dot(b(s,:),a(r,:))) < 1e-5
                o = [b(s,:);cross(a(r,:),b(s,:));a(r,:)];
                add_o = true;
                for t = 1:size(c,1)
                    if add_o && max(max(abs(o\c{t}-eye(3)))) < 1e-5
                        add_o = false;
                    end
                end
                if add_o
                    for t = 1:size(OL,1)
                        for u = 1:size(OR,1)
                            o = OL{t}*[b(s,:);cross(a(r,:),b(s,:));a(r,:)]*OR{u};
                            add_o = true;
                            for v = 1:size(c,1)
                                if add_o && max(max(abs(o\c{v}-eye(3)))) < 1e-5
                                    add_o = false;
                                end
                            end
                            if add_o
                                c{end+1,1} = o;
                            end
                        end
                    end
                end
            end
        end
    end
    for r = 1:mult(q)
        rax = randn(1,3);
        ran = randn*noise;
        nse = real_rotation(rax,ran);
        for s = 1:size(c,1)
            C{end+1,1} = nse*c{s};
        end
    end
end

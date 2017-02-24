function pole_figures(C)
% pole_figures(C) - constructs pole figures for the texture with crystal
%   orientations stored in C. Plots the three figures for crystal
%   directions pa, pb, and pc. Currently configured for rolling textures,
%   though may be configured for other textures simply by changing the axis
%   labels.
%
% pa, pb, pc - crystal directions used in the construction of the pole
%   figures
% Pa, Pb, Pc - orientations of the corresponding crystal directions with
%   respect to the sample frame. Third coordinate is forced to be positive
%   since a stereographic projection is used.
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

% Crystal directions
pa = [1 0 0];
pb = [1 1 0];
pc = [1 1 1];

pa = pa/sqrt(pa*pa');
pb = pb/sqrt(pb*pb');
pc = pc/sqrt(pc*pc');

% Construct orientations of poles
Pa = zeros(size(C,1),3);
Pb = zeros(size(C,1),3);
Pc = zeros(size(C,1),3);
for q = 1:size(C,1)
    Pa(q,:) = (C{q}*pa')';
    Pb(q,:) = (C{q}*pb')';
    Pc(q,:) = (C{q}*pc')';
end
Pa(Pa(:,3)<0,:) = -Pa(Pa(:,3)<0,:);
Pb(Pb(:,3)<0,:) = -Pb(Pb(:,3)<0,:);
Pc(Pc(:,3)<0,:) = -Pc(Pc(:,3)<0,:);

% Reduction to remove redundant points
swap = [];
while size(Pa,1)>0
    swap(end+1,:) = Pa(1,:);
    Pa = Pa(~all(abs(Pa-ones(size(Pa,1),1)*Pa(1,:))<1e-5,2),:);
end
Pa = swap;

swap = [];
while size(Pb,1)>0
    swap(end+1,:) = Pb(1,:);
    Pb = Pb(~all(abs(Pb-ones(size(Pb,1),1)*Pb(1,:))<1e-5,2),:);
end
Pb = swap;

swap = [];
while size(Pc,1)>0
    swap(end+1,:) = Pc(1,:);
    Pc = Pc(~all(abs(Pc-ones(size(Pc,1),1)*Pc(1,:))<1e-5,2),:);
end
Pc = swap;

aa = [2*Pa(:,1)./(Pa(:,3)+1),2*Pa(:,2)./(Pa(:,3)+1),ones(size(Pa(:,3)))];
bb = [2*Pb(:,1)./(Pb(:,3)+1),2*Pb(:,2)./(Pb(:,3)+1),ones(size(Pb(:,3)))];
cc = [2*Pc(:,1)./(Pc(:,3)+1),2*Pc(:,2)./(Pc(:,3)+1),ones(size(Pc(:,3)))];

% Plotting in stereographic projection
figure(1), clf, hold on;
title('(100) Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',34);
scatter3(aa(:,1),aa(:,2),aa(:,3),30,[0 0 1],'filled');

figure(2), clf, hold on;
title('(110) Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',34);
scatter3(bb(:,1),bb(:,2),bb(:,3),30,[0 0 1],'filled');

figure(3), clf, hold on;
title('(111) Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',34);
scatter3(cc(:,1),cc(:,2),cc(:,3),30,[0 0 1],'filled');

x = linspace(-2,2,100);
for q = 1:3
    figure(q);
    plot(x,sqrt(4-x.^2),'k','linewidth',2);
    plot(x,-sqrt(4-x.^2),'k','linewidth',2);
    
    set(gcf,'position',[100,100,575,625]);
    axis([-2,2,-2,2,0,2]), axis equal, axis off;
    text('position',[0,0.2,2],'string','ND','fontname','Times New Roman','fontsize',28)
    text('position',[2.18,0.2,2],'string','RD','fontname','Times New Roman','fontsize',28)
    view(-90,90);
end

function inverse_pole_figures(C)
% inverse_pole_figures(C) - constructs inverse pole figures for the texture
%   with crystal orientations stored in C. Plots three complete figures
%   and three reduced figures for sample directions pa, pb, and pc.
%   Currently configured for rolling textures, though may be configured
%   for other textures simply by changing the titles.
%
% pa, pb, pc - sample directions used in the construction of the inverse 
%   pole figures
% Pa, Pb, Pc - orientations of the corresponding sample directions with
%   respect to the crystal frame. Third coordinate is forced to be
%   positive since a stereographic projection is used.
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

% Sample directions
pa = [1 0 0];
pb = [0 1 0];
pc = [0 0 1];

pa = pa/sqrt(pa*pa');
pb = pb/sqrt(pb*pb');
pc = pc/sqrt(pc*pc');

% Construct orientations of poles
Pa = zeros(size(C,1),3);
Pb = zeros(size(C,1),3);
Pc = zeros(size(C,1),3);
for q = 1:size(C,1)
    Pa(q,:) = (C{q}'*pa')';
    Pb(q,:) = (C{q}'*pb')';
    Pc(q,:) = (C{q}'*pc')';
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
title('RD Inverse Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',30);
scatter(aa(:,1),aa(:,2),30,[0 0 1],'filled');

figure(3), clf, hold on;
title('TD Inverse Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',30);
scatter(bb(:,1),bb(:,2),30,[0 0 1],'filled');

figure(5), clf, hold on;
title('ND Inverse Pole Figure','position',[-2.5,0,0],'fontname','Times New Roman','fontsize',30);
scatter(cc(:,1),cc(:,2),30,[0 0 1],'filled');

% Generation of guide lines
pts = 100;
b = linspace(-pi/2,pi/2,pts)';

e = zeros(pts,3,8);
e(:,:,1) = [sin(b) zeros(pts,1) cos(b)];
e(:,:,2) = (real_rotation([0 0 1],pi/4)*e(:,:,1)')';
e(:,:,3) = (real_rotation([0 0 1],pi/2)*e(:,:,1)')';
e(:,:,4) = (real_rotation([0 0 1],-pi/4)*e(:,:,1)')';
e(:,:,5) = (real_rotation([1 0 0],pi/4)*e(:,:,1)')';
e(:,:,6) = (real_rotation([1 0 0],-pi/4)*e(:,:,1)')';
e(:,:,7) = (real_rotation([0 1 0],pi/4)*e(:,:,3)')';
e(:,:,8) = (real_rotation([0 1 0],-pi/4)*e(:,:,3)')';
for q = 1:size(e,3)
    e(:,:,q) = [2*e(:,1,q)./(e(:,3,q)+1),2*e(:,2,q)./(e(:,3,q)+1),e(:,3,q)];
end
    
b = linspace(-2,2,100);
for q = 1:2:5
    figure(q);
    for r = 1:size(e,3)
        plot(e(:,1,r),e(:,2,r),'k','linewidth',2);
    end
    plot(b,sqrt(4-b.^2),'k','linewidth',2);
    plot(b,-sqrt(4-b.^2),'k','linewidth',2);
    
    set(gcf,'position',[100,100,575,625]);
    axis([-2,2,-2,2,0,2]), axis equal, axis off;
    text('position',[0.03,0.21,2],'string','001','fontname','Times New Roman','fontsize',30)
    text('position',[2.20,0.21,2],'string','100','fontname','Times New Roman','fontsize',30)
    text('position',[0.03,2.45,2],'string','010','fontname','Times New Roman','fontsize',30)
    view(-90,90);
end

% Plotting in stereographic projection
figure(2), clf, hold on;
title('RD Inverse Pole Figure','position',[-0.12,-0.42,0],'fontname','Times New Roman','fontsize',30);
filter = (aa(:,3)>=0 & aa(:,1)>=0 & aa(:,1)<=-aa(:,2) & aa(:,2)>2-sqrt(8-aa(:,1).^2));
scatter3(aa(filter,1),aa(filter,2),aa(filter,3),30,[0 0 1],'filled');

figure(4), clf, hold on;
title('TD Inverse Pole Figure','position',[-0.12,-0.42,0],'fontname','Times New Roman','fontsize',34);
filter = (bb(:,3)>=0 & bb(:,1)>=0 & bb(:,1)<=-bb(:,2) & bb(:,2)>2-sqrt(8-bb(:,1).^2));
scatter3(bb(filter,1),bb(filter,2),bb(filter,3),30,[0 0 1],'filled');

figure(6), clf, hold on;
title('ND Inverse Pole Figure','position',[-0.12,-0.42,0],'fontname','Times New Roman','fontsize',34);
filter = (cc(:,3)>=0 & cc(:,1)>=0 & cc(:,1)<=-cc(:,2) & cc(:,2)>2-sqrt(8-cc(:,1).^2));
scatter3(cc(filter,1),cc(filter,2),cc(filter,3),30,[0 0 1],'filled');

b = linspace(0,-1+sqrt(3),100);
for q = 2:2:6
    figure(q);
    plot(b,2-sqrt(8-b.^2),'k','linewidth',2);
    plot([-1+sqrt(3) 0 0],[1-sqrt(3) 0 2-2*sqrt(2)],'k','linewidth',2);
    
    set(gcf,'position',[100,100,575,625]);
    axis([-.1,.85,-.95,.1]), axis equal, axis off;
    text('position',[0,0.11,0],'string','001','fontname','Times New Roman','fontsize',28)
    text('position',[0,-.84,0],'string','101','fontname','Times New Roman','fontsize',28)
    text('position',[.78,-.68,0],'string','111','fontname','Times New Roman','fontsize',28)
    view(-90,90);
end

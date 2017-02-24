function euler(T)
% euler(T) - plots the crytal orientations stored as normalized quaternions
%   in the struct C using the Euler angle representation, in degrees. The
%   reduced Euler angle space for cubic crystal symmetry is shown, along 
%   with sections perpendicular to the \phi_2 axis.
%
% sections - number of sections perpendicular to the \phi_2 axis to plot
% T - crystal orientations represented as normalized quaternions
% e - triplets of Euler angles corresponding to crystal orientations, in 
%   degrees
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

sections = 19;

T = cell2mat(T);
e = zeros(size(T,1),3);
e(:,2) = acos(2*(T(:,1).^2+T(:,4).^2)-1);
e(:,1) = acos(2*(T(:,1).*T(:,2)-T(:,3).*T(:,4))./sin(e(:,2)));
e(2*(T(:,4).*T(:,2)+T(:,1).*T(:,3))./sin(e(:,2))<0,1) = 2*pi-e(2*(T(:,4).*T(:,2)+T(:,1).*T(:,3))./sin(e(:,2))<0,1);
e(:,3) = acos(2*(T(:,1).*T(:,2)+T(:,3).*T(:,4))./sin(e(:,2)));
e(2*(T(:,4).*T(:,2)-T(:,1).*T(:,3))./sin(e(:,2))<0,3) = 2*pi-e(2*(T(:,4).*T(:,2)-T(:,1).*T(:,3))./sin(e(:,2))<0,3);

figure(1), clf;
set(gcf,'position',[65,390,380,410]);
hold on, box on;

filter = all(e<=pi/2,2);
scatter3(e(filter,2),e(filter,1),e(filter,3),30,[0 0 1],'filled');

axis equal, axis([0 pi/2 0 pi/2 0 pi/2]);
set(gca,'zdir','reverse');
xlabel('\Phi'),ylabel('\phi_1'),zlabel('\phi_2');
view(140,30);

E = cell(sections,1);
e_edges = linspace(0,pi/2-eps,sections);

E{1} = e((e_edges(1)<=e(:,3) & e(:,3)<(e_edges(1)+e_edges(2))/2) | ((e_edges(end-1)+e_edges(end))/2<=e(:,3) & e(:,3)<e_edges(end)),:);
for q = 2:(sections-1)
    E{q} = e((e_edges(q-1)+e_edges(q))/2<=e(:,3) & e(:,3)<(e_edges(q)+e_edges(q+1))/2,:);
end
E{sections} = E{1};

for j = 1:sections
    figure(j+1), clf;
    set(gcf,'position',[20*j+45,410-20*j,380,410]);
    hold on, box on;
    
    plot([0,0,pi/2,pi/2,0],[0,pi/2,pi/2,0,0],'color',[0 0 0],'linewidth',1);
    scatter(E{j}(:,1),E{j}(:,2),30,[0 0 1],'filled');
    
    axis equal, axis([0,pi/2,0,pi/2]);
    set(gca,'xaxislocation','top','ydir','reverse');
    xlabel('\phi_1'),ylabel('\Phi');
    text(0.7,0.8,[num2str(e_edges(j)*180/pi),'\circ'],'fontname','Times New Roman','fontsize',30);
    view(0,90);
end

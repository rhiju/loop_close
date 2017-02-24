function scaled_axis(T,type)
% scaled_axis(T,type) - plots the crytal orientations stored as normalized
%   quaternions in the struct T using representations based on scaling the
%   rotation axis by a monotonic function of the rotation angle. Setting 
%   type to 'axis-angle' gives the full axis-angle space, setting type to
%   'rodrigues' gives the full Rodrigues vector space and individual cross
%   sections (for cubic crystal symmetry), and setting type to 'quaternion'
%   gives an orthographic projection of the normalized quaternions.
%
% pts - resolution of the numerical representation
% sections - number of surfaces of constant beta to plot
% T - crystal orientations represented as normalized quaternions
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

pts = 200;
sections = 9;

T = cell2mat(T);

figure(1), clf;
set(gcf,'position',[65,390,380,410]);
hold on, grid on, box on;

switch type
    case 'axis-angle'
        b = zeros(pts,pts);
        c = zeros(pts,pts);

        B = linspace(0,pi,pts);
        C = linspace(0,2*pi,pts);
        for j = 1:pts;
            b(:,j) = B(j);
            c(j,:) = C(j);
        end

        e3 = cos(b);
        e2 = sin(b).*cos(c);
        e1 = sin(b).*sin(c);
        
        n = diag(2*acos(T(:,1))./sqrt(sum(T(:,2:4).^2,2)))*T(:,2:4);
        scatter3(n(:,1),n(:,2),n(:,3),30,[0 0 1],'filled');
        surf(pi*e2,pi*e1,pi*e3,'edgealpha',.1,'facealpha',0)

        axis equal, axis([-pi pi -pi pi -pi pi]);
    case 'rodrigues'
        a = 3-2*realsqrt(2);
        b = realsqrt(2)-1;
        
        n = diag(tan(acos(T(:,1)))./sqrt(sum(T(:,2:4).^2,2)))*T(:,2:4);
        r = n(all(abs(n)<b,2) & abs(n(:,1)+n(:,2)+n(:,3))<1 & abs(n(:,1)+n(:,2)-n(:,3))<1 & abs(n(:,1)-n(:,2)+n(:,3))<1 & abs(-n(:,1)+n(:,2)+n(:,3))<1,:);
        scatter3(r(:,1),r(:,2),r(:,3),30,[0 0 1],'filled');
                
        c = [b b b b b b b b b]';
        d = [-a a b b a -a -b -b -a]';
        e = [b b a -a -b -b -a a b]';
        plot3(c,d,e,'k','linewidth',1);
        plot3(-c,d,e,'k','linewidth',1);
        plot3(e,d,c,'k','linewidth',1);
        plot3(e,d,-c,'k','linewidth',1);
        plot3(d,c,e,'k','linewidth',1);
        plot3(d,-c,e,'k','linewidth',1);
        
        axis equal, axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
        set(gca,'xtick',-0.5:0.5:0.5,'ytick',-0.5:0.5:0.5,'ztick',-0.5:0.5:0.5);
        
        R = cell(1,sections);
        z_edges = linspace(-b,b,sections+1);
        for n = 1:sections
            R{n} = r(z_edges(n)<=r(:,3) & r(:,3)<z_edges(n+1),:);
        end
        
        for m = 1:sections
            figure(m+1), clf;
            set(gcf,'position',[20*m+45,410-20*m,380,410]);
            hold on, box on;
            
            z = (z_edges(m)+z_edges(m+1))/2;
            d = 1-abs(z)-b;
            plot([b d -d -b -b -d d b b],[d b b d -d -b -b -d d],'color',[0 0 0],'linewidth',1);
            scatter(R{m}(:,1),R{m}(:,2),30,[0 0 1],'filled');
            
            e = b+0.005;
            axis equal, axis([-e e -e e]), axis off;
            text(-0.2,-b-0.06,['z=',num2str(z,'%4.3f')],'fontname','Times New Roman','fontsize',30);
        end
    case 'quaternion'
        b = zeros(pts,pts);
        c = zeros(pts,pts);

        B = linspace(0,pi,pts);
        C = linspace(0,2*pi,pts);
        for j = 1:pts;
            b(:,j) = B(j);
            c(j,:) = C(j);
        end

        e3 = cos(b);
        e2 = sin(b).*cos(c);
        e1 = sin(b).*sin(c);
        
        scatter3(T(:,2),T(:,3),T(:,4),30,[0 0 1],'filled');
        surf(e2,e1,e3,'edgealpha',.1,'facealpha',0)

        axis equal, axis([-1 1 -1 1 -1 1]);
    otherwise
        error([type,' is not a known scaled_axis parameterization.'])
end

figure(1);
xlabel('X'),ylabel('Y'),zlabel('Z');
view(140,30);

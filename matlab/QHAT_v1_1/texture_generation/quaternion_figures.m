function quaternion_figures(T,bounded,crystal_type,sectioning,plot_type)
% quaternion_figures(T,bounded,crystal_type,sectioning,plotting) - Various
%   projections of the crystal orientations stored as normalized
%   quaternions in the struct T. When plot_type is 'stereographic', the 
%   orientation space is sectioned into concentric spherical shells and 
%   projected stereographically into the plane, with theta restricted to be
%   less than pi/2. When plot_type is 'area_preserve', the orientation 
%   space is sectioned into concentric spherical sheels and an 
%   area-preserving projection is used. When plot_type is 'fundamental', 
%   Figure 1 gives the orthographic projection of the entire orientation
%   space and Figure 2 gives a rough representation of the boundary of the
%   fundamental zone.
% 
% bounded - constrains orientations to the fundamental zone when true. When
%   false, all orientations is shown.
% crystal_type - when 'cubic', fundamental zone is for cubic crystal and
%   orthorhombic sample symmetry. When 'orthorhomic', fundamental zone is
%   for orthorhomic crystal and orthorhomic sample symmetry.
% sectioning - when 'angle', sectioning of the fundamental zone is in equal
%   angle intervals. When 'volume', sectioning is by equal volume intervals
% plot_type - when 'stereographic', sections the fundamental zone into
%   concentric spherical shells and projects them stereographically into
%   the plane. When 'area_preserve', sections the fundamental zone into 
%   concentric spherical shells and uses an area-preserving projection into
%   the plane. When 'fundamental', the fundamental zone is projected
%   orthographically.
% pts - resolution of the numerical representation
% sections - number of surfaces of constant beta to plot
%
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% A_max - maximum value of beta in the fundamental zone. Should be
%   restricted to less than pi/2.
% A_edges - controls the sectioning type as described in the code. Stores
%   the values of beta bounding concentric spherical shells.
% e4, e3, e2, e1 - w, z, x, y coordinates, respectively
% filter - array of ones and nans used to restrict plotting of surfaces
%   to regions within the fundamental zone. Currently implemented for
%   orthorhomic and cubic crystal symmetries, as described in the code.
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
sections = 12;

T = cell2mat(T);

a = zeros(pts,pts/2,sections,3);
b = zeros(pts,pts/2,sections,3);
c = zeros(pts,pts/2,sections,3);

if bounded
    switch crystal_type
        case 'cubic'
            A_max = pi/5.733;
        case 'orthorhombic'
            A_max = pi/3;
        otherwise
            error([crystal_type,' is not a known crystal type.'])
    end
else
    A_max = pi/2;
end

% Selects sectioning by equal angle intervals when sectioning is 'angle' 
%   and by equal volume intervals when sectioning is 'volume'
switch sectioning
    case 'angle'
        A_edges = linspace(0,A_max,sections+1);
    case 'volume'
        A_edges = linspace(0,A_max-sin(A_max),sections+1);
        for j = 1:size(A_edges,2)
           A_edges(j) = fminsearch(@(x) abs(A_edges(j)-(x-sin(x))),1);
        end
    otherwise
        error([sectioning,' is not a known sectioning type.'])
end

A = (A_edges(2:end)+A_edges(1:(end-1)))/2;
B = linspace(0,pi,pts/2);
C = linspace(0,2*pi,pts);
for j = 1:sections
    a(:,:,j,1) = A_edges(j);
    a(:,:,j,2) = A(j);
    a(:,:,j,3) = A_edges(j+1);
end
for j = 1:pts/2
    b(:,j,:,:) = B(j);
end
for j = 1:pts
    c(j,:,:,:) = C(j);
end

e4 = cos(a);
e3 = sin(a).*cos(b);
e2 = sin(a).*sin(b).*cos(c);
e1 = sin(a).*sin(b).*sin(c);

% Selects a filter for cubic crystals when crystal_type is 'cubic' and for
%   orthorhomic crystals when crystal_type is 'orthorhombic'
filter = nan(pts,pts/2,sections,3);
if bounded
    filter(e3<0 | e2<0) = 1;
    switch crystal_type
        case 'cubic'
            filter(e4<abs(e1)/(sqrt(2)-1) | e4<abs(e2)/(sqrt(2)-1) | e4<abs(e3)/(sqrt(2)-1) | e4<abs(e1+e2+e3) | e4<abs(-e1+e2+e3) | e4<abs(e1-e2+e3) | e4<abs(e1+e2-e3)) = 1;
        case 'orthorhombic'
            filter(e4<abs(e1) | e4<abs(e2) | e4<abs(e3)) = 1;
        otherwise
            error([crystal_type,' is not a known crystal type.'])
    end
end

% Performs sectioning and stereographic projection of the fundamental zone
%   when 'stereographic', sectioning and area preserving projection of the 
%   fundamental zone when 'area_preserve', and orthographic projection of 
%   the fundamental zone when 'fundamental'
switch plot_type
    case 'stereographic'
        t = cell(1,sections);
        for r = 1:sections
            filter_t = (T(:,1)>=0 & T(:,4)>=0 & A_edges(r)<=acos(T(:,1)) & acos(T(:,1))<A_edges(r+1));
            t{r} = diag(1./sqrt(sum(T(filter_t,2:4).^2,2)))*T(filter_t,2:4);
        end
        
        half = floor(pts/4);
        e2 = 2*sin(b(:,1:half,1)).*cos(c(:,1:half,1))./(cos(b(:,1:half,1))+1);
        e1 = 2*sin(b(:,1:half,1)).*sin(c(:,1:half,1))./(cos(b(:,1:half,1))+1);
        x = linspace(-2,2,100);
        for j = 1:sections
            figure(j), clf;
            set(gcf,'position',[20*j+45,410-20*j,380,410]);
            hold on, grid on, box on;
            scatter3(2*t{j}(:,1)./(t{j}(:,3)+1),2*t{j}(:,2)./(t{j}(:,3)+1),ones(size(t{j}(:,3))),30,[0 0 1],'filled');
            for k = 1:size(filter,4)
                surf(e2.*filter(:,1:half,j,k),e1.*filter(:,1:half,j,k),filter(:,1:half,j,k),'edgealpha',0,'facealpha',0.1,'facecolor',[0 0 0]);
            end
            plot(x,sqrt(4-x.^2),'k','linewidth',1);
            plot(x,-sqrt(4-x.^2),'k','linewidth',1);

            axis equal, axis([-2,2,-2,2,0,2]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            title([num2str(A_edges(j)*360/pi),' to ',num2str(A_edges(j+1)*360/pi)]);
            view(0,90);
        end
    case 'area_preserve'
        a_q = real(acos(T(:,1)));
        b_q = real(acos(T(:,4)./sin(a_q)));
        b_q(abs(sin(a_q))<1e-8) = 0;
        c_q = real(acos(T(:,2)./(sin(a_q).*sin(b_q))));
        c_q(T(:,3)./(sin(a_q).*sin(b_q))<0) = 2*pi-c_q(T(:,3)./(sin(a_q).*sin(b_q))<0);
        c_q(abs(sin(a_q).*sin(b_q))<1e-8) = 0;
        
        t = cell(1,sections);
        for r = 1:sections
            f_t = (T(:,1)>=0 & T(:,4)>=0 & A_edges(r)<=acos(T(:,1)) & acos(T(:,1))<A_edges(r+1));
            t{r} = diag((3/4*(2*a_q(f_t)-sin(2*a_q(f_t)))).^(1/3).*sqrt(2*(1-cos(b_q(f_t)))))*[cos(c_q(f_t)) sin(c_q(f_t))];
        end
        
        half = floor(pts/4);
        e2 = sqrt(2*(1-cos(b(:,1:half,1)))).*cos(c(:,1:half,1));
        e1 = sqrt(2*(1-cos(b(:,1:half,1)))).*sin(c(:,1:half,1));
        r_max = (3/4*(2*A_max-sin(2*A_max)))^(1/3)*sqrt(2);
        for j = 1:sections
            figure(j), clf;
            set(gcf,'position',[20*j+45,410-20*j,380,410]);
            hold on, grid on, box on;
            
            shade = [A_edges(j) A(j) A_edges(j+1)];
            r = (3/4*(2*shade-sin(2*shade))).^(1/3);
            x = linspace(-sqrt(2)*r(3),sqrt(2)*r(3),100);
            
            scatter3(t{j}(:,1),t{j}(:,2),ones(size(t{j}(:,1))),30,[0 0 1],'filled');
            for k = 1:size(filter,4)
                surf(r(k)*e2.*filter(:,1:half,j,k),r(k)*e1.*filter(:,1:half,j,k),filter(:,1:half,j,k),'edgealpha',0,'facealpha',0.1,'facecolor',[0 0 0]);
            end
            plot3(x,real(sqrt(2*r(3)^2-x.^2)),2*ones(size(x)),'k','linewidth',1);
            plot3(x,-real(sqrt(2*r(3)^2-x.^2)),2*ones(size(x)),'k','linewidth',1);
            
            axis equal, axis([-r_max,r_max,-r_max,r_max,0,2]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            title([num2str(A_edges(j)*360/pi),' to ',num2str(A_edges(j+1)*360/pi)]);
            view(0,90);
        end
    case 'fundamental'
        swap = isnan(filter);
        filter = nan(size(filter));
        filter(swap) = 1;

        figure(1), clf;
        set(gcf,'position',[65,390,380,410]);
        hold on, grid on, box on;
        t = T(T(:,1)>=0,2:4);
        scatter3(t(:,1),t(:,2),t(:,3),30,[0 0 1],'filled');
        spr = (1:2:pts);
        surf(sin(b(spr,spr,1,1)).*cos(c(spr,spr,1,1)),sin(b(spr,spr,1,1)).*sin(c(spr,spr,1,1)),cos(b(spr,spr,1,1)),'edgealpha',.1,'facealpha',0)

        figure(2), clf;
        set(gcf,'position',[65,40,380,410]);
        hold on, grid on, box on;
        for k = find(A < pi/2 & squeeze(any(any(isnan(filter(:,:,:,2)),1),2))')
            surf(e2(:,:,k,2).*filter(:,:,k,2),e1(:,:,k,2).*filter(:,:,k,2),e3(:,:,k,2).*filter(:,:,k,2),'edgealpha',0.4,'facecolor',[.7 .7 .7]);
        end

        for j = 1:2
            figure(j);
            axis equal, axis([-1,1,-1,1,-1,1]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            view(140,30);
        end
    otherwise
        error([plot_type,' is not a known plotting type']);
end

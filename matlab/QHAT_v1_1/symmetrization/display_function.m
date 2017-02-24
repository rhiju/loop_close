% Various projections of the rth symmetrized harmonic of principal index n,
%   as calculated by make_symm.m or make_symm_loop.m. Symmetrized harmonic
%   is constructed from harmonics with indices given by lindex, mindex,
%   and iindex with coefficients given by weights. When plot_type is
%   'stereographic', the symmetrized harmonic is sectioned into concentric
%   spherical shells and projected stereographically into the plane, with
%   theta restricted to be less than pi/2. When plot_type is 
%   'area_preserve', the symmetrized harmonic is sectioned into concentric
%   spherical sheels and an area-preserving projection is used. When 
%   plot_type is 'fundamental', Figure 1 gives the portion of the 
%   symmetrized harmonic inside the fundamental zone and Figure 2 gives a 
%   rough representation of the boundary of the fundamental zone. When 
%   plot_type is 'isosurface', Figure 1 gives the same as for 'fundamental'
%   and Figure 2 gives the nodal surface inside the fundamental zone. 
%   Coloring set by the variable clr.
% 
% N - order of the representation, or equivalently index N
% r - number of the symmetrized harmonic for the given value of N
% bounded - constrains the function to the fundamental zone when true. When
%   false, the entire function is shown.
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
%   orthographically.  When 'isosurface', the nodal surface in the
%   fundamental zone is projected orthographically.
% pts - resolution of the numerical representation
% sections - number of surfaces of constant beta to plot
% clr - colormap that sets positive to yellow, negative to blue
%
% weight - nonzero coefficients of the real harmonics to be included in the
%   expansion of the input function
% lindex, mindex, iindex - l, m and i of the real harmonics with nonzero
%   coefficients to be included in the expansion of the input function. 
%   i = 0 corresponds to c, i = 1 corresponds to s.
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% P, p - associated Legendre polynomials
% G, g - normalized Gegenbauer polynomials
% A_max - maximum value of beta in the fundamental zone. Should be
%   restricted to less than pi/2.
% A_edges - controls the sectioning type as described in the code. Stores
%   the values of beta bounding concentric spherical shells.
% e4, e3, e2, e1 - w, z, x, y coordinates, respectively
% input - function to be visualized, defined by the expansion coefficients
%   in weights
% filter - array of ones and nans used to restrict plotting of surfaces
%   to regions within the fundamental zone. Currently implemented for
%   orthorhomic and cubic crystal symmetries, as described in the code.
% scale - maximum absolute values of the input function, used for
%   appropriate scaling of the color map to maximize contrast
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

clear; %clc;

N = input('Principal index: ');
r = input('Harmonic number: ');

bounded = false;
crystal_type = 'cubic';
sectioning = 'angle';
plot_type = 'area_preserve';

pts = 200;
sections = 12;
clr = [linspace(0,1,50)',ones(50,1),ones(50,1);1,1,1;ones(50,1),ones(50,1),linspace(1,0,50)'];

load(['eigenvectors/',num2str(N),'.mat']);
indices = (abs(sym(:,r+3)) > 1e-8);

weight = sym(indices,r+3)';
lindex = sym(indices,1)';
mindex = sym(indices,2)';
iindex = sym(indices,3)';

a = zeros(pts,pts/2,sections);
b = zeros(pts,pts/2,sections);
c = zeros(pts,pts/2,sections);
p = zeros(pts,pts/2,sections);
g = zeros(pts,pts/2,sections);

if bounded
    switch crystal_type
        case 'cubic'
            A_max = pi/5.733;
        case 'orthorhomic'
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
for j = 1:sections;
    a(:,:,j) = A(j);
end
for j = 1:pts/2;
    b(:,j,:) = B(j);
end
for j = 1:pts;
    c(j,:,:) = C(j);
end

e4 = cos(a);
e3 = sin(a).*cos(b);
e2 = sin(a).*sin(b).*cos(c);
e1 = sin(a).*sin(b).*sin(c);

input = zeros(pts,pts/2,sections);

% Construct input
for j = 1:size(weight,2)
    L = lindex(j);
    M = mindex(j);
    
    P = legendre(L,cos(B));
    G = 2^(L+1/2)*sqrt(((N+1)*factorial(N-L))/(pi*factorial(N+L+1)))*factorial(L)*sin(A).^L.*mfun('G',N-L,L+1,cos(A));
    for k = 1:sections;
        g(:,:,k) = G(k);
    end
    for k = 1:pts/2;
        p(:,k,:) = P(M+1,k);
    end
    if M > 0.5
        if iindex(j) < 0.5
            input = input+weight(j)*(-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*cos(M*c);
        else
            input = input+weight(j)*(-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*sin(M*c);
        end
    else
        input = input+weight(j)*g*sqrt((2*L+1)/(4*pi)).*p;
    end
end

% Selects a filter for cubic crystals when cubic_crystal is true and
%   for orthorhomic crystal symmetry when cubic_crystal is false
filter = nan(pts,pts/2,sections);
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

scale = eps+max(max(max(abs(input))));
% Performs sectioning and stereographic projection of the fundamental zone
%   when 'stereographic', orthographic projection of the fundamental zone
%   when 'fundamental', and nodal surface of the fundamental zone when
%   'isosurface'
switch plot_type
    case 'stereographic'
        half = floor(pts/4);
        e2 = 2*sin(b(:,1:half,1)).*cos(c(:,1:half,1))./(cos(b(:,1:half,1))+1);
        e1 = 2*sin(b(:,1:half,1)).*sin(c(:,1:half,1))./(cos(b(:,1:half,1))+1);
        for j = 1:sections
            figure(j), clf;
            set(gcf,'position',[20*j+45,410-20*j,380,410]);
            hold on, grid on, box on;

            surf(e2,e1,zeros(pts,half),input(:,1:half,j)+eps,'edgealpha',0);
            surf(e2.*filter(:,1:half,j),e1.*filter(:,1:half,j),filter(:,1:half,j),'edgealpha',0,'facealpha',0.3,'facecolor',[0 0 0]);

            pos = round(50*max(max(max(input(:,:,j))))/scale);
            neg = round(-50*min(min(min(input(:,:,j))))/scale);
            colormap(clr(51-neg:51+pos,:));

            axis equal, axis([-2,2,-2,2,0,2]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            title([num2str(A_edges(j)*360/pi),' to ',num2str(A_edges(j+1)*360/pi)]);
            view(0,90);
        end
    case 'area_preserve'
        half = floor(pts/4);
        e2 = sqrt(2*(1-cos(b(:,1:half,1)))).*cos(c(:,1:half,1));
        e1 = sqrt(2*(1-cos(b(:,1:half,1)))).*sin(c(:,1:half,1));
        r_max = (3/4*(2*A_max-sin(2*A_max)))^(1/3)*sqrt(2);
        for j = 1:sections
            figure(j), clf;
            set(gcf,'position',[20*j+45,410-20*j,380,410]);
            hold on, grid on, box on;
            
            r = (3/4*(2*A(j)-sin(2*A(j))))^(1/3);
            surf(r*e2,r*e1,zeros(pts,half),input(:,1:half,j)+eps,'edgealpha',0);
            surf(r*e2.*filter(:,1:half,j),r*e1.*filter(:,1:half,j),filter(:,1:half,j),'edgealpha',0,'facealpha',0.3,'facecolor',[0 0 0]);
            
            pos = round(50*max(max(max(input(:,:,j))))/scale);
            neg = round(-50*min(min(min(input(:,:,j))))/scale);
            colormap(clr(51-neg:51+pos,:));
            
            axis equal, axis([-r_max,r_max,-r_max,r_max,0,1]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            title([num2str(A_edges(j)*360/pi),' to ',num2str(A_edges(j+1)*360/pi)]);
            view(0,90);
        end
    case 'fundamental'
        swap = isnan(filter);
        filter = nan(size(filter));
        filter(swap) = 1;

        figure(2), clf;
        set(gcf,'position',[65,40,380,410]);
        hold on, grid on;
        for k = find(A < pi/2 & squeeze(any(any(isnan(filter),1),2))')
            surf(e2(:,:,k).*filter(:,:,k),e1(:,:,k).*filter(:,:,k),e3(:,:,k).*filter(:,:,k),'edgealpha',0.4,'facecolor',[.7 .7 .7]);
        end
        
        filter(abs(input)<0.1*scale) = nan;
        
        figure(1), clf;
        set(gcf,'position',[65,390,380,410]);
        hold on, grid on, box on;
        for k = find(A < pi/2)
            surf(e2(:,:,k).*filter(:,:,k),e1(:,:,k).*filter(:,:,k),e3(:,:,k).*filter(:,:,k),(input(:,:,k)+eps),'edgealpha',0,'facealpha',0.2);
        end
        pos = round(50*max(max(max(input(:,:,A < pi/2))))/scale);
        neg = round(-50*min(min(min(input(:,:,A < pi/2))))/scale);
        colormap(clr(51-neg:51+pos,:));

        for j = 1:2
            figure(j);
            axis equal, axis([-1,1,-1,1,-1,1]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            view(140,30);
        end
    case 'isosurface'
        swap = isnan(filter);
        filter = nan(size(filter));
        filter(swap) = 1;
        
        figure(2), clf;
        set(gcf,'position',[65,40,380,410]);
        hold on, grid on;
        p = patch(isosurface(e2,e1,e3,input.*filter,0));
        set(p,'FaceColor',[.6 .6 .6],'EdgeAlpha',0);
        daspect([1 1 1]), view(140,30);
        camlight, lighting gouraud;
        
        filter(abs(input)<0.1*scale) = nan;
        
        figure(1), clf;
        set(gcf,'position',[65,390,380,410]);
        hold on, grid on, box on;
        for k = find(A < pi/2)
            surf(e2(:,:,k).*filter(:,:,k),e1(:,:,k).*filter(:,:,k),e3(:,:,k).*filter(:,:,k),(input(:,:,k)+eps),'edgealpha',0,'facealpha',0.2);
        end
        pos = round(50*max(max(max(input(:,:,A < pi/2))))/scale);
        neg = round(-50*min(min(min(input(:,:,A < pi/2))))/scale);
        colormap(clr(51-neg:51+pos,:));

        for j = 1:2
            figure(j);
            axis equal, axis([-1,1,-1,1,-1,1]);
            xlabel('X'),ylabel('Y'),zlabel('Z');
            view(140,30);
        end
    otherwise
        error([plot_type,' is not a known plotting type']);
end

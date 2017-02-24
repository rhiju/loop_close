% Plots the real hyperspherical harmonics for the indices N, L, and M.
%   Figures 1 and 2 give portions of ZC and ZS, respectively, for rotation
%   angles less than pi. Figures 3 and 4 give portions of ZC and ZS,
%   respectively, for rotation angles greater than pi. Functions shown in 
%   orthographic projection from the unit sphere in four dimensions. 
%   Coloring set by the variable clr.
% 
% N, L, M - indices of the hyperspherical harmonic
% pts - resolution of the numerical representation
% clr - colormap that sets positive to yellow, negative to blue
% 
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% P, p - associated Legendre polynomials
% G, g - normalized Gegenbauer polynomials
% e4, e3, e2, e1 - w, z, x, y coordinates, respectively
% ZC, ZS - real hyperspherical harmonics for the indicated indices
% ZC_filter, ZS_filter - filters to avoid plotting function regions with
%   small absolute values
% ZC_scale, ZS_scale - maximum absolute values of the functions, used for
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

clear, clc;

N = 1;
L = 1;
M = 1;

pts = 50;
clr = [linspace(0,1,50)',ones(50,1),ones(50,1);1,1,1;ones(50,1),ones(50,1),linspace(1,0,50)'];

a = zeros(pts,pts,pts);
b = zeros(pts,pts,pts);
c = zeros(pts,pts,pts);
p = zeros(pts,pts,pts);
g = zeros(pts,pts,pts);

A = linspace(0,pi,pts);
B = linspace(0,pi,pts);
P = legendre(L,cos(B));
C = linspace(0,2*pi,pts);
G = 2^(L+1/2)*sqrt(((N+1)*factorial(N-L))/(pi*factorial(N+L+1)))*factorial(L)*sin(A).^L.*mfun('G',N-L,L+1,cos(A));
for j = 1:pts;
    a(:,:,j) = A(j);
    g(:,:,j) = G(j);
    b(:,j,:) = B(j);
    p(:,j,:) = P(M+1,j);
    c(j,:,:) = C(j);
end

e4 = cos(a);
e3 = sin(a).*cos(b);
e2 = sin(a).*sin(b).*cos(c);
e1 = sin(a).*sin(b).*sin(c);

if M > 0.5
    ZC = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*cos(M*c);
    ZS = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*sin(M*c);
else
    ZC = g*sqrt((2*L+1)/(4*pi)).*p;
    ZS = zeros(size(p));
end

ZC_filter = ones(pts,pts,pts);
ZC_filter(abs(ZC) < (eps+max(max(max(ZC)))-min(min(min(ZC))))*0.075) = nan;
ZC_scale = eps+max(max(max(abs(ZC))));

ZS_filter = ones(pts,pts,pts);
ZS_filter(abs(ZS) < (eps+max(max(max(ZS)))-min(min(min(ZS))))*0.075) = nan;
ZS_scale = eps+max(max(max(abs(ZS))));

figure(1), clf;
set(gcf,'position',[20,530,380,410]);
hold on, grid on;
for k = find(A < pi/2)
    surf(e2(:,:,k).*ZC_filter(:,:,k),e1(:,:,k).*ZC_filter(:,:,k),e3(:,:,k).*ZC_filter(:,:,k),(ZC(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(ZC(:,:,A < pi/2))))/ZC_scale);
neg = round(-50*min(min(min(ZC(:,:,A < pi/2))))/ZC_scale);
colormap(clr(51-neg:51+pos,:));

figure(2), clf;
set(gcf,'position',[20,40,380,410]);
hold on, grid on;
for k = find(A < pi/2)
    surf(e2(:,:,k).*ZS_filter(:,:,k),e1(:,:,k).*ZS_filter(:,:,k),e3(:,:,k).*ZS_filter(:,:,k),(ZS(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(ZS(:,:,A < pi/2))))/ZS_scale);
neg = round(-50*min(min(min(ZS(:,:,A < pi/2))))/ZS_scale);
colormap(clr(51-neg:51+pos,:));

figure(3), clf;
set(gcf,'position',[410,530,380,410]);
hold on, grid on;
for k = find(A > pi/2)
    surf(e2(:,:,k).*ZC_filter(:,:,k),e1(:,:,k).*ZC_filter(:,:,k),e3(:,:,k).*ZC_filter(:,:,k),(ZC(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(ZC(:,:,A > pi/2))))/ZC_scale);
neg = round(-50*min(min(min(ZC(:,:,A > pi/2))))/ZC_scale);
colormap(clr(51-neg:51+pos,:));

figure(4), clf;
set(gcf,'position',[410,40,380,410]);
hold on, grid on;
for k = find(A > pi/2)
    surf(e2(:,:,k).*ZS_filter(:,:,k),e1(:,:,k).*ZS_filter(:,:,k),e3(:,:,k).*ZS_filter(:,:,k),(ZS(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(ZS(:,:,A < pi/2))))/ZS_scale);
neg = round(-50*min(min(min(ZS(:,:,A < pi/2))))/ZS_scale);
colormap(clr(51-neg:51+pos,:));

for j = 1:4
    figure(j);
    axis equal, axis([-1,1,-1,1,-1,1]);
    xlabel('X'),ylabel('Y'),zlabel('Z');
    view(140,30);
end

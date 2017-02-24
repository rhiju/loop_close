% Visualization of the effect of left and right isoclinic rotations on an
%   arbitrary function on the unit sphere in four dimensions, performed by
%   expanding the function using the real harmonics and transforming the
%   expansion coefficients. Notice that the input function is defined only
%   after definition of the quantities a, b and c. Rotations are defined in
%   lax, rax, lan and ran. Figures 1 and 2 give portions of the input and
%   output functions, respectively, for rotation angles less than pi.
%   Figures 3 and 4 give portions of the input and output functions, 
%   respectively, for rotation angle greater than pi. Functions shown in
%   orthographic projection from the unit sphere in four dimensions.
%   Coloring set by the variable clr. 
% 
% N_max - maximum value of N used in representation of input function
% lax, rax - rotation axes for the left and right isoclinic rotations.
%   Left and right isoclinic rotations correspond to left and right 
%   multiplication by a quaternion.
% lan, ran - rotation angles for the left and right isoclinic rotations.
%   Left and right isoclinic rotations correspond to left and right 
%   multiplication by a quaternion.
% pts - resolution of the numerical representation
% clr - colormap that sets positive to yellow, negative to blue
% 
% lq, rq - quaternion modulus for the left and right isoclinic rotations
% lQ, rQ - pure quaternion for the left and right isoclinic rotations
% lR, rR - 2l+1 dimensional representations of the left and right isoclinic 
%   rotations, where l = N/2. Calculated using the complex representation.
%   Rows and columns ordered by decreasing m. Final index gives increasing
%   values of N, up to N_max.
% mpr, mup - primed row index and unprimed column index, respectively
% R - the (N+1)^2 dimensional representations of SO(4) given by combining
%   lR and rR using CG. Initially given for the complex representation with
%   increasing L, and decreasing M for a given L. Transformed using U to
%   be in the real representation with increasing L, by decreasing M for a
%   given L, and by C before S for a given M. Final index gives increasing
%   values of N, up to N_max.
% CG - Clebsh-Gordan coefficient matrix required to couple lR and rR.
%   Calculated using the complex representation. Columns ordered by 
%   increasing L, and decreasing M for a given L. Rows ordered for R to be
%   in the standard form for a four dimensional representation of SO(4)
%   when N = 1.
% U - similarity transformation to bring R to the real representation.
%   Rows ordered by increasing L, and by decreasing M for a given L.
%   Columns ordered by increasing L, by decreasing M for a given L, and by
%   C before S for a given M.
%
% A, a - beta, or half of the rotation angle
% B, b - theta, or the polar angle in spherical coordinates
% C, c - phi, or the azimuthal angle in spherical coordinates
% P, p - associated Legendre polynomials
% G, g - normalized Gegenbauer polynomials
% dp - differential volume element, required for calculation of expansion
%   coefficients for input function
% e4, e3, e2, e1 - w, z, x, y coordinates, respectively
% input - function to be rotated, defined as an explicit function of a, b
%   and c
% output - rotated function, defined by applying the rotation R to the
%   calculated expansion coefficients
% cc, cs - expansion coefficients of the input function, calculated by
%   integration with the real harmonics ZC and ZS, respectively. First
%   index corresponds to increasing N, second index corresponds to
%   increasing L, third index corresponds to increasing M.
% r - vector of real expansion coefficients for a given N. Ordered by
%   increasing L, decreasing M for a given L, and C before S for a given M.
% in_filter, out_filter - filters to avoid plotting function regions with
%   small absolute values
% in_scale, out_scale - maximum absolute values of the functions, used for
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

clc, clear;

N_max = 4;

lax = [1 1 1];
lan = pi/3;

rax = [1 1 1];
ran = -pi/3;

pts = 50;
clr = [linspace(0,1,50)',ones(50,1),ones(50,1);1,1,1;ones(50,1),ones(50,1),linspace(1,0,50)'];

% Construct quaternions for left and right isoclinic rotations
lax = lax/sqrt(lax*lax');
lq = cos(lan/2);
lQ = sin(lan/2).*lax;

rax = rax/sqrt(rax*rax');
rq = cos(ran/2);
rQ = sin(ran/2).*rax;

% Construct lR and rR
lR = zeros(N_max+1,N_max+1,N_max+1);
rR = zeros(N_max+1,N_max+1,N_max+1);
for N = 0:N_max
    l = N/2;
    for mpr = -l:l
        for mup = -l:l
            for k = max([0,mup-mpr]):min([l-mpr,l+mup])
                lR(l-mpr+1,l-mup+1,N+1) = lR(l-mpr+1,l-mup+1,N+1)+((lq-1i*lQ(3))^(l+mup-k)*(lq+1i*lQ(3))^(l-mpr-k)*(-lQ(2)-1i*lQ(1))^(mpr-mup+k)*(lQ(2)-1i*lQ(1))^(k))/(factorial(l+mup-k)*factorial(l-mpr-k)*factorial(mpr-mup+k)*factorial(k));
                rR(l-mpr+1,l-mup+1,N+1) = rR(l-mpr+1,l-mup+1,N+1)+((rq+1i*rQ(3))^(l+mup-k)*(rq-1i*rQ(3))^(l-mpr-k)*(rQ(2)+1i*rQ(1))^(mpr-mup+k)*(-rQ(2)+1i*rQ(1))^(k))/(factorial(l+mup-k)*factorial(l-mpr-k)*factorial(mpr-mup+k)*factorial(k));
            end
            lR(l-mpr+1,l-mup+1,N+1) = lR(l-mpr+1,l-mup+1,N+1)*sqrt(factorial(l+mpr)*factorial(l-mpr)*factorial(l+mup)*factorial(l-mup));
            rR(l-mpr+1,l-mup+1,N+1) = rR(l-mpr+1,l-mup+1,N+1)*sqrt(factorial(l+mpr)*factorial(l-mpr)*factorial(l+mup)*factorial(l-mup));
        end
    end
end

% Construct CG and R
R = zeros((N_max+1)^2,(N_max+1)^2,N_max+1);
for N = 0:N_max
    l = N/2;
    CG = zeros((N+1)^2,(N+1)^2);
    for L = 0:N
        for M = -L:L
            for lm = -l:l
                for rm = -l:l
                    if lm+rm == M
                        tau = 0;
                        for r = max([0,-L+l-lm,-L+l+rm]):min([N-L,l-lm,l+rm])
                            tau = tau+(-1)^r/(factorial(l-lm-r)*factorial(l+rm-r)*factorial(L-l+lm+r)*factorial(L-l-rm+r)*factorial(N-L-r)*factorial(r));
                        end
                        CG((N+1)*(l-lm)+(l-rm+1),(L+1)^2-L-M) = tau*sqrt(factorial(N-L)*factorial(L)^2*(2*L+1)/factorial(N+L+1))*sqrt(factorial(L+M)*factorial(L-M)*factorial(l+lm)*factorial(l-lm)*factorial(l+rm)*factorial(l-rm));
                    end
                end
            end
        end
    end
    R(1:(N+1)^2,1:(N+1)^2,N+1) = CG'*kron(rR(1:N+1,1:N+1,N+1),lR(1:N+1,1:N+1,N+1))*CG;
end

% Construct U
for N = 0:N_max
    U = zeros((N+1)^2,(N+1)^2);
    for L = 0:N
        U((L+1)^2-L,(L+1)^2) = (1i)^L;
        for M = 1:L
            U((L+1)^2-L-M,(L+1)^2-2*M) = (-1)^M*(1i)^L/sqrt(2);
            U((L+1)^2-L+M,(L+1)^2-2*M) = (1i)^L/sqrt(2);
            U((L+1)^2-L-M,(L+1)^2-2*M+1) = (-1)^M*(1i)^(L-1)/sqrt(2);
            U((L+1)^2-L+M,(L+1)^2-2*M+1) = -(1i)^(L-1)/sqrt(2);
        end
    end
    R(1:(N+1)^2,1:(N+1)^2,N+1) = U'*R(1:(N+1)^2,1:(N+1)^2,N+1)*U;
    R(1:(N+1)^2,1:(N+1)^2,N+1) = real(R(1:(N+1)^2,1:(N+1)^2,N+1));
end

a = zeros(pts,pts,pts);
b = zeros(pts,pts,pts);
c = zeros(pts,pts,pts);
p = zeros(pts,pts,pts);
g = zeros(pts,pts,pts);

A = linspace(0,pi,pts);
B = linspace(0,pi,pts);
C = linspace(0,2*pi,pts);
for j = 1:pts;
    a(:,:,j) = A(j);
    b(:,j,:) = B(j);
    c(j,:,:) = C(j);
end
dp = sin((a(1:pts-1,1:pts-1,1:pts-1)+a(2:pts,2:pts,2:pts))/2).^2*diff(A(1:2)).*sin((b(1:pts-1,1:pts-1,1:pts-1)+b(2:pts,2:pts,2:pts))/2)*diff(B(1:2))*diff(C(1:2));

e4 = cos(a);
e3 = sin(a).*cos(b);
e2 = sin(a).*sin(b).*cos(c);
e1 = sin(a).*sin(b).*sin(c);

input = sin(b).^2-1/2;
output = zeros(pts,pts,pts);
cc = zeros(N+1,N+1,N+1);
cs = zeros(N+1,N+1,N+1);

% Calculate expansion of input
for N = 0:N_max
    for L = 0:N
        for M = 0:L
            P = legendre(L,cos(B));     
            G = 2^(L+1/2)*sqrt(((N+1)*factorial(N-L))/(pi*factorial(N+L+1)))*factorial(L)*sin(A).^L.*mfun('G',N-L,L+1,cos(A));
            for j = 1:pts
                g(:,:,j) = G(j);
                p(:,j,:) = P(M+1,j);
            end
            if M > 0.5
                ZC = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*cos(M*c);
                ZS = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*sin(M*c);
            else
                ZC = g*sqrt((2*L+1)/(4*pi)).*p;
                ZS = zeros(size(p));
            end
            zc = input.*ZC;
            zs = input.*ZS;
            cc(N+1,L+1,M+1) = sum(sum(sum(dp.*(zc(1:pts-1,1:pts-1,1:pts-1)+zc(2:pts,2:pts,2:pts))/2)));
            cs(N+1,L+1,M+1) = sum(sum(sum(dp.*(zs(1:pts-1,1:pts-1,1:pts-1)+zs(2:pts,2:pts,2:pts))/2)));
        end
    end
end

% Perform rotation
for N = 0:N_max
    r = [];
    for L = 0:N
        sr = cc(N+1,L+1,1);
        for M = 1:L;
            sr = [cc(N+1,L+1,M+1) cs(N+1,L+1,M+1) sr];
        end
        r = [r sr];
    end
    r = R(1:(N+1)^2,1:(N+1)^2,N+1)*r';
    for L = 0:N
        sr = r((L+1)^2-2*L:(L+1)^2);
        cc(N+1,L+1,1) = sr(2*L+1);
        for M = 1:L;
            cc(N+1,L+1,M+1) = sr(2*L-2*M+1);
            cs(N+1,L+1,M+1) = sr(2*L-2*M+2);
        end
    end
end

% Construct output
for N = 0:N_max
    for L = 0:N
        for M = 0:L
            P = legendre(L,cos(B));      
            G = 2^(L+1/2)*sqrt(((N+1)*factorial(N-L))/(pi*factorial(N+L+1)))*factorial(L)*sin(A).^L.*mfun('G',N-L,L+1,cos(A));
            for j = 1:pts
                g(:,:,j) = G(j);
                p(:,j,:) = P(M+1,j);
            end
            if M > 0.5
                ZC = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*cos(M*c);
                ZS = (-1)^M*g*sqrt((2*L+1)/(2*pi)*factorial(L-M)/factorial(L+M)).*p.*sin(M*c);
            else
                ZC = g*sqrt((2*L+1)/(4*pi)).*p;
                ZS = zeros(size(p));
            end
            output = output+cc(N+1,L+1,M+1)*ZC+cs(N+1,L+1,M+1)*ZS;
        end
    end
end

% Plotting
in_filter = ones(pts,pts,pts);
in_filter(abs(input) < (eps+max(max(max(input)))-min(min(min(input))))*0.075) = nan;
in_scale = eps+max(max(max(abs(input))));

out_filter = ones(pts,pts,pts);
out_filter(abs(output) < (eps+max(max(max(output)))-min(min(min(output))))*0.075) = nan;
out_scale = eps+max(max(max(abs(output))));

figure(1), clf;
set(gcf,'position',[20,530,380,410]);
hold on, grid on;
for k = find(A < pi/2)
    surf(e2(:,:,k).*in_filter(:,:,k),e1(:,:,k).*in_filter(:,:,k),e3(:,:,k).*in_filter(:,:,k),(input(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(input(:,:,A < pi/2))))/in_scale);
neg = round(-50*min(min(min(input(:,:,A < pi/2))))/in_scale);
colormap(clr(51-neg:51+pos,:));

figure(2), clf;
set(gcf,'position',[20,40,380,410]);
hold on, grid on;
for k = find(A < pi/2)
    surf(e2(:,:,k).*out_filter(:,:,k),e1(:,:,k).*out_filter(:,:,k),e3(:,:,k).*out_filter(:,:,k),(output(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(output(:,:,A < pi/2))))/out_scale);
neg = round(-50*min(min(min(output(:,:,A < pi/2))))/out_scale);
colormap(clr(51-neg:51+pos,:));

figure(3), clf;
set(gcf,'position',[410,530,380,410]);
hold on, grid on;
for k = find(A > pi/2)
    surf(e2(:,:,k).*in_filter(:,:,k),e1(:,:,k).*in_filter(:,:,k),e3(:,:,k).*in_filter(:,:,k),(input(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(input(:,:,A > pi/2))))/in_scale);
neg = round(-50*min(min(min(input(:,:,A > pi/2))))/in_scale);
colormap(clr(51-neg:51+pos,:));

figure(4), clf;
set(gcf,'position',[410,40,380,410]);
hold on, grid on;
for k = find(A > pi/2)
    surf(e2(:,:,k).*out_filter(:,:,k),e1(:,:,k).*out_filter(:,:,k),e3(:,:,k).*out_filter(:,:,k),(output(:,:,k)+eps),'edgealpha',0,'facealpha',0.1);
end
pos = round(50*max(max(max(output(:,:,A < pi/2))))/out_scale);
neg = round(-50*min(min(min(output(:,:,A < pi/2))))/out_scale);
colormap(clr(51-neg:51+pos,:));

for j = 1:4
    figure(j);
    axis equal, axis([-1,1,-1,1,-1,1]);
    xlabel('X'),ylabel('Y'),zlabel('Z');
    view(140,30);
end

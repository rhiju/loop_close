function real_euler(rhsh)
% real_euler(shsh) - visualization in Euler angle space of the function
%   stored in rhsh by the expansion coefficients over the real harmonics. 
%   rhsh is the return value of real_coeffs.m. Figure 1 gives the 
%   isosurface of 1/4 the maximum value in the reduced Euler angle space 
%   for cubic materials, and the remaining figures give individual sections
%   of this space perpendicular to the phi_2 axis.
%
% rhsh - a cell array storing the coefficients of the expansion. The index
%   of the cell corresponds to N/2+1 and the coefficients within the cell
%   are ordered by increasing L, decreasing M for a given L, and with C
%   before S for a given M.
% N_max - maximum principal order of the real harmonics used
% pts_q - resolution of the real harmonics used for interpolation
% pts_e - desired resolution of the representation in Euler space
% sections - number of surfaces of constant beta to plot
% clr - colormap that sets positive to yellow, negative to blue
% A_e, a_e - phi_1, or the first Euler angle
% B_e, b_e - Phi, or the second Euler angle
% C_e, c_e - phi_e, or the third Euler angle
% q - quaternions components corresponding to Euler angle mesh points
% a_q - beta corresponding to Euler angle mesh points
% b_q - theta corresponding to Euler angle mesh points
% c_q - phi corresponding to Euler angle mesh points
% A_r - beta of the real harmonics used for interpolation
% B_r - theta of the real harmonics used for interpolation
% C_r, c_r - phi of the real harmonics used for interpolation
% input - function stored in rhsh, represented in the space of normalized
%   quaternions
% output - function stored in rhsh, represented in the space of Euler
%   angles
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

N_max = 2*(size(rhsh,1)-1);

pts_q = 200;
pts_e = 200;
sections = 19;
clr = [ones(50,1),linspace(0,1,50)',linspace(0,1,50)';1,1,1;linspace(1,0,50)',linspace(1,0,50)',ones(50,1)];

A_e = linspace(0,pi/2,pts_e);
B_e = linspace(0,pi/2,pts_e);
C_e = linspace(0,pi/2,sections);
a_e = zeros(pts_e,pts_e,sections);
b_e = zeros(pts_e,pts_e,sections);
c_e = zeros(pts_e,pts_e,sections);
for j = 1:pts_e
    a_e(j,:,:) = A_e(j);
    b_e(:,j,:) = B_e(j);
end
for j = 1:sections
    c_e(:,:,j) = C_e(j);
end

q = zeros(pts_e,pts_e,sections,4);
q(:,:,:,1) = cos(b_e/2).*cos((a_e+c_e)/2);
q(:,:,:,2) = sin(b_e/2).*cos((a_e-c_e)/2);
q(:,:,:,3) = sin(b_e/2).*sin((a_e-c_e)/2);
q(:,:,:,4) = cos(b_e/2).*sin((a_e+c_e)/2);

a_q = real(acos(q(:,:,:,1)));
b_q = real(acos(q(:,:,:,4)./sin(a_q)));
b_q(abs(sin(a_q))<1e-8) = 0;
c_q = real(acos(q(:,:,:,2)./(sin(a_q).*sin(b_q))));
c_q(q(:,:,:,3)./(sin(a_q).*sin(b_q))<0) = 2*pi-c_q(q(:,:,:,3)./(sin(a_q).*sin(b_q))<0);
c_q(abs(sin(a_q).*sin(b_q))<1e-8) = 0;

A_r = linspace(0,pi,pts_q/2);
B_r = linspace(0,pi,pts_q/2);
C_r = linspace(0,2*pi,pts_q);

c_r = zeros(pts_q,pts_q/2,pts_q/2);
g_r = zeros(pts_q,pts_q/2,pts_q/2);
p_r = zeros(pts_q,pts_q/2,pts_q/2);
for j = 1:pts_q
    c_r(j,:,:) = C_r(j);
end

input = zeros(pts_q,pts_q/2,pts_q/2);
for N = 0:2:N_max
    for L = 0:N
        for M = 0:L
            G_r = zeros(size(A_r));
            for s = 0:N-L
                G_r = G_r+((prod((L+1):(L+s))*prod((N-L-s+1):(N-s)))/prod(1:s))*cos((2*s-N+L)*A_r);
            end
            G_r = 2^(L+1/2)*realsqrt((N+1)/(pi*prod((N-L+1):(N+L+1))))*sin(A_r).^L.*G_r;
            P_r = legendre(L,cos(B_r));
            for s = 1:pts_q/2
                g_r(:,:,s) = G_r(s);
                p_r(:,s,:) = P_r(M+1,s);
            end
            
            if M > 0.5
                input = input+rhsh{N/2+1}((L+1)^2-2*M)*(-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*g_r.*p_r.*cos(M*c_r);
                input = input+rhsh{N/2+1}((L+1)^2-2*M+1)*(-1)^(M)*realsqrt((2*L+1)/(2*pi*prod((L-M+1):(L+M))))*g_r.*p_r.*sin(M*c_r);
            else
                input = input+rhsh{N/2+1}((L+1)^2)*(-1)^L*realsqrt((2*L+1)/(4*pi))*g_r.*p_r;
            end
        end
    end
end

output = interp3(B_r,C_r,A_r,input,b_q,c_q,a_q,'linear');
scale = eps+max(max(max(abs(output))));

figure(1), clf;
set(gcf,'position',[20,40,380,410]);
hold on, box on;
p = patch(isosurface(b_e,a_e,c_e,output,0.25*scale));
set(p,'FaceColor',[.6 .6 .6],'EdgeAlpha',0);
axis equal, axis([0,pi/2,0,pi/2,0,pi/2]);
set(gca,'zdir','reverse');
xlabel('\Phi'),ylabel('\phi_1'),zlabel('\phi_2');
view(140,30);
camlight(140,30), lighting gouraud;

for j = 1:sections
    figure(j+1), clf;
    set(gcf,'position',[20*j+45,410-20*j,380,410]);
    hold on, box on;
    
    surf(a_e(:,:,j),b_e(:,:,j),-ones(pts_e,pts_e),output(:,:,j)+eps,'edgealpha',0);
    plot([0,0,pi/2,pi/2,0],[0,pi/2,pi/2,0,0],'color',[0 0 0],'linewidth',1);
    pos = round(50*max(max(max(output(:,:,j))))/scale);
    neg = round(-50*min(min(min(output(:,:,j))))/scale);
    colormap(clr(51-neg:51+pos,:));
    
    axis equal, axis([0,pi/2,0,pi/2]);
    set(gca,'xaxislocation','top','ydir','reverse');
    xlabel('\phi_1'),ylabel('\Phi');
    text(0.7,0.8,[num2str(C_e(j)*180/pi),'\circ'],'fontname','Times New Roman','fontsize',30);
    view(0,90);
end

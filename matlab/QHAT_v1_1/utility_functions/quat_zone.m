function quat_zone(fig_num)
% quat_zone(fig_num) - plots an orthographic projection of the fundamental
%   zones of the orientation and misorientation spaces of the normalized
%   quaternion space for cubic crystal symmetry
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

rt = sqrt(1:100);

figure(fig_num), clf;
hold on, grid on;

c = linspace(-(2-rt(2))/4,(2-rt(2))/4,100);
b = sqrt((5-2*rt(2))*(1-c.^2))/rt(17);
a = sqrt((5-2*rt(2))*(1-c.^2))/rt(17);
for m = 1:2
    for n = 1:2
        plot3((-1)^m*a,(-1)^n*b,c,'k','linewidth',1);
        plot3(c,(-1)^m*a,(-1)^n*b,'k','linewidth',1);
        plot3((-1)^n*b,c,(-1)^m*a,'k','linewidth',1);
    end
end

a = linspace((2-rt(2))/4,rt(2)/4,100);
b = (-(4+rt(2))*a+sqrt(7*(3-rt(2))+(-31+8*rt(2))*a.^2))/7;
c = (a+sqrt(3+rt(2)-(5+2*rt(2))*a.^2))/(2+3*rt(2));
for m = 1:2
    for n = 1:2
        for o = 1:2
            plot3((-1)^m*a,(-1)^n*b,(-1)^o*c,'k','linewidth',1);
            plot3((-1)^o*c,(-1)^m*a,(-1)^n*b,'k','linewidth',1);
            plot3((-1)^n*b,(-1)^o*c,(-1)^m*a,'k','linewidth',1);
        end
    end
end

plot3([0,1/sqrt(2*(2+rt(2)))],[0,0],[0,0],'r','linewidth',1)
plot3([0,sqrt((5-2*rt(2))/17)],[0,sqrt((5-2*rt(2))/17)],[0,0],'r','linewidth',1)
plot3([0,rt(3)/6],[0,rt(3)/6],[0,rt(3)/6],'r','linewidth',1)
c = linspace(0,(2-rt(2))/4,100);
b = sqrt((5-2*rt(2))*(1-c.^2))/rt(17);
a = sqrt((5-2*rt(2))*(1-c.^2))/rt(17);
plot3(a,b,c,'r','linewidth',1)
b = linspace(0,sqrt((5-2*rt(2))/17),100);
a = sqrt((2-rt(2))*(1-b.^2))/2;
c = zeros(100,1);
plot3(a,b,c,'r','linewidth',1)
c = linspace((2-rt(2))/4,rt(3)/6,100);
b = (-2*c+sqrt(6-8*c.^2))/6;
a = (-2*c+sqrt(6-8*c.^2))/6;
plot3(a,b,c,'r','linewidth',1)
a = linspace(rt(3)/6,sqrt((5-2*rt(2))/17),100);
b = (-2*a+sqrt(6-8*a.^2))/6;
c = (-2*a+sqrt(6-8*a.^2))/6;
plot3(a,b,c,'r','linewidth',1)
b = linspace(sqrt((5-2*rt(2))/34),rt(2)/4,100);
a = (b+sqrt(3+rt(2)-(5+2*rt(2))*b.^2))/(2+3*rt(2));
c = (-(4+rt(2))*b+sqrt(7*(3-rt(2))+(-31+8*rt(2))*b.^2))/7;
plot3(a,b,c,'r','linewidth',1)
a = linspace(1/sqrt(2*(2+rt(2))),sqrt((5-2*rt(2))/17),100);
b = sqrt(1/2-(2+rt(2))*a.^2);
c = sqrt(1/2-(2+rt(2))*a.^2);
plot3(a,b,c,'r','linewidth',1)

axis equal, axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
xlabel('x'), ylabel('y'), zlabel('z');
set(gca,'xtick',-0.5:0.5:0.5,'ytick',-0.5:0.5:0.5,'ztick',-0.5:0.5:0.5);
view(140,30);

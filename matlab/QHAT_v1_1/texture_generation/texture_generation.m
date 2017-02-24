% Construction and visualization of various standard crystallographic
%   textures. The desired texture is specified in texture, and the
%   visualization in plot_type. Able to handle various sample and crystal
%   symmetries as defined by the symmetry group generators. Number of
%   crystals is specified in mult, and deviation from the ideal texture in
%   noise.
%
% texture - string that specifies the texture. Allowed values are 'cube',
%   'fiber', 'copper', 'brass'.
% fiber_axis - when texture is 'fiber', this specifies the common
%   crystallogrpahic axis of the drawing direction
% plot_type - string that specifies the visualization. Allowed values are
%   'pole' for pole figures, 'inverse' for inverse pole figures,
%   'stereographic' for a sectioning of the normalized quaterion space with
%   stereogrpahic projection, 'area_preserve' for a sectioning of the
%   normalized quaternion space with area-preserving projection,
%   'fundamental' or 'quaternion' for a projection of the normalized 
%   quaternion space, 'axis-angle' for the full axis-angle space,
%   'rodrigues' for the full Rodrigues vector space and individual cross 
%   sections, and 'euler' for the reduced Euler angle space and individual
%   cross sections.
% bounded - when plot_type is 'fundamental' or 'stereographic', set to
%   true to show only the fundamental zone and to false to show the full
%   space
% crystal_type - when plot_type is 'fundamental' or 'stereographic', set
%   to 'cubic' for cubic crystal and orthorhombic sample symmetry or to
%   'orthorhombic' for orthorhomic crystal and orthorhomic sample symmetry
% sxa, sxb - axes of rotations that generate sample symmetry group
% sna, snb - angles of rotations that generate sample symmetry group
% cxa, cxb - axes of rotations that generate crystal symmetry group
% cna, cnb - angles of rotations that generate crystal symmetry group
% mult - number of crystals to be included in the texture
% noise - angle of random rotation applied to oriented crystals. Angle is
%   sampled from a Gaussian distribution with specified standard deviation
%   for each crystal.
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

texture = 'fiber';
fiber_axis = [0 0 1];

plot_type = 'area_preserve';
bounded = false;
crystal_type = 'cubic';

% Sample symmetry generators
sxa = [0 0 1];
sna = pi;

sxb = [1 0 0];
snb = pi;

% Crystal symmetry generators
cxa = [0 0 1];
cna = pi/2;

cxb = [1 0 0];
cnb = pi/2;

% Symmetry groups
OL = group_symm(real_rotation(sxa,sna),real_rotation(sxb,snb));
OR = group_symm(real_rotation(cxa,cna),real_rotation(cxb,cnb));

mult = 100;
noise = pi/18;

switch texture
    case 'random'
        C = random(OL,OR,mult);
    case 'cube'
        C = cube(OL,OR,noise,mult);
    case 'fiber'
        C = fiber(OR,fiber_axis,noise,mult);
    case 'copper'
        C = rolling(OL,OR,true,noise,mult);
    case 'brass'
        C = rolling(OL,OR,false,noise,mult);
    otherwise
        error([texture,' is not a known texture type.'])
end

T = convert(C);

switch plot_type
    case 'pole'
        pole_figures(C);
    case 'inverse'
        inverse_pole_figures(C);
    case 'stereographic'
        quaternion_figures(T,bounded,crystal_type,'angle',plot_type);
    case 'area_preserve'
        quaternion_figures(T,bounded,crystal_type,'angle',plot_type);
    case 'fundamental'
        quaternion_figures(T,bounded,crystal_type,'angle',plot_type);
    case 'axis-angle'
        scaled_axis(T,plot_type);
    case 'rodrigues'
        scaled_axis(T,plot_type);
    case 'quaternion'
        scaled_axis(T,plot_type);
    case 'euler'
        euler(T);
    otherwise
        error([plot_type,' is not a known plot type.'])
end

--- Contents ---

1: Introduction

2: Citation Details

3: Mathematical Conventions

4: Usage
   4.0: Learning Scripts
   4.1: Symmetrization
      4.1.0: Eigenvectors
      4.1.1: Eigenfunctions
   4.2: Texture Generation
   4.3: Orientation Function
   4.4: Utility Functions

5: Correspondence

6: License

--- 1: Introduction ---

Quaternionic Harmonic Analysis of Texture is a collection of MATLAB scripts 
and functions used to manipulate and visualize crystallographic texture
information represented as a probability distribution of normalized
quaternions on the unit 3-sphere.

Quaternionic Harmonic Analysis of Texture is open-source software, and is
distributed under the GPLv2.  It was produced at Lawrence Livermore
National Laboratory by Jeremy Mason, reachable at jkylemason@gmail.com.

--- 2: Citation Details ---

Please cite the following PhD dissertation if you use Quaternionic Harmonic
Analysis of Texture in your research and/or software.

Mason, Jeremy K. 2009. Analysis of crystallographic texture information by
the hyperspherical harmonic expansion. PhD Thesis, Massachusetts Institute
of Technology, Boston. 230 p.

--- 3: Mathematical Conventions ---

Given the substantial opportunity for conflicting rotation conventions,
this software strives to use consistent conventions throughout.  The
following conventions should always be assumed unless explicitly stated
otherwise.

Rotations are interpreted in the active sense.  That is, the coordinate
system is aligned with the laboratory frame, and a rotation is a physical 
operation applied to an object with respect to the invariant frame.

Rotations act to the right.  That is, when a product of rotations is
applied to a column vector of coordinates, the rightmost remaining
operation is always applied first.  This convention applies equally to
rotations represented as matrices and as normalized quaternions.

A quaternion is always normalized.  The four components of a quaternion are
constructed from the axis $\overline a$ and angle $\omega$ of a rotation as
$[\cos(\omega/2), \overline a \cdot \sin(\omega/2)]$.  The first component
is occasionally referred to as the scalar part, and the final three
components are occasionally referred to as the vector part.

There is a two-to-one homomorphism between the group of proper three
dimensional rotations and the group of quaternions.  For the purposes of
visualizing a probability distribution of quaternions, our convention is to
often use only the quaternion with a positive scalar part.

The ideas implemented in this software are more thoroughly explained in
"Analysis of crystallographic texture information by the hyperspherical
harmonic expansion", available at http://hdl.handle.net/1721.1/53251.

--- 4: Usage ---

The files included in this release are organized into the following
collections.  More complete documentation of the use and implementation of
the individual scripts and functions is provided in the files themselves. 
Currently, the relative location for many of the scripts is rigid.

--- 4.0: Learning Scripts ---

The four scripts in ./learning_scripts are intended as an introduction to
the notation, the variable naming conventions and the visualization format. 
Our recommendation is to approach them in the following order:

plot_harmonic.m: plots one of the real hyperspherical harmonics

plot_combination.m: plots a linear combination of the real hyperspherical
harmonics for a given principal index

rotate_combination.m: rotates a linear combination of the real
hyperspherical harmonics for a given principal index

rotate_general.m: rotates an arbitrary function on the unit sphere in four
dimensions

--- 4.1: Symmetrization ---

The scripts in ./symmetrization are used for the calculation and
visualization of the symmetrized hyperspherical harmonics.  This is
generally considered to be a prerequisite for texture generation or the
calculation of orientation distribution functions.

display_coeffs.m: gives the expansion coefficients of one of the
symmetrized hyperspherical harmonics as human readable output.  Called
after execution of scripts in ./symmetrization/eigenvectors.

display function.m: gives various visualizations of projections of the
symmetrized hyperspherical harmonics.  Called after execution of scripts
in ./symmetrization/eigenvectors.

--- 4.1.0: Eigenvectors ---

The files in ./symmetrization/eigenvectors are used to calculate the
expansion coefficients of the symmetrized hyperspherical harmonics in the
basis provided by the real hyperspherical harmonics.

make_symm.m: calculates the expansion coefficients of the symmetrized
hyperspherical harmonics for a given value of the principal index.  Often
the only script that requires manipulation from the user.

make_symm_loop.m: calculates the expansion coefficients of the symmetrized
hyperspherical harmonics up to a maximum value of the principal index. 
Often the only script that requires manipulation from the user.

clebscgordan.m: constructs a sparse form of the Clebsch-Gordan coefficient
matrix

complexreal.m: constructs a sparse matrix that transforms from the complex
hyperspherical harmonics to the real hyperspherical harmonics

rotation.m: constructs the irreducible representation of SO(3) for the
indicated dimension

clean.m: searches for linear combinations of the symmetrized hyperspherical
harmonics that minimize the number of nonzero expansion coefficients

--- 4.1.1: Eigenfunctions ---

symm_construct.m: explicitly calculates the values of the symmetrized
hyperspherical harmonics defined by the expansion coefficients in 
./symmetrization/eigenvectors

--- 4.2: Texture Generation ---

The thirteen files in ./texture_generation construct and allow the
visualization of various standard crystallographic textures.  These may be
used as input for the functions in ./orientation_function.

texture_generation.m: a script that orchestrates the calling of various
other functions within ./texture_generation.  Often the only script that
requires manipulation from the user.

real_rotation.m: constructs the canonical representation of SO(3) for the
indicated rotation axis and angle

group_symm.m: constructs the full symmetry group generated by the input
matrices

random.m: returns a discrete random texture with the specified symmetries

cube.m: returns a discrete cube texture with the specified symmetries and
tolerance.  Generally used with cubic material symmetry.

fiber.m: return a discrete fiber texture with the specified symmetries,
fiber axis and tolerance

rolling.m: return a discrete copper or brass texture with the specified
symmetries and tolerance.  Generally used with cubic material and
orthorhombic sample symmetry.

convert.m: converts a cell array of rotation matrices to a cell array of
quaternions

pole_figures.m: constructs pole figures for the input texture

inverse_pole_figures.m: constructs inverse pole figures for the input
texture

quaternion_figures.m: constructs projections of the quaternion
representation of the input texture

scaled_axis.m: constructs visualizations of the input texture represented
as axis-angle pairs, Rodrigues vectors, or as quaternions

euler.m: constructs a visualization of the Euler angle representation of
the input texture

--- 4.3: Orientation Function ---

Given a discrete crystallographic texture, the functions in 
./orientation_function perform various manipulations of the corresponding
orientation distribution function (ODF).  The functions used for the
symmetrized hyperspherical harmonic expansion generally require the prior
calculation of the corresponding symmetrized hyperspherical harmonics in 
./symmetrization/eigenfunctions.

real_coeffs.m: calculates the expansion coefficients of the discrete
texture over the real hyperspherical harmonics

real_ghost.m: modifies the expansion coefficients from real_coeff.m to
suppress regions of negative probability density in the ODF

real_plot.m: projection and visualization from the space of quaternions of
the ODF expressed by a collection of expansion coefficients over the real
hyperspherical harmonics

real_euler.m: visualization in the space of Euler angles of the ODF
expressed by a collection of expansion coefficients over the real
hypserpherical harmonics

symm_coeffs.m: calculates the expansion coefficients of the discrete
texture over the symmetrized hyperspherical harmonics

symm_ghost.m: modifies the expansion coefficients from symm_coeff.m to
suppress regions of negative probability density in the ODF.

symm_plot.m: projection and visualization from the space of quaternions of
the ODF expressed by a collection of expansion coefficients over the
symmetrized hyperspherical harmonics.

symm_euler.m: visualization in the space of Euler angles of the ODF
expressed by a collection of expansion coefficients over the symmetrized
hypserpherical harmonics.

--- 4.4: Utility Functions ---

The six functions in ./utility_functions are intended to be used for
debugging purposes, or to follow the development in the original PhD
dissertation.

reduced_factorial.m: alternate calculation of products of factorials that
reduces the occurrence of overflow errors

clebschgordan.m: calculates the value of a Clebsch-Gordan coefficient,
following current quantum mechanical conventions

sixj.m: calculates the value of a Wigner 6-j symbol

axis_to_matr.m: calculates the three-dimensional rotation matrix of a
rotation expressed by a rotation axis and angle

matr_to_quat.m: robustly calculates the quaternion and the rotation axis
and angle corresponding to a given three-dimensional rotation matrix

quat_zone.m: plots an orthographic projection of the fundamental zones in
the space of normalized quaternions for cubic crystal symmetry

--- 5: Correspondence ---

Please direct all inquiries to Jeremy K. Mason, reachable at
jkylemason@gmail.com.

--- 6: License ---

Quaternionic Harmonic Analysis of Texture is open-source software, and is
distributed under the GPLv2.  Please refer to the "LICENSE.txt" file for
more complete license information.

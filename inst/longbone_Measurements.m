## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {csg-toolkit} {@var{measurements} =} longbone_Measurements ()
##
## Measurement names of long bone geometry data.
##
## @code{@var{measurements} = longbone_Measurements ()} returns a cell array of
## character vectors with the names of the longbone measurements in the same
## order as they appear in the @var{DATA} output argument returned by the
## @code{longbone_Geometry} function.
##
## @seealso{longbone_Geometry}
## @end deftypefn

function M = longbone_Measurements ()

  ## Check for valid number of input variables
  if (nargin > 0)
    error ("longbone_Measurements: does accept input arguments.");
  endif

  ## Check if input argument is provided
  M = {'Max Distance', ...
       'Area 20%', 'Perimeter 20%', 'ArPerIndex 20%', 'Ix 20%', 'Iy 20%', ...
       'Ixy 20%', 'Ix/Iy 20%', 'Imin 20%', 'Imax 20%', 'Imax/Imin 20%', ...
       'theta 20%', 'Dihedral angle 20_35', ......
       'Area 35%', 'Perimeter 35%', 'ArPerIndex 35%', 'Ix 35%', 'Iy 35%', ...
       'Ixy 35%', 'Ix/Iy 35%', 'Imin 35%', 'Imax 35%', 'Imax/Imin 35%', ...
       'theta 35%', 'Dihedral angle 35_50', ...
       'Area 50%', 'Perimeter 50%', 'ArPerIndex 50%', 'Ix 50%', 'Iy 50%', ...
       'Ixy 50%', 'Ix/Iy 50%', 'Imin 50%', 'Imax 50%', 'Imax/Imin 50%', ...
       'theta 50%', 'Dihedral angle 50_65', ...
       'Area 65%', 'Perimeter 65%', 'ArPerIndex 65%', 'Ix 65%', 'Iy 65%', ...
       'Ixy 65%', 'Ix/Iy 65%', 'Imin 65%', 'Imax 65%', 'Imax/Imin 65%', ...
       'theta 65%', 'Dihedral angle 65_80', ...
       'Area 80%', 'Perimeter 80%', 'ArPerIndex 80%', 'Ix 80%', 'Iy 80%', ...
       'Ixy 80%', 'Ix/Iy 80%', 'Imin 80%', 'Imax 80%', 'Imax/Imin 80%', ...
       'theta 80%', 'Diaphyseal Bending'};

endfunction

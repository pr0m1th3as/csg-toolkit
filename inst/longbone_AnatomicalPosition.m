## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {@var{v_rot} =} longbone_AnatomicalPosition (@var{v}, @var{f})
## @deftypefnx {csg-toolkit} {[@var{v_rot}, @var{vn_rot}] =} longbone_AnatomicalPosition (@var{v}, @var{f}, @var{vn})
##
## Rotate into anatomical posistion.
##
## This function automatically identifies the longbone defined by its vertices,
## @var{v}, and faces, @var{f}, rotates it into anatomical position and
## it translates its barycenter into the origin of the Cartesian coordinate
## axes, i.e. @qcode{[0, 0, 0]}.  The transformed (rotated/translated) vertices
## are returned in @var{v_rot}.
##
## The function can also take a third input argument, @var{vn}, defining the
## normals of the 3D mesh, which is rotated accordingly and returned as the
## second output argument @var{vn_rot}.
##
## @seealso{longbone_Geometry, meshRotation, readObj}
## @end deftypefn

function [varargout] = longbone_AnatomicalPosition (varargin)

  ## Check input output
  if (nargin < 2 || nargin > 3)
    print_usage;
  endif

  ## Get input arguments
  v = varargin{1};
  f = varargin{2};
  normals = false;
  if (nargin == 3)
    vn = varargin{3};
    normals = true;
  endif

  ## Check input arguments
  if (size (v, 2) != 3)
    error ("longbone_AnatomicalPosition: vertices must be an Nx3 matrix.");
  endif
  if (size (f, 2) != 3)
    error ("longbone_AnatomicalPosition: faces must be an Nx3 matrix.");
  endif
  if (normals)
    if (size (vn, 2) != 3)
      error ("longbone_AnatomicalPosition: normals must be an Nx3 matrix.");
    endif
  endif

  ## Translate barycenter to origin
  ## (it takes some iterations to minimize round off errors)
  for i = 1:10
    v = v - meshBarycenter (v, f);
  endfor

  ## Compute anatomical plane normals
  [TransPlane_normal, ~] = longbone_AnatomicalNormals (v, f);

  ## Rotate longtitudinal axis to [0,0,1]
  v_rot = meshRotation (v, TransPlane_normal, [0,0,1]);
  if (normals)
    vn_rot = meshRotation (vn, TransPlane_normal, [0,0,1]);
  endif

  ## Compute anatomical plane normals
  [~, CorPlane_normal] = longbone_AnatomicalNormals (v_rot, f);

  ## Rotate anteroposterior axis to [1,0,0]
  v_rot = meshRotation (v_rot, CorPlane_normal, [1,0,0]);
  if (normals)
    vn_rot = meshRotation (vn_rot, CorPlane_normal, [1,0,0]);
  endif

  ## Return arguments
  varargout{1} = v_rot;
  if (nargout > 1)
    if (normals)
      varargout{2} = vn_rot;
    else
      varargout{2} = [];
  endif

endfunction

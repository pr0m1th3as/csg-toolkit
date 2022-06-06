## Copyright (C) 2018-2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn{Function} [@var{GEOM}, @var{SMoA}, @var{polyline}] = polygon3D (@var{section}, @var{normal})
## @deftypefnx{Function} [@var{GEOM}, @var{SMoA}, @var{polyline}] = polygon3D (@var{section}, @var{normal}, @var{direction})
## @deftypefnx{Function} @var{GEOM} = polygon3D (@var{...})
## @deftypefnx{Function} [@var{GEOM}, @var{SMoA}] = polygon3D (@var{...})
##
## This function takes the 3D coordinates of a planar cross section along with
## with its normal and calculates some of its geometric properties such as cross
## sectional area, perimeter and principal second moments of area.  The points
## of the planar section may be unsorted, since the function will sort them in
## counter-clockwise order automatically, but need to constitute a simple
## polygon (a flat shape consisting of straight, non-intersecting line segments
## or "sides" that are joined pair-wise to form a closed path).
##
## If called with only 2 input arguments as in the following example
##
## @code{[GEOM, SMoA, polyline] = simple_polygon3D(section, normal)} @end code
## 
## Then,
##
## @multitable @columnfractions .25 .75
## @item @code{section} @tab is an Nx3 matrix containing the (x,y,z) point
## coordinates of a planar cross section.
## @item @code{normal} @tab is an 1x3 vector defining the unit normal of the
## section plane.
## @end multitable
##
## and,
##
## @multitable @columnfractions .25 .75
## @item @code{GEOM} @tab is a scalar structure containing the fields:
## @end multitable
##
## @multitable @columnfractions .05 .3 .65
## @item @tab @code{GEOM.Area} @tab calculated in mm^2
## @item @tab @code{GEOM.Perimeter} @tab calculated in mm
## @item @tab @code{GEOM.Centroid} @tab x,y,z centroid coordinates of the cross
## section.
## @end multitable
##
## @multitable @columnfractions .25 .75
## @item @code{SMoA} @tab is a scalar structure containing the fields:
## @end multitable
##
## @multitable @columnfractions .05 .3 .65
## @item @tab @code{SMoA.Ixy} @tab product of 2nd moments of area in mm^4.
## @item @tab @code{SMoA.Imin} @tab minimum 2nd moment of area in mm^4.
## @item @tab @code{SMoA.Imax} @tab maximum 2nd moment of area in mm^4.
## @end multitable
##
## @multitable @columnfractions .25 .75
## @item @code{polyline} @tab is a scalar structure containing the fields:
## @end multitable
##
## @multitable @columnfractions .05 .3 .65
## @item @tab @code{polyline.poly2D} @tab Nx2 matrix containing the (x,y)
## coordinates of the cross section on the 2D local axes of the sectioning
## plane ordered counter.
## @item @tab @code{polyline.poly3D} @tab Nx3 matrix containing the original 3D
## coordinates of the cross section ordered counter clockwise.
## @end multitable
##
## If called with 3 input arguments as in the following example
##
## @code{[GEOM, SMoA, polyline] = simple_polygon3D(section, normal, direction)}
## @end code
## 
## Then, @var{direction} is an 1x3 vector defining the unit normal of the
## coronal plane that provides alignment of the bone to anatomical position,
## while @var{normal} is assumed to point upwards.  The @code{SMoA} structure
## contains the additional fields, as follows:
##
## @multitable @columnfractions .05 .3 .65
## @item @tab @code{SMoA.Ix} @tab 2nd moment of area with respect to x axis,
## which is collinear with the intersection of the coronal plane and the
## sectioning plane.  Calculated in mm^4.
## @item @tab @code{SMoA.Iy} @tab 2nd moment of area with respect to y axis,
## which is collinear with the intersection of the sagital plane and the
## sectioning plane.  Calculated in mm^4.
## @item @tab @code{SMoA.theta} @tab angle of rotation of the principal axis
## Imax of 2nd moment of area with respect to x axis expressed in degrees
## in the range (-90, 90).
## @end multitable
## @seealso{meshArea, meshBarycenter, meshSection, longbone_CustomGeometry}
## @end deftypefn

function [varargout] = simple_polygon3D (varargin)
  ## check the number of input variables
  if nargin < 2 || length (varargin) > 3
    error 'wrong number of input arguments';
  endif
  ## for 2 input variables define matrix of planar section point coordinates and
  ## normal vector of slicing plane
  if nargin == 2
    section = varargin{1};
    normal = varargin{2};
  endif
  ## for 3 input variables define matrix of planar section point coordinates,
  ## normal vector of slicing plane and direction vector for local
  ## 2D coordinates of the slicing plane
  if nargin == 3
    section = varargin{1};
    normal = varargin{2};
    CorPlane = varargin{3};
  endif
  loc0 = section(1,:);
  ## find furthest point of polygon from loc0 to avoid rounding errors
  ## in order to calculate local x axis.
  [loc1_maxD, loc1_maxD_index] = max (sqrt (sum ((loc0 - section).^2, 2)));
  locx = section(loc1_maxD_index,:) - loc0;
  ## if anatomical position vector is present check if local x axis
  ## points towards the right hand side by find the normal vector of the plane
  ## perpendicular to the slicing plane which is also parallel to the coronal
  ## plane defined by the anatomical position vector
  if nargin==3
    locx_aligned = cross (CorPlane, normal);
    if sum(locx .* locx_aligned) < 0      ## loxc should point to the right
      ## inverse local x axis and swap local origin
      locx = locx * -1;
      loc0 = section(loc1_maxD_index,:);
    endif
  endif
  locy = cross (normal, locx);
  ## normalize local axes
  locx = locx ./ sqrt (sum (locx .^ 2));
  locy = locy ./ sqrt (sum (locy .^ 2));
  ## transform cross section to 2D coordinates
  polygon_2D = [sum((section - loc0) .* locx, 2), ...
                sum((section - loc0) .* locy, 2)];
  ## keep only unique points
  polygon_2D = unique (polygon_2D, "rows");
  ## sort points of polygon in counter clockwise order
  Midpoint = mean (polygon_2D);
  A1 = atan2 (polygon_2D(:,2) - Midpoint(2), polygon_2D(:,1) - Midpoint(1));
  polygon_2Ds = sortrows ([A1, polygon_2D], 1);
  polygon_2Ds = polygon_2Ds(:,[2:3]);
  ## iterate through the points and check for nearest neighbor between each
  ## pair of consecutive points to compensate for concave polygons that might
  ## produce strange polygon boundaries due to the initial arctan sorting
  len = length (polygon_2Ds);
  c_index = [1:len,1:len];
  for i = 1:len
    current_point = polygon_2Ds(i,:);
    next_point = polygon_2Ds(c_index(i+1),:);
    dist_next = distancePoints (current_point, next_point);
    compare_vector = polygon_2Ds(c_index([i+2:len*0.8+i]),:);
    distance2point = distancePoints (compare_vector, current_point);
    [minD, minD_index] = min (distance2point);
    minD_point = polygon_2Ds(c_index(i+1+minD_index),:);
    if (minD < dist_next)
      polygon_2Ds(c_index(i + 1),:) = minD_point;
      polygon_2Ds(c_index(i + 1 + minD_index),:) = next_point;
    endif
  endfor
  ## add the first point at the end to close the polygon
  poly = [polygon_2Ds; polygon_2Ds(1,:)];
  ## calculate the centroid, area and perimeter of the polygon
  Cx = 0; Cy = 0; Area = 0; Perimeter = 0;
  for i=1:length (polygon_2Ds)
    temp = ((poly(i,1) * poly(i + 1,2)) - (poly(i + 1,1) * poly(i,2)));
    Cx += (poly(i,1) + poly(i + 1,1)) * temp;
    Cy += (poly(i,2) + poly(i + 1,2)) * temp;
    Area += temp;
    Perimeter += sqrt((poly(i,1)-poly(i+1,1))^2 + (poly(i,2)-poly(i+1,2))^2);
  endfor
  Area = 0.5 * Area;
  Cx = (1 / (6 * Area)) * Cx;
  Cy = (1 / (6 * Area)) * Cy;
  ## transform centroid to original 3D coordinates
  Centroid = loc0 + (Cx * locx) + (Cy * locy);
  GEOM.Area = Area;
  GEOM.Perimeter = Perimeter;
  GEOM.Centroid = Centroid;
  ## check if second moments of area are required before proceeding
  if nargout == 1
    varargout{1} = GEOM;
  elseif nargout > 1
    ## move polygon's centroid to origin
    poly =  poly - [Cx, Cy];
    ## if anatomical position vector is present rotate the local coordinate
    ## axes to align locx axis with the coronal plane
    if (nargin == 3)
      ## find the angle between locx and locx_aligned to determine the
      ## required 2D rotation matrix for the polygon
      temp = cross (locx, locx_aligned);
      rot_A = atan2 (norm (temp), dot (locx, locx_aligned));
      ## check the direction of the angle
      sign = dot (normal, temp);
      if sign > 0
        rot_A = -rot_A;
      endif
      rot_M = [cos(rot_A), -sin(rot_A); sin(rot_A), cos(rot_A)];
      poly = (rot_M * poly')';
    endif
    ## calculate the parameters of second moments of area
    Ix = 0; Iy = 0; Ixy = 0;
    for i=1:length(polygon_2Ds)
      temp = (poly(i,1) * poly(i+1,2)) - (poly(i+1,1) * poly(i,2));
      Ix += ((poly(i,2)^2) + (poly(i,2) * poly(i+1,2)) + (poly(i+1,2)^2)) * temp;
      Iy += ((poly(i,1)^2) + (poly(i,1) * poly(i+1,1)) + (poly(i+1,1)^2)) * temp;
      Ixy += ((poly(i,1) * poly(i+1,2)) + (2 * poly(i,1) * poly(i,2)) + ...
             (2 * poly(i+1,1) * poly(i+1,2)) + (poly(i+1,1) * poly(i,2))) * temp;
    endfor
    Ix = (1 / 12) * Ix;
    Iy = (1 / 12) * Iy;
    Ixy = (1 / 24) * Ixy;
    ## if anatomical positioning vector is present return Ix and Iy values
    if (nargin == 3)
      SMoA.Ix = Ix;
      SMoA.Iy = Iy;
    endif
    SMoA.Ixy = Ixy;
    SMoA.Imin = ((Ix + Iy) * 0.5) - sqrt((((Ix - Iy) * 0.5) ^ 2) + (Ixy ^ 2));
    SMoA.Imax = ((Ix + Iy) * 0.5) + sqrt((((Ix - Iy) * 0.5) ^ 2) + (Ixy ^ 2));
    ## if anatomical positioning vector is present return the rotation angle
    ## theta of the principle axis Imax with respect to local x axis
    if (nargin == 3)
      if Ix >= Iy
        SMoA.theta = 0.5 * rad2deg (atan (-Ixy / ((Ix - Iy) / 2)));
      elseif (Ix < Iy) && (Ixy < 0)
        SMoA.theta = 0.5 * rad2deg (atan (-Ixy / ((Ix - Iy) / 2))) + 90;
      else
        SMoA.theta = 0.5 * rad2deg (atan (-Ixy / ((Ix - Iy) / 2))) - 90;
      endif
    endif
    varargout{1} = GEOM;
    varargout{2} = SMoA;
  endif
  if (nargout==3)
    ## store the cross sectional points in 2D coordinates with centroid
    ## centered to origin and properly oriented in direction vector was given
    polyline.poly2D = poly([1:end-1],:);
    ## transform cross sectional points back to original 3D coordinates
    polyline.poly3D = loc0 + (polygon_2Ds(:,1) * locx) + (polygon_2Ds(:,2) * locy);
    varargout{3} = polyline;
  endif
endfunction
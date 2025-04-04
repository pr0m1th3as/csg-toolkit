## Copyright (C) 2023-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} longbone_FragmentGeometry (@var{filename}, @var{points})
## @deftypefnx {csg-toolkit} {} longbone_FragmentGeometry (@var{folder}, @var{filename}, @var{points})
##
## This function analyzes the cross-sectional geometry of any diaphyseal long
## bone fragment at custom intervals defined by user defined points on the
## bone's surface.
##
## The fragment must be saved in Wavefront OBJ 3D model format, specified in
## @var{filename}, and the user defined points saved in a Meshlab PickedPoints
## file format, specified in @var{points}.  Optionally, the folder containing
## the 3D model can be defined in @var{folder}, if the function is invoked from
## a different working directory.  All thre input arguments must be strings.
##
## @code{longbone_FragmentGeometry} will process an arbitrary number of cross
## sections, according to the user defined points, as specified in @var{points},
## as long as they lie between 10% and 90% of the fragment's maximum length.
## Points beyond this range are ignored.
##
## @code{longbone_FragmentGeometry} does not return any output arguments, but it
## saves all geometric properties in CSV files as described in the following
## section.  3D models must be pure triangular meshes and their coordinate units
## are assumed to be in @math{mm}.  Each OBJ file must explicitly contain a
## single 3D model of any long bone fragment.
##
## Assuming a 3D model named @qcode{"bone_ID.obj"}, the following files are
## generated:
##
## @multitable @columnfractions .33 .02 .65
## @headitem Filename @tab @tab Description
## @item @qcode{Fgeometry-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{area} (measured in @math{mm^2}), @qcode{perimeter} (measured in
## @math{mm}), @qcode{centroid} (returned as [@var{x},@var{y},@var{z}]
## coordinates in @math{mm} units), and the sectioning @qcode{normal} (returned
## as an [@var{x},@var{y},@var{z}] vector). Each row of the CSV file corresponds
## to the user defined points at the same order as they appear in the relevant
## Meshlab PickedPoints file.  The aforementioned properties for each cross
## section are stored as a row vector with the respective order in columns
## [[2],[3],[4:6],[7:9]].  The first column contains the index of the
## corresponding cross-sectioning point in increasing order.
##
## @item @qcode{Finertia-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{Ixy}, @qcode{Imin}, and @qcode{Imax}, all measured in @math{mm^4}.
## Accordingly, they are saved in a row vector [2:4] with the first column
## containing the index of each cross section.
##
## @item @qcode{Fpolyline2D-bone_ID.csv} @tab @tab It contains the 2D
## coordinates (on an arbitrary @var{x},@var{y} local axis) for each cross
## section as an @math{Nx2} matrix, where @math{N} is the number of points for
## each cross-sectional polygon.  These polygons are ordered as column duplets,
## i.e. [[1:2],...,[end-1:end]], in the same order their corresponding points
## are stored in the Meshlab PickedPoints file.  Polygon coordinates start from
## the second row ([2:end],:), while the first row (1,[1,3,...,end-1])
## contains the index of the corresponding cross-sectioning points.
##
## @item @qcode{Fpolyline3D-bone_ID.csv} @tab @tab It contains the 3D model
## coordinates for each cross section as an @math{Nx3} matrix, where @math{N} is
## the number of points for each cross-sectional polygon.  These polygons are
## ordered as column triplets, i.e. [[1:3],...,[end-2:end]], in the same order
## their corresponding points are stored in the Meshlab PickedPoints file.
## Polygon coordinates start from starting from the second row ([2:end],:),
## while the first row (1,[1,4,...,end-2]) contains the index of the
## corresponding cross-sectioning points.
## @end multitable
##
## @seealso{longbone_Geometry, longbone_CustomGeometry, visualize_CrossSections}
## @end deftypefn

function longbone_FragmentGeometry (varargin)

  ## Check input
  if (nargin < 2 || nargin > 3)
    print_usage;
  endif

  ## Parse arguments
  if (nargin == 2)
    folder = pwd;
    filename = varargin{1};
    points = varargin{2};
  else
    folder = varargin{1};
    filename = varargin{2};
    points = varargin{3};
  endif

  ## Check input arguments
  if (! ischar (folder))
    error ("longbone_FragmentGeometry: FOLDER must be a string.");
  elseif (! (exist (folder) == 7))
    error ("longbone_FragmentGeometry: FOLDER is not a valid directory.");
  endif

  if (! ischar (filename))
    error ("longbone_FragmentGeometry: FILENAME must be a string.");
  elseif (! (strcmpi (filename([end - 3:end]), ".obj")))
    error ("longbone_FragmentGeometry: 3D model must be in OBJ file format.");
  elseif (! (exist (fullfile (folder, filename)) == 2))
    error ("longbone_FragmentGeometry: 3D model does not exist.");
  endif

  if (! ischar (points))
    error ("longbone_FragmentGeometry: POINTS must be a string.");
  elseif (exist (fullfile (folder, points)) == 2)
    ## Load Meshlab points for custom sectioning points from pp file
    CSP_data = read_MeshlabPoints (fullfile (folder, points));
    CSP_data(:,1) = [];
  else
    error (strcat (["longbone_FragmentGeometry: Meshlab Point file for"], ...
                   [" custom sectioning points does not exist."]));
  endif

  ## Load vertices and faces of triangular mesh from obj file
  [v,f] = readObj (fullfile (folder, filename));

  ## Keep only vertex coordinates (just in case color information is present)
  v = v(:,[1:3]);

  ## Find the maximum distance of the fragment
  [maxDistance, maxd_V1, maxd_V2] = longbone_maxDistance (v);

  ## Print bone's max distance
  page_screen_output (0);
  page_output_immediately (1);
  printf ("\nFragment in %s has a maximum distance of %f mm\n\n", ...
          filename, maxDistance);

  ## Calculate the normal vector of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt (sum ((maxd_V2 - maxd_V1) .^ 2));

  ## Calculate planes at 20 and 80% of the maximum length of the fragment
  point_1 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.20;
  point_2 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.35;
  point_3 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.50;
  point_4 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.65;
  point_5 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.80;

  ## Calculate the sectioning points for each plane
  plane_1 = meshSection (v, f, point_1, normal);
  plane_2 = meshSection (v, f, point_2, normal);
  plane_3 = meshSection (v, f, point_3, normal);
  plane_4 = meshSection (v, f, point_4, normal);
  plane_5 = meshSection (v, f, point_5, normal);

  ## Calculate centroid and area for each cross section
  CS_Geometry(1) = simple_polygon3D (plane_1, normal);
  CS_Geometry(2) = simple_polygon3D (plane_2, normal);
  CS_Geometry(3) = simple_polygon3D (plane_3, normal);
  CS_Geometry(4) = simple_polygon3D (plane_4, normal);
  CS_Geometry(5) = simple_polygon3D (plane_5, normal);

  ## Store the centroids of the initial cross-sectional areas
  Centroid_1 = CS_Geometry(1).Centroid;
  Centroid_2 = CS_Geometry(2).Centroid;
  Centroid_3 = CS_Geometry(3).Centroid;
  Centroid_4 = CS_Geometry(4).Centroid;
  Centroid_5 = CS_Geometry(5).Centroid;

  ## Calculate the projection of the user defined points to the longitudinal
  ## axis defined by maxd_V1 and maxd_V2 and keep the point that lie within
  ## the range of 10% to 90% of maximum length
  long_axis = maxd_V1 - maxd_V2;
  p = 0;
  for i = 1: size (CSP_data, 1)
    dist_CSP = CSP_data(i,:) - maxd_V2;
    point_proj = maxd_V2 + (dot (dist_CSP, long_axis) / ...
                            dot (long_axis, long_axis)) * long_axis;
    ## Check range
    dist = distancePoints (maxd_V2, point_proj);
    if (dist < maxDistance * 0.9 && dist > maxDistance * 0.1)
      p += 1;
      ratio = dist / maxDistance;
      ratio = round (ratio * 1000) / 1000;
      Centroid_N(p,:) = [ratio, point_proj];
    endif
  endfor

  ## Calculate new normals for each final cross section at 20, 35, 50, 65, 80%
  n1 = (Centroid_1 - Centroid_2) ./ sqrt (sum ((Centroid_1 - Centroid_2) .^2 ));
  n2 = (Centroid_1 - Centroid_3) ./ sqrt (sum ((Centroid_1 - Centroid_3) .^2 ));
  n3 = (Centroid_2 - Centroid_4) ./ sqrt (sum ((Centroid_2 - Centroid_4) .^2 ));
  n4 = (Centroid_3 - Centroid_5) ./ sqrt (sum ((Centroid_3 - Centroid_5) .^2 ));
  n5 = (Centroid_4 - Centroid_5) ./ sqrt (sum ((Centroid_4 - Centroid_5) .^2 ));

  ## Store centroids and normal in matrices
  CENTROIDS = [Centroid_1; Centroid_2; Centroid_3; Centroid_4; Centroid_5];
  NORMALS = [n1; n2; n3; n4; n5];

  ## Calculate new cross sections for every point in centroid_N by using the
  ## normal of the nearest default centroid (20%, 35%, 50%, 65%, 80%)
  for i = 1:size (Centroid_N, 1)
    ## Find nearest default centroid
    dist2centroid = distancePoints (Centroid_N(i,[2:4]), CENTROIDS);
    index = find (dist2centroid == min (dist2centroid));

    ## Calculate cross section with the normal of the nearest centroid
    section = meshSection (v, f, Centroid_N(i,[2:4]), NORMALS(index,:));
    [CS_Geometry, SMoA, polyline] = simple_polygon3D (section, NORMALS(index,:));

    ## Print results for centroids and cross sectional areas
    printf (strcat (["\nCross section at point %d has an area of %f mm^2,"], ...
                    [" perimeter of %f mm \n and centroid coordinates"], ...
                    [" are: x:%f y:%f z:%f\n"]), Centroid_N(i,1), ...
            CS_Geometry.Area, CS_Geometry.Perimeter, CS_Geometry.Centroid);

    ## Arrange to matrices for csv files
    geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, ...
                     CS_Geometry.Perimeter, CS_Geometry.Centroid, ...
                     NORMALS(index,:)];
    inertia(i,:) = [Centroid_N(i,1), SMoA.Ixy, SMoA.Imin, SMoA.Imax];
    polygon2D(1,i*2-1) = Centroid_N(i,1);
    polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
    polygon3D(1,i*3-2) = Centroid_N(i,1);
    polygon3D([2:length(polyline.poly3D)+1],[i*3-2:i*3]) = polyline.poly3D;
  endfor

  ## Save to files
  starting = "Fgeometry-";
  name = filename([1:length(filename) - 4]);
  endfile = ".csv";
  filename = strcat (starting, name, endfile);
  csvwrite (filename, geometry);
  starting = "Finertia-";
  filename = strcat (starting, name, endfile);
  csvwrite (filename, inertia);
  starting = "Fpolyline2D-";
  filename = strcat (starting, name, endfile);
  csvwrite (filename, polygon2D);
  starting = "Fpolyline3D-";
  filename = strcat (starting, name, endfile);
  csvwrite (filename, polygon3D);

endfunction

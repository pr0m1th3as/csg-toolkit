## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} longbone_CustomGeometry (@var{filename}, @var{points})
## @deftypefnx {csg-toolkit} {} longbone_CustomGeometry (@var{filename}, @var{bones}, @var{points})
## @deftypefnx {csg-toolkit} {} longbone_CustomGeometry (@var{folder}, @var{filename}, @var{points})
## @deftypefnx {csg-toolkit} {} longbone_CustomGeometry (@var{folder}, @var{filename}, @var{bones}, @var{points})
##
## This function analyzes the cross-sectional geometry of an intact humerus,
## ulna, femur, or tibia bone along its longitudinal axis at custom intervals.
##
## @code{longbone_CustomGeometry (@var{filename}, @var{points})} will analyze
## the 3D model in @var{filename}, which must be a char string, at custom
## intervals defined in@var{points}.  @var{points} can be either a char string
## or a numerical vector.  The string must specify a Meshlab PickedPoints file
## containing an arbitrary number of points on the bone's surface, each one of
## which specify a cross-sectional plane along the bones longitudinal axis.
## Alternatively, the numerical vector specifies an arbitrary number of
## cross-sectional planes as the ratios of the sectioning points' distance from
## the proximal end to the bone's maximum length in the range @math{[0.1, 0.9]}.
## Values beyond this range are ignored.
##
## @var{folder}, which must be a char string, defines the relative or absolute
## path to the directory containing the 3D model in @var{filename}.  When
## omitted, the current working directory is assumed.
##
## @var{bones} must be cell array of strings specifying one or more long bones
## that should be analyzed.  @code{longbone_CustomGeometry} will automatically
## determine what bone is represented in the 3D model, but it will analyze it
## only if it matches one of the bones named in @var{bones}.  Valid options are:
##
## @itemize
## @item @qcode{"Humerus"}
## @item @qcode{"Ulna"}
## @item @qcode{"Femur"}
## @item @qcode{"Tibia"}
## @item @qcode{"All"}
## @end itemize
##
## @code{longbone_CustomGeometry} does not return any output arguments, but it
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
## @item @qcode{Cgeometry-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{area} (measured in @math{mm^2}), @qcode{perimeter} (measured in
## @math{mm}), @qcode{centroid} (returned as [@var{x},@var{y},@var{z}]
## coordinates in @math{mm} units), as well as the sectioning and orientation
## @qcode{normals} (returned as an [@var{x},@var{y},@var{z}] vectors).  Each row
## of the CSV file corresponds to a different cross section.  The first and last
## rows contain the default cross sections at 20% and 80% respectively.  The
## intermediate rows follow the order of the sectioning points in the Meshlab
## PickedPoints file or in the numerical vector parsed in @var{points}.
## The aforementioned properties for each cross section are stored as a row
## vector with the respective order in columns [[2],[3],[4:6],[7:9],[10:12]].
## The first column contains the ratios of the sectioning points' distance from
## the proximal end to the bone's maximum length.
##
## @item @qcode{Cinertia-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{Ix}, @qcode{Iy} @qcode{Ixy}, @qcode{Imin}, and @qcode{Imax}, all
## measured in @math{mm^4}, as well as @qcode{theta}, measured in degrees.
## Accordingly, they are saved as row vectors [2:7] with the first column
## containing the ratios of the sectioning points' distance from the proximal
## end to the bone's maximum length.
##
## @item @qcode{Cpolyline2D-bone_ID.csv} @tab @tab It contains the 2D
## coordinates (on an arbitrary @var{x},@var{y} local axis) for each cross
## section as an @math{Nx2} matrix, where @math{N} is the number of points for
## each cross-sectional polygon.  These polygons are ordered as column duplets,
## i.e. [[1:2],...,[end-1:end]].  The first and the last column duplets contain
## the default cross sections at 20% and 80%, respectively.  The intermediate
## column duplets follow the order of the sectioning points in the Meshlab
## PickedPoints file or in the numerical vector parsed in @var{points}.
## Polygon coordinates start from the second row ([2:end],:), while the first
## row contains the relevant ratios at columns (1, [1,3,...,end-1]).
##
## @item @qcode{Cpolyline3D-bone_ID.csv} @tab @tab It contains the 3D model
## coordinates for each cross section as an @math{Nx3} matrix, where @math{N} is
## the number of points for each cross-sectional polygon.  These polygons are
## ordered as column triplets, i.e. [[1:3],...,[end-2:end]].  The first and the
## last column triplets contain the default cross sections at 20% and 80%,
## respectively.  The intermediate column triplets follow the order of the
## sectioning points in the Meshlab PickedPoints file or in the numerical vector
## parsed in @var{points}.  Polygon coordinates start from the second row
## ([2:end],:), while the first row contains the relevant ratios at columns
## (1, [1,3,...,end-1]).
## @end multitable
##
## For any analyzed 3D model, the initial alignment points are either
## read from the corresponding Meshlab PickedPoints file (e.g. 'bone_ID.pp'), if
## present in the same folder with the OBJ, or they are automatically registered
## with the @code{longbone_Registration} function.  The optimized alignment
## points are appended in the existing MPP file or saved in newly created,
## accordingly.
##
## @seealso{longbone_Geometry, longbone_FragmentGeometry,
## visualize_CrossSections}
## @end deftypefn

function longbone_CustomGeometry (varargin)

  ## Check input
  if (nargin < 2 || nargin > 4)
    print_usage;
  endif

  ## 1st arg is filename
  if (nargin == 2 && ischar (varargin{1}))
    folder = "";
    filename = varargin{1};
    bone = {"All"};
  elseif (nargin == 2)
    error ("longbone_CustomGeometry: FILENAME must be a string.");
  endif

  ## 1st and 2nd args are either filename and bone or folder and filename
  if (nargin == 3 && ischar (varargin{1}) && ischar (varargin{2}))
    folder = varargin{1};
    filename = varargin{2};
    bone = {"All"};
  elseif (nargin == 3 && ischar (varargin{1}) && iscellstr (varargin{2}))
    folder = "";
    filename = varargin{1};
    bone = varargin{2};
  elseif (nargin == 3)
    if (! ischar (varargin{1}) && ischar (varargin{2}))
      error ("longbone_CustomGeometry: FOLDER must be a string.");
    endif
    if (! ischar (varargin{1}) && iscellstr (varargin{2}))
      error ("longbone_CustomGeometry: FILENAME must be a string.");
    endif
    if (ischar (varargin{1}) && ! (ischar (varargin{2}) ||
                                   iscellstr (varargin{2})))
      error (strcat (["longbone_CustomGeometry: second argument must be"], ...
                     [" either a string for FILENAME or a cell array of"], ...
                     [" strings for selected BONES."]));
    endif
  endif

  ## 1st, 2nd and 3rd args are folder, filename, and bone
  if (nargin == 4 && ischar (varargin{1}) && ischar (varargin{2})
                  && iscell (varargin{3}))
    folder = varargin{1};
    filename = varargin{2};
    bone = varargin{3};
  elseif (nargin == 4)
    if (! ischar (varargin{1}))
      error ("longbone_CustomGeometry: FOLDER must be a string.");
    endif
    if (! ischar (varargin{2}))
      error ("longbone_CustomGeometry: FILENAME must be a string.");
    endif
    if (! iscellstr (varargin{3}))
      error ("longbone_CustomGeometry: BONES must be a cell array of strings.");
    endif
  endif

  ## Check last argument being a string or a numerical vector
  if (ischar (varargin{end}))
    filenamePP = varargin{end}(:)';
    if (exist (fullfile (folder, filenamePP)) == 2)
      ## Load Meshlab points for custom sectioning points from pp file
      CSP_data = read_MeshlabPoints (fullfile (folder, filenamePP));
      CSP_data(:,1) = [];
      CSP_points = true;
    else
      error (strcat (["longbone_CustomGeometry: Meshlab Point file for"], ...
                     [" custom sectioning points does not exist."]));
    endif
  elseif (isnumeric (varargin{end}))
    CSP_data = varargin{end}(:)';
    CSP_points = false;
  else
    error (strcat (["longbone_CustomGeometry: last input argument must"], ...
                   [" be a string or a numerical vector."]));
  endif

  ## Check filename has a valid .obj extension
  if (! (strcmpi (filename([end - 3:end]), ".obj")))
    error ("longbone_CustomGeometry: 3D model must be in OBJ file format.");
  endif

  ## Load vertices and faces of triangular mesh from obj file
  [v,f] = readObj (fullfile (folder, filename));

  ## Keep only vertex coordinates (just in case color information is present)
  v = v(:,[1:3]);

  ## Check if orientation points are available in corresponding pp file
  ## otherwise automatically register new ones on the bone surface
  ## for Meshlab point files
  filenamePP = filename([1:length(filename) - 4]);
  extension = ".pp";
  filenamePP = strcat (filenamePP, extension);
  if (exist (fullfile (folder, filenamePP)) == 2)
    ## Load Meshlab points for mediolateral axis from pp file
    MLA_points = read_MeshlabPoints (fullfile (folder, filenamePP));
    MLA_points(:,1) = [];
    find_bone = true;
  else
    [bonesel, MLA_points] = longbone_Registration (v, f);
    find_bone = false;
  endif

  ## Find bone if necessary
  if (find_bone)
    bonesel = longbone_Registration (v, f);
  endif

  ## Check bone selection
  if (any (strcmpi (bone, "All")))
    bone = bonesel;
  else
    if (numel (bone) == 1 && ! strcmpi (bone, bonesel))
      printf ("Model %s is not a %s\n", filename, bone);
      return;
    endif
    if (numel (bone) == 2 && ! (any (strcmpi (bone, bonesel))))
      printf ("Model %s is neither a %s nor a %s\n", ...
             filename, bone{1}, bone{2});
      return;
    endif
    if (numel (bone) == 3 && ! (any (strcmpi (bone, bonesel))))
      printf ("Model %s is not a %s, a %s or a %s\n", ...
              filename, bone{1}, bone{2}, bone{3});
      return;
    endif
    bone = bonesel;
  endif

  ## Check valid bone selection
  if (! (strcmpi (bone, "Humerus") || strcmpi (bone, "Ulna") ||
         strcmpi (bone, "Femur") || strcmpi (bone, "Tibia")))
    printf ("Bone should be a Humerus, an Ulna, a Femur or a Tibia\n");
    return;
  endif

  ## Find the maximum distance of the bone
  [maxDistance, maxd_V1, maxd_V2] = longbone_maxDistance (v);

  ## Print bone's max distance
  page_screen_output (0);
  page_output_immediately (1);
  printf ("\n%s in %s has a maximum distance of %f mm\n\n", ...
          bone, filename, maxDistance);

  ## Calculate the normal vector of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt (sum ((maxd_V2 - maxd_V1) .^ 2));

  ## Calculate planes at 20 and 80% of the maximum length of the bone
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

  ## Check the centroids' locations from proximal (20%) to distal (80%) and
  ## if required reverse their order to comply with this standard.
  proximal = distancePoints (Centroid_1, MLA_points(1,:));
  distal = distancePoints (Centroid_5, MLA_points(1,:));
  prox_Edge = maxd_V2;
  dist_Edge = maxd_V1;
  if ((strcmpi (bone, "Humerus") || strcmpi (bone, "Femur")) &&
      (proximal < distal))
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
    prox_Edge = maxd_V1;
    dist_Edge = maxd_V2;
  elseif ((strcmpi (bone, "Tibia") || strcmpi (bone, "Ulna")) &&
          (proximal > distal))
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
    prox_Edge = maxd_V1;
    dist_Edge = maxd_V2;
  endif

  ## Check for custom sectioning points from pp file
  if (CSP_points)
    ## Calculate the projection of the user defined points to the longitudinal
    ## axis defined by maxd_V1 and maxd_V2 and keep the point that lie within
    ## the range of 10% to 90% of maximum length
    long_axis = prox_Edge - dist_Edge;
    p = 1;
    ## First centroid at 20%
    Centroid_N(p,:) = [0.2, Centroid_1];
    for i = 1: size (CSP_data, 1)
      dist_CSP = CSP_data(i,:) - dist_Edge;
      point_proj = dist_Edge + (dot (dist_CSP, long_axis) / ...
                                dot (long_axis, long_axis)) * long_axis;
      ## Check range
      dist = distancePoints (dist_Edge, point_proj);
      if (dist < maxDistance * 0.9 && dist > maxDistance * 0.1)
        p += 1;
        ratio = dist / maxDistance;
        ratio = round (ratio * 1000) / 1000;
        Centroid_N(p,:) = [ratio, point_proj];
      endif
    endfor
    ## Last centroid at 80%
    p += 1;
    Centroid_N(p,:) = [0.8, Centroid_5];
  else
    ## Keep the points that lie within the range of 10% to 90% of maximum length
    p = 1;
    ## First centroid at 20%
    Centroid_N(p,:) = [0.2, Centroid_1];
    for i = 1: size (CSP_data, 2)
      if (CSP_data(i) <= 0.9 && CSP_data(i) >= 0.1)
        p += 1;
        point_N = maxd_V1 + (maxd_V2 - maxd_V1) .* CSP_data(i);
        plane_N = meshSection (v, f, point_N, normal);
        CS_Geometry = simple_polygon3D (plane_N, normal);
        Centroid_N(p,:) = [CSP_data(i), CS_Geometry.Centroid];
      endif
    endfor
    ## Last centroid at 80%
    p += 1;
    Centroid_N(p,:) = [0.8, Centroid_5];
  endif

  ## Calculate normals of anatomical planes
  if (! strcmpi (bone, "Ulna"))         # for humerus, femur, and tibia
    ## Calculate the mediolateral axis vector
    MLA_vector = MLA_points(1,:) - MLA_points(2,:);
    ## Calculate longitudinal axis vector
    PDA_vector = Centroid_1 - Centroid_5;
    ## Calculate the unit normal of the transverse plane facing upwards
    TransPlane_normal = PDA_vector ./ sqrt (sum (PDA_vector .^ 2));
    ## Calculate the vector of the coronal plane
    CorPlane_normal = cross (PDA_vector, MLA_vector);
    ## Normalize the vector of the coronal plane
    CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
  else                                  # for ulna
    ## No optimization of user defined point for ulna
    ## Calculate the projection of user defined point to the longitudinal axis
    C5_MLA = MLA_points(1,:) - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    point_proj = Centroid_1 + (dot (C5_MLA, C5_C1) / dot (C5_C1, C5_C1)) * C5_C1;
    ## Recalculate the coronal plane unit vector by normalizing the vector from
    ## midpoint_proj to MLA_midpoint
    CorPlane_normal = MLA_points(1,:) - point_proj;
    CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
    ## Calculate the unit normal of the transverse plane facing upwards
    TransPlane_normal = C5_C1 ./ sqrt (sum (C5_C1 .^ 2));
    ## Calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    MLA_opt_point_A = MLA_points(1,:);
    MLA_opt_point_B = MLA_points(2,:);
  endif

  ## Find if coronal plane normal points in the right direction (towards the
  ## front) by checking the dot product of the normal with the vector between
  ## the first point of MLA vector and the nearest centroid. If not, reverse
  ## the normal so that it points to the front

  ## HUMERUS
  if (strcmpi (bone, "Humerus"))
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if (sum (MLP1_C5 .* CorPlane_normal) > 0)   ## MLA_points should be in front
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## Find the points that lie on or below the transverse plane at C5
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find (TP_dotP <= 0),:);
    ## Find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 1 and 5 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + ...
                    (dot (C1_midpoint, C1_C5) / dot (C1_C5, C1_C5)) * C1_C5;
    ## Recalculate the coronal plane unit vector by normalizing the vector from
    ## midpoint_proj to MLA_midpoint
    CorPlane_normal = MLA_midpoint - midpoint_proj;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    ## Calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## Find the points that lie in front of the coronal plane at midpoint_proj
    CP_dotP = sum (CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_front = v_below(find (CP_dotP > 0),:);
    ## Find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_below_front - MLA_midpoint), 2);
    v_side1 = v_below_front(find (ML_dotP <= 0),:);
    v_side2 = v_below_front(find (ML_dotP >= 0),:);
    ## Optimize alignment points (5 iterations)
    for i = 1:5
      ## Calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj), 2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj), 2);
      ## Find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane
      [d1,i1] = max (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = max (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## Calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## Normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt( sum (CorPlane_normal .^ 2));
      ## Check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if (sum (MLP1_C5 .* CorPlane_normal) > 0) ## MLA_points should be in front
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor

  ## FEMUR
  elseif (strcmpi (bone, "Femur"))
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if (sum (MLP1_C5 .* CorPlane_normal) < 0)     ## MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## Find the points that lie on or below the transverse plane at C5
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find (TP_dotP <= 0),:);
    ## Find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 1 and 5 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + ...
                    (dot (C1_midpoint, C1_C5) / dot (C1_C5, C1_C5)) * C1_C5;
    ## Recalculate the coronal plane unit vector by normalizing the vector from
    ## MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    ## Calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## Find the points that lie behind the coronal plane of the midpoint_proj
    CP_dotP = sum (CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_behind = v_below(find (CP_dotP < 0),:);
    ## Find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_below_behind - MLA_midpoint), 2);
    v_side1 = v_below_behind(find (ML_dotP <= 0),:);
    v_side2 = v_below_behind(find (ML_dotP >= 0),:);
    ## Optimize alignment points (5 iterations)
    for i = 1:5
      ## Calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj), 2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj), 2);
      ## Find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane
      [d1,i1] = min (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = min (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## Calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## Normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
      ## Check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if (sum (MLP1_C5 .* CorPlane_normal) < 0)   ## MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor

  ## TIBIA
  elseif (strcmp (bone, "Tibia"))
    MLP1_C1 = Centroid_1 - MLA_points(1,:);
    if (sum (MLP1_C1 .* CorPlane_normal) < 0)     ## MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## Find the points that lie on or above the transverse plane at C1
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_1), 2);
    v_above = v(find (TP_dotP >= 0),:);
    ## Find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 5 and 1 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C5_midpoint = MLA_midpoint - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    midpoint_proj = Centroid_5 + ...
                    (dot (C5_midpoint, C5_C1) / dot (C5_C1, C5_C1)) * C5_C1;
    ## Recalculate the coronal plane unit vector by normalizing the vector from
    ## MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal .^ 2));
    ## Calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## Find the points that lie behind the coronal plane of the midpoint_pro
    CP_dotP = sum (CorPlane_normal .* (v_above - midpoint_proj), 2);
    v_above_behind = v_above(find (CP_dotP < 0),:);
    ## Find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_above_behind - MLA_midpoint), 2);
    v_side1 = v_above_behind(find (ML_dotP <= 0),:);
    v_side2 = v_above_behind(find (ML_dotP >= 0),:);
    ## Optimize alignment points (5 iterations)
    for i = 1:5
      ## Calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj),2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj),2);
      ## Find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane
      [d1,i1] = min (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = min (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## Calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## Normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal .^ 2));
      ## Check the optimized CorPlane_normal points towards the right direction
      MLP1_C1 = Centroid_1 - MLA_optimal_point1;
      if (sum (MLP1_C1 .* CorPlane_normal) < 0)   ## MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  endif

  ## Print information about alignment points
  if (! strcmp (bone, "Ulna"))         # for humerus, femur, and tibia
    if (! register)
      ## Print initial user defined MLA_points
      printf ("\nUser defined MLA points A and B are:\n");
      printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
              MLA_points(1,:), MLA_points(2,:));
    else
      ## Print initial user defined MLA_points
      printf ("\nAutomatically registered MLA points A and B are:\n");
      printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
              MLA_points(1,:), MLA_points(2,:));
    endif
    ## Find which optimized MLA points correspond to the user defined pair
    d1 = distancePoints (MLA_points(1,:), MLA_optimal_point1);
    d2 = distancePoints (MLA_points(1,:), MLA_optimal_point2);
    if d1 > d2
      MLA_opt_point_A = MLA_optimal_point2;
      MLA_opt_point_B = MLA_optimal_point1;
    else
      MLA_opt_point_A = MLA_optimal_point1;
      MLA_opt_point_B = MLA_optimal_point2;
    endif
    ## Print optimized MLA_points
    printf ("\nOptimized MLA points A and B are:\n");
    printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
            MLA_opt_point_A, MLA_opt_point_B);
  else                                  # for ulna
    if (! register)
      ## Print initial user defined MLA_points
      printf ("\nUser defined point for orienting Ulna is:\n");
      printf ("x: %f y: %f z: %f\n", MLA_points(1,:));
    else
      ## Print initial user defined MLA_points
      printf ("\nAutomatically registered point for orienting Ulna is:\n");
      printf ("x: %f y: %f z: %f\n", MLA_points(1,:));
    endif
  endif

  ## Save user defined and optimized MLA point back to the Meshlab Point file
  MLP = [MLA_points([1:2],:); MLA_opt_point_A; MLA_opt_point_B];
  write_MeshlabPoints (fullfile (folder, filenamePP), filename, MLP);

  ## Calculate new normals for each final cross section at 20, 35, 50, 65, 80%
  ## All normals should be pointing upwards, i.e. from distal towards proximal
  n1 = (Centroid_1 - Centroid_2) ./ sqrt (sum ((Centroid_1 - Centroid_2) .^2 ));
  n2 = (Centroid_1 - Centroid_3) ./ sqrt (sum ((Centroid_1 - Centroid_3) .^2 ));
  n3 = (Centroid_2 - Centroid_4) ./ sqrt (sum ((Centroid_2 - Centroid_4) .^2 ));
  n4 = (Centroid_3 - Centroid_5) ./ sqrt (sum ((Centroid_3 - Centroid_5) .^2 ));
  n5 = (Centroid_4 - Centroid_5) ./ sqrt (sum ((Centroid_4 - Centroid_5) .^2 ));

  ## Calculate new cross sections for every point in centroid_N by using the
  ## normal of the nearest default centroid (20%, 35%, 50%, 65%, 80%)

  ## Always start with cross section at 20%
  section_1 = meshSection (v, f, Centroid_1, n1);
  [CS_Geometry, SMoA, polyline] = simple_polygon3D (section_1, n1, ...
                                                    CorPlane_normal);
  ## Print results for centroids and cross sectional areas
  printf (strcat (["Cross section at 20%% has an area of %f mm2,"], ...
                  [" perimeter of %f mm \n and centroid coordinates"], ...
                  [" are: x:%f y:%f z:%f\n\n"]), ...
          CS_Geometry.Area, CS_Geometry.Perimeter, CS_Geometry.Centroid);
  ## Arrange to matrices for csv files
  i = 1;
  geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter, ...
                   CS_Geometry.Centroid, n1, CorPlane_normal];
  inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                  SMoA.Imax, SMoA.theta];
  polygon2D(1,i*2-1) = Centroid_N(i,1);
  polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
  polygon3D(1,i*3-2) = Centroid_N(i,1);
  polygon3D([2:length(polyline.poly3D)+1],[i*3-2:i*3]) = polyline.poly3D;

  ## Continue with user defined points
  CENTROIDS = [Centroid_1; Centroid_2; Centroid_3; Centroid_4; Centroid_5];
  NORMALS = [n1; n2; n3; n4; n5];
  for i = 2:size (Centroid_N, 1) - 1
    ## Find nearest default centroid
    dist2centroid = distancePoints (Centroid_N(i,[2:4]), CENTROIDS);
    index = find (dist2centroid == min (dist2centroid));
    ## Calculate cross section with the normal of the nearest centroid
    section = meshSection (v, f, Centroid_N(i,[2:4]), NORMALS(index,:));
    [CS_Geometry, SMoA, polyline] = simple_polygon3D (section, ...
                                    NORMALS(index,:), CorPlane_normal);
    ## Print results for centroids and cross sectional areas
    printf (strcat (["Cross section at %f has an area of %f mm2,"], ...
                    [" perimeter of %f mm \n and centroid coordinates"], ...
                    [" are: x:%f y:%f z:%f\n\n"]), Centroid_N(i,1), ...
            CS_Geometry.Area, CS_Geometry.Perimeter, CS_Geometry.Centroid);
    ## Arrange to matrices for csv files
    geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter, ...
                     CS_Geometry.Centroid, NORMALS(index,:), CorPlane_normal];
    inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                    SMoA.Imax, SMoA.theta];
    polygon2D(1,i*2-1) = Centroid_N(i,1);
    polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
    polygon3D(1,i*3-2) = Centroid_N(i,1);
    polygon3D([2:length(polyline.poly3D)+1],[i*3-2:i*3]) = polyline.poly3D;
  endfor

  ## Always end with cross section at 80%
  section_5 = meshSection (v, f, Centroid_5, n5);
  [CS_Geometry, SMoA, polyline] = simple_polygon3D (section_5, n5, ...
                                                    CorPlane_normal);
  ## Print results for centroids and cross sectional areas
  printf (strcat (["Cross section at 80%% has an area of %f mm2,"], ...
                  [" perimeter of %f mm \n and centroid coordinates"], ...
                  [" are: x:%f y:%f z:%f\n\n"]), ...
          CS_Geometry.Area, CS_Geometry.Perimeter, CS_Geometry.Centroid);
  ## Arrange to matrices for csv files
  i += 1;
  geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter, ...
                   CS_Geometry.Centroid, n5, CorPlane_normal];
  inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                  SMoA.Imax, SMoA.theta];
  polygon2D(1,i*2-1) = Centroid_N(i,1);
  polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
  polygon3D(1,i*3-2) = Centroid_N(i,1);
  polygon3D([2:length(polyline.poly3D)+1],[i*3-2:i*3]) = polyline.poly3D;

  ## Save to files
  starting = "Cgeometry-";
  name = filename([1:length(filename) - 4]);
  endfile = ".csv";
  filename = strcat (starting, name, endfile);
  csvwrite (fullfile (folder, filename), geometry);
  starting = "Cinertia-";
  filename = strcat (starting, name, endfile);
  csvwrite (fullfile (folder, filename), inertia);
  starting = "Cpolyline2D-";
  filename = strcat (starting, name, endfile);
  csvwrite (fullfile (folder, filename), polygon2D);
  starting = "Cpolyline3D-";
  filename = strcat (starting, name, endfile);
  csvwrite (fullfile (folder, filename), polygon3D);

endfunction

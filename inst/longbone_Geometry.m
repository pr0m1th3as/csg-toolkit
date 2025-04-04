## Copyright (C) 2018-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} longbone_Geometry (@var{filename})
## @deftypefnx {csg-toolkit} {} longbone_Geometry (@var{filename}, @var{bones})
## @deftypefnx {csg-toolkit} {} longbone_Geometry (@var{folder}, @var{filename})
## @deftypefnx {csg-toolkit} {} longbone_Geometry (@var{folder}, @var{filename}, @var{bones})
## @deftypefnx {csg-toolkit} {[@var{GEOM}, @var{SMoA}, @var{bone}, @var{EXTRA}, @var{data}] =} longbone_Geometry (@dots{})
##
## This function analyzes the cross-sectional geometry of an intact humerus,
## ulna, femur, or tibia bone at 20%, 35%, 50%, 65% and 80% along the bone's
## maximum length.
##
## @code{longbone_Geometry (@var{filename})} will analyze the 3D model specified
## in @var{filename} provided that it conforms to the Wavefront OBJ file format
## and it is a pure tringular mesh.  The function will automatically determine
## the type of long bone in @var{filename} and register the required initial
## alignment points according to the bone type.  The 3D model is assumed to be
## in @math{mm} units and present in the working directory.
##
## @var{folder}, which must be a character vector, defines the relative or
## absolute path to the directory containing the 3D model in @var{filename}.
## When omitted, the current working directory is assumed.
##
## @var{bones} must be cell array of character vectors specifying one or more
## long bones that should be analyzed.  @code{longbone_Geometry} will only
## analyze the 3D model if it matches one of the bones named in @var{bones},
## unless the @qcode{'Force'} option is also selected, in which case the user
## must specify a single bone such as in @code{@{'Force', 'Humerus'@}}.  Valid
## options are:
##
## @enumerate
## @item @qcode{'Humerus'}
## @item @qcode{'Ulna'}
## @item @qcode{'Femur'}
## @item @qcode{'Tibia'}
## @item @qcode{'All'}
## @item @qcode{'Force'}
## @end enumerate
##
## Note: When @qcode{'Force'} is included, initial alignment points are also
## required and it is the user's responsibility to ensure the appropriate bone
## is being analyzed.  This option is provided so that damaged or heaviliy
## deformed bones, which may turn out as undefined, can be analyzed.
##
## If the user wants to define the initial alignment points, this can be done
## with a side-car Meshlab PickedPoints file with the same base filename and in
## same directory as the 3D model.  In such case, @code{longbone_Geometry} skips
## the automatic initial point registration and uses the first two points in the
## side-car file for optimizing the bone's mediolateral axis.
##
## @enumerate
## @item For the humerus, these points must be positioned anteriorly on the
## trochlea and capitulum at the distal end of the humerus.
## @item For the ulna, only a single initial alignment point is required that
## must be placed at the anterior tip of the coronoid process.
## @item For the femur, these points must be located posteriorly on the
## medial and lateral condyles at the distal end of the femur.
## @item For the tibia, these points must be positioned posteriorly of the
## medial and lateral condyles at the proximal end of the tibia.
## @end enumerate
##
## Note: The function automatically optimizes the local extremal points of the
## mediolateral axis based on the initial points and subsequently uses the
## optimized coronal plane orientation for computing the second moments of area.
## However, the initial alignment point for ulna is not optimized.  The
## optimized alignment points are appended in the existing MPP file or saved in
## newly created, accordingly.
##
## When @code{longbone_Geometry} is called without any output arguments, it
## stores all geometric properties in CSV files as described in the following
## section.  Each OBJ file must explicitly contain a single 3D bone model.
##
## Assuming a 3D model named @qcode{'bone_ID.obj'}, the following files are
## generated:
##
## @multitable @columnfractions .33 .02 .65
## @headitem Filename @tab @tab Description
## @item @qcode{Dgeometry-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{area} (measured in @math{mm^2}), @qcode{perimeter} (measured in
## @math{mm}), @qcode{centroid} (returned as [@var{x},@var{y},@var{z}]
## coordinates in @math{mm} units), as well as the sectioning and orientation
## @qcode{normals} (returned as an [@var{x},@var{y},@var{z}] vectors).  Each row
## of the CSV file corresponds to a different cross section at 20%, 35%, 50%,
## 65%, and 80% along the bone's maximum length.  The aforementioned properties
## for each cross section are stored as a row vector with the respective order
## in columns [:,[2],[3],[4:6],[7:9],[10:12]]. The first column contains the
## ratios of the sectioning points' distance from the proximal end to the bone's
## maximum length.
##
## @item @qcode{Dinertia-bone_ID.csv} @tab @tab It contains the properties of
## @qcode{Ix}, @qcode{Iy} @qcode{Ixy}, @qcode{Imin}, and @qcode{Imax}, all
## measured in @math{mm^4}, as well as @qcode{theta}, measured in degrees.
## Accordingly, they are saved as row vectors [2:7] with the first column
## containing the ratios of the sectioning points' distance from the proximal
## end to the bone's maximum length.
##
## @item @qcode{Dpolyline2D-bone_ID.csv} @tab @tab It contains the 2D
## coordinates (on an arbitrary @var{x},@var{y} local axis) for each cross
## section as an @math{Nx2} matrix, where @math{N} is the number of points for
## each cross-sectional polygon at 20%, 35%, 50%, 65%, and 80% along the bone's
## maximum length.  These polygons are ordered as 5 column duplets,
## i.e. [[1:2],[3:4],[5:6],[7:8],[9:10]].  Polygon coordinates start from the
## second row ([2:end],:), while the first row contains the relevant ratios at
## columns (1, [1,3,5,7,9]).
##
## @item @qcode{Dpolyline3D-bone_ID.csv} @tab @tab It contains the 3D model
## coordinates for each cross section as an @math{Nx3} matrix, where @math{N} is
## the number of points for each cross-sectional polygon at 20%, 35%, 50%, 65%,
## and 80% along the bone's maximum length.  These polygons are ordered as 5
## column triplets, i.e. [[1:3],[4:6],[7:9],[10:12],[13:15]].  Polygon
## coordinates start from the second row ([2:end],:), while the first row
## contains the relevant ratios at columns (1, [1,4,7,10,13]).
## @end multitable
##
## @qcode{[@var{GEOM}, @var{SMoA}, @var{bone}, @var{EXTRA}, @var{data}] =
## longbone_Geometry (@dots{})} may also return up to five output arguments.
## @var{GEOM} and @var{SMoA} are @math{5x1} structure arrays containing the
## cross-sectrional geometric properties and the second moments of area,
## respectively.  Each element of the structure arrays corresponds to 20%, 35%,
## 50%, 65%, and 80% cross sections.  For more information about the contents of
## the returned structure arrays, see @code{simple_polygon3D}.  @var{bone} is a
## character vector containing the name of the bone as inferred by the
## @code{longbone_Registration} function.  If a 3D model is not recognized, then
## @var{bone} is returned as @qcode{'undefined'}, unless the @qcode{'Force'} has
## been opted by the user.
##
## @var{EXTRA} is an additional scalar structure with the following fields:
##
## @enumerate
## @item @qcode{maxDistance}
## @item @qcode{maxd_V1}
## @item @qcode{maxd_V2}
## @item @qcode{MLA_opt_point_A}
## @item @qcode{MLA_opt_point_B}
## @item @qcode{diaphyseal_bending}
## @item @qcode{angle_20_35}
## @item @qcode{angle_35_50}
## @item @qcode{angle_50_65}
## @item @qcode{angle_65_80}
## @item @qcode{ArPerIndex20}
## @item @qcode{ArPerIndex35}
## @item @qcode{ArPerIndex50}
## @item @qcode{ArPerIndex65}
## @item @qcode{ArPerIndex80}
## @end enumerate
##
## @qcode{maxDistance} refers to the bone's maximum distance measurement, as
## defined by the 3D points in @qcode{maxd_V1} and @qcode{maxd_V2}.  Optimized
## mediolateral axis points are returned in the @qcode{MLA_opt_point_*} fields.
## @qcode{diaphyseal_bending} contains the sum of dihedral angles (in degrees)
## between each pair of consecutive cross-sectional planes, which are also
## returned in the @qcode{angle_$$_$$} fields.  @qcode{ArPerIndex$$} the
## Area-Perimeter Index for each consecutive cross section.
##
## @var{data} is a numeric row vector containing all the measurements computed
## by the @code{longbone_Geometry} in tabular form to facilitate further
## processing.  Use the @code{longbone_Measurements} function to obtain a cell
## array of character vectors with the definitions of the measurements in
## @var{data}.
##
## When @code{longbone_Geometry} is called with output arguments, then no CSV
## files are generated.  When called with more than three output arguments, then
## saving the optimized alignment points to an existing or new MPP file is also
## skipped.
##
## @seealso{longbone_CustomGeometry, longbone_FragmentGeometry,
## longbone_Measurements, longbone_Registration, simple_polygon3D}
## @end deftypefn

function [varargout] = longbone_Geometry (varargin)

  ## Check input output
  if (nargin < 1 || nargin > 3 || nargout > 5)
    print_usage;
  endif

  ## Only 1 arg: filename
  if (nargin == 1 && ischar (varargin{1}))
    folder = '';
    filename = varargin{1};
    bone = {'All'};
  elseif (nargin == 1)
    error ("longbone_Geometry: FILENAME must be a string.");
  endif

  ## Only 2 args: filename and bone or folder and filename
  if (nargin == 2 && ischar (varargin{1}) && ischar (varargin{2}))
    folder = varargin{1};
    filename = varargin{2};
    bone = {'All'};
  elseif (nargin == 2 && ischar (varargin{1}) && iscellstr (varargin{2}))
    folder = '';
    filename = varargin{1};
    bone = varargin{2};
  elseif (nargin == 2)
    if (! ischar (varargin{1}) && ischar (varargin{2}))
      error ("longbone_Geometry: FOLDER must be a string.");
    endif
    if (! ischar (varargin{1}) && iscellstr (varargin{2}))
      error ("longbone_Geometry: FILENAME must be a string.");
    endif
    if (ischar (varargin{1}) && ! (ischar (varargin{2}) ||
                                   iscellstr (varargin{2})))
      error (strcat ("longbone_Geometry: second argument must be either", ...
                     " a string for FILENAME or a cell array of strings", ...
                     " for selected BONES."));
    endif
  endif

  ## Only 3 args: folder, filename, and bone
  if (nargin == 3 && ischar (varargin{1}) && ischar (varargin{2})
                  && iscell (varargin{3}))
    folder = varargin{1};
    filename = varargin{2};
    bone = varargin{3};
  elseif (nargin == 3)
    if (! ischar (varargin{1}))
      error ("longbone_Geometry: FOLDER must be a string.");
    endif
    if (! ischar (varargin{2}))
      error ("longbone_Geometry: FILENAME must be a string.");
    endif
    if (! iscellstr (varargin{3}))
      error ("longbone_Geometry: BONES must be a cell array of strings.");
    endif
  endif

  ## Check filename has a valid .obj extension
  if (! (strcmpi (filename([end - 3:end]), '.obj')))
    error ("longbone_Geometry: 3D model must be in OBJ file format.");
  endif

  ## Check bone contains valid options
  valid_names = {'Humerus', 'Ulna', 'Femur', 'Tibia', 'All', 'Force'};
  fcn =  @(x) any (strcmpi (x, valid_names));
  if (! all (cellfun (fcn, bone)))
    error ("longbone_Geometry: invalid bone selection.");
  endif
  if (any (strcmpi (bone, 'Force')))
    do_force = true;
    if (numel (bone) != 2 || any (strcmpi (bone, 'All')))
      error (strcat ("longbone_Geometry: only one bone name can", ...
                     " be defined when using the 'Force' option."));
    endif
    bone(strcmpi (bone, 'Force')) = [];
  else
    do_force = false;
  endif

  ## Prepare return list
  varargout = cell (1, nargout);

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
    register = false;
  elseif (! do_force)
    [bonesel, MLA_points] = longbone_Registration (v, f);
    find_bone = false;
    register = true;
  else
    error (strcat ("longbone_Geometry: cannot use 'Force' option if", ...
                   " a side-car Meshlab PP file is not be available."));
  endif

  ## Find bone if necessary
  if (find_bone && ! do_force)
    bonesel = longbone_Registration (v, f);
  endif

  ## Check bone selection
  if (! do_force)
    if (any (strcmpi (bone, "All")))
      bone = bonesel;
    else
      if (numel (bone) == 1 && ! strcmpi (bone, bonesel))
        printf ("Model %s is not a %s.\n", filename, bone{1});
        varargout{3} = bonesel;
        return;
      endif
      if (numel (bone) == 2 && ! (any (strcmpi (bone, bonesel))))
        printf ("Model %s is neither a %s nor a %s.\n", ...
               filename, bone{1}, bone{2});
        varargout{3} = bonesel;
        return;
      endif
      if (numel (bone) == 3 && ! (any (strcmpi (bone, bonesel))))
        printf ("Model %s is not a %s, a %s or a %s.\n", ...
                filename, bone{1}, bone{2}, bone{3});
        varargout{3} = bonesel;
        return;
      endif
      bone = bonesel;
    endif
  endif

  ## Check valid bone selection
  if (! (strcmpi (bone, 'Humerus') || strcmpi (bone, 'Ulna') ||
         strcmpi (bone, 'Femur') || strcmpi (bone, 'Tibia')))
    printf ("Bone should be either a Humerus, an Ulna, a Femur, or a Tibia.\n");
    printf ("3D model in %s is undefined.\n", filename);
    varargout{3} = bonesel;
    return;
  endif

  ## Find the maximum distance of the bone
  [maxDistance, maxd_V1, maxd_V2] = longbone_maxDistance (v);

  ## Print bone's max distance
  page_screen_output (0);
  page_output_immediately (1);
  printf ("\n%s in %s has a maximum distance of %0.1f mm\n", ...
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
  GEOM(1) = simple_polygon3D (plane_1, normal);
  GEOM(2) = simple_polygon3D (plane_2, normal);
  GEOM(3) = simple_polygon3D (plane_3, normal);
  GEOM(4) = simple_polygon3D (plane_4, normal);
  GEOM(5) = simple_polygon3D (plane_5, normal);

  ## Store the centroids of the initial cross-sectional areas
  Centroid_1 = GEOM(1).Centroid;
  Centroid_2 = GEOM(2).Centroid;
  Centroid_3 = GEOM(3).Centroid;
  Centroid_4 = GEOM(4).Centroid;
  Centroid_5 = GEOM(5).Centroid;

  ## Check the centroids' locations from proximal (20%) to distal (80%) and
  ## if required reverse their order to comply with this standard.
  proximal = distancePoints (Centroid_1, MLA_points(1,:));
  distal = distancePoints (Centroid_5, MLA_points(1,:));
  if ((strcmpi (bone, 'Humerus') || strcmpi (bone, 'Femur')) &&
      (proximal < distal))
    Centroid_1 = GEOM(5).Centroid;
    Centroid_2 = GEOM(4).Centroid;
    Centroid_3 = GEOM(3).Centroid;
    Centroid_4 = GEOM(2).Centroid;
    Centroid_5 = GEOM(1).Centroid;
  elseif ((strcmpi (bone, 'Tibia') || strcmpi (bone, 'Ulna')) &&
          (proximal > distal))
    Centroid_1 = GEOM(5).Centroid;
    Centroid_2 = GEOM(4).Centroid;
    Centroid_3 = GEOM(3).Centroid;
    Centroid_4 = GEOM(2).Centroid;
    Centroid_5 = GEOM(1).Centroid;
  endif

  ## Calculate normals of anatomical planes
  if (! strcmpi (bone, 'Ulna'))         # for humerus, femur, and tibia
    ## Calculate the mediolateral axis vector
    MLA_vector = MLA_points(1,:) - MLA_points(2,:);
    ## Calculate longitudinal axis vector
    PDA_vector = Centroid_1 - Centroid_5;
    ## Calculate the unit normal of the transverse plane facing upwards
    TransPlane_normal = PDA_vector ./ sqrt (sum (PDA_vector .^ 2));
    ## Calculate the vector of the coronal plane
    CorPlane_normal = cross(PDA_vector, MLA_vector);
    ## Normalize the vector of the coronal plane
    CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
  else                                  # for ulna
    ## No optimization of user defined points for ulna
    ## Calculate the projection of user defined point to the longitudinal axis
    C5_MLA = MLA_points(1,:) - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    point_proj = Centroid_1 + (dot (C5_MLA, C5_C1) / dot (C5_C1, C5_C1)) * C5_C1;
    ## Calculate the coronal plane unit vector by normalizing the vector from
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
  if (strcmpi (bone, 'Humerus'))
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
  elseif (strcmpi (bone, 'Femur'))
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
  elseif (strcmpi (bone, 'Tibia'))
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
  if (! strcmpi (bone, 'Ulna'))         # for humerus, femur, and tibia
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
  if (nargout < 4)
    write_MeshlabPoints (fullfile (folder, filenamePP), filename, MLP);
  endif

  ## Calculate new normals for each final cross section at 20, 35, 50, 65, 80%
  ## All normals should be pointing upwards, i.e. from distal towards proximal
  n1 = (Centroid_1 - Centroid_2) ./ sqrt (sum ((Centroid_1 - Centroid_2) .^2 ));
  n2 = (Centroid_1 - Centroid_3) ./ sqrt (sum ((Centroid_1 - Centroid_3) .^2 ));
  n3 = (Centroid_2 - Centroid_4) ./ sqrt (sum ((Centroid_2 - Centroid_4) .^2 ));
  n4 = (Centroid_3 - Centroid_5) ./ sqrt (sum ((Centroid_3 - Centroid_5) .^2 ));
  n5 = (Centroid_4 - Centroid_5) ./ sqrt (sum ((Centroid_4 - Centroid_5) .^2 ));

  ## Recalculate the sectioning points for each cross section
  section_1 = meshSection (v, f, Centroid_1, n1);
  section_2 = meshSection (v, f, Centroid_2, n2);
  section_3 = meshSection (v, f, Centroid_3, n3);
  section_4 = meshSection (v, f, Centroid_4, n4);
  section_5 = meshSection (v, f, Centroid_5, n5);

  ## Recalculate CSG properties for each cross section
  [GEOM(1), SMoA(1), polyline(1)] = simple_polygon3D (section_1, n1, ...
                                                             CorPlane_normal);
  [GEOM(2), SMoA(2), polyline(2)] = simple_polygon3D (section_2, n2, ...
                                                             CorPlane_normal);
  [GEOM(3), SMoA(3), polyline(3)] = simple_polygon3D (section_3, n3, ...
                                                             CorPlane_normal);
  [GEOM(4), SMoA(4), polyline(4)] = simple_polygon3D (section_4, n4, ...
                                                             CorPlane_normal);
  [GEOM(5), SMoA(5), polyline(5)] = simple_polygon3D (section_5, n5, ...
                                                             CorPlane_normal);

  ## Arrange to matrices for csv files
  cs = [20, 35, 50, 65, 80];
  ns = [n1; n2; n3; n4; n5];
  for i = 1:5
    ## Print results for centroids and cross sectional areas
    printf (strcat (["\nCross section at %d%% has an area of %f mm^2,"], ...
                    [" perimeter of %f mm \n and centroid coordinates"], ...
                    [" are: x:%f y:%f z:%f\n"]), ...
            cs(i), GEOM(i).Area, GEOM(i).Perimeter, GEOM(i).Centroid);

    ## Arrange to matrices for csv files
    geometry(i,:) = [cs(i), GEOM(i).Area, GEOM(i).Perimeter, ...
                     GEOM(i).Centroid, ns(i,:), CorPlane_normal];
    inertia(i,:) = [cs(i), SMoA(i).Ix, SMoA(i).Iy, SMoA(i).Ixy, ...
                    SMoA(i).Imin, SMoA(i).Imax, SMoA(i).theta];
    polygon2D(1,i*2-1) = cs(i) / 100;
    polygon2D([2:length(polyline(i).poly2D)+1],[i*2-1:i*2]) = ...
                                                            polyline(i).poly2D;
    polygon3D(1,i*3-2) = cs(i) / 100;
    polygon3D([2:length(polyline(i).poly3D)+1],[i*3-2:i*3]) = ...
                                                            polyline(i).poly3D;
  endfor

  ## Compute dihedral angles and their total sum and append them into
  ## a structure along with maxDistance, maxd_V1, maxd_V2, and optimized
  ## points MLA_opt_point_A and MLA_opt_point_B
  if (nargout > 3)
    extra.maxDistance = maxDistance;
    extra.maxd_V1 = maxd_V1;
    extra.maxd_V2 = maxd_V2;
    extra.MLA_opt_point_A = MLA_opt_point_A;
    extra.MLA_opt_point_B = MLA_opt_point_B;
    extra.diaphyseal_bending = 0;
    extra.angle_20_35 = atan2d (norm (cross (n1, n2)), dot (n1, n2));
    extra.diaphyseal_bending += extra.angle_20_35;
    extra.angle_35_50 = atan2d (norm (cross (n2, n3)), dot (n2, n3));
    extra.diaphyseal_bending += extra.angle_35_50;
    extra.angle_50_65 = atan2d (norm (cross (n3, n4)), dot (n3, n4));
    extra.diaphyseal_bending += extra.angle_50_65;
    extra.angle_65_80 = atan2d (norm (cross (n4, n5)), dot (n4, n5));
    extra.diaphyseal_bending += extra.angle_65_80;
    extra.ArPerIndex20 = ((GEOM(1).Area) * 4 * pi) ./ ((GEOM(1).Perimeter) .^ 2);
    extra.ArPerIndex35 = ((GEOM(2).Area) * 4 * pi) ./ ((GEOM(2).Perimeter) .^ 2);
    extra.ArPerIndex50 = ((GEOM(3).Area) * 4 * pi) ./ ((GEOM(3).Perimeter) .^ 2);
    extra.ArPerIndex65 = ((GEOM(4).Area) * 4 * pi) ./ ((GEOM(4).Perimeter) .^ 2);
    extra.ArPerIndex80 = ((GEOM(5).Area) * 4 * pi) ./ ((GEOM(5).Perimeter) .^ 2);
  endif

  ## Aggregate measurements in numeric arrays
  if (nargout > 4)
    data = [maxDistance, ...
            GEOM(1).Area, GEOM(1).Perimeter, extra.ArPerIndex20, SMoA(1).Ix, ...
            SMoA(1).Iy, SMoA(1).Ixy, SMoA(1).Ix/SMoA(1).Iy, SMoA(1).Imin, ...
            SMoA(1).Imax, SMoA(1).Imax/SMoA(1).Imin, SMoA(1).theta, ...
            extra.angle_20_35, ...
            GEOM(2).Area, GEOM(2).Perimeter, extra.ArPerIndex35, SMoA(2).Ix, ...
            SMoA(2).Iy, SMoA(2).Ixy, SMoA(2).Ix/SMoA(2).Iy, SMoA(2).Imin, ...
            SMoA(2).Imax, SMoA(2).Imax/SMoA(2).Imin, SMoA(2).theta, ...
            extra.angle_35_50, ...
            GEOM(3).Area, GEOM(3).Perimeter, extra.ArPerIndex50, SMoA(3).Ix, ...
            SMoA(3).Iy, SMoA(3).Ixy, SMoA(3).Ix/SMoA(3).Iy, SMoA(3).Imin, ...
            SMoA(3).Imax, SMoA(3).Imax/SMoA(3).Imin, SMoA(3).theta, ...
            extra.angle_50_65, ...
            GEOM(4).Area, GEOM(4).Perimeter, extra.ArPerIndex65, SMoA(4).Ix, ...
            SMoA(4).Iy, SMoA(4).Ixy, SMoA(4).Ix/SMoA(4).Iy, SMoA(4).Imin, ...
            SMoA(4).Imax, SMoA(4).Imax/SMoA(4).Imin, SMoA(4).theta, ...
            extra.angle_65_80, ...
            GEOM(5).Area, GEOM(5).Perimeter, extra.ArPerIndex80, SMoA(5).Ix, ...
            SMoA(5).Iy, SMoA(5).Ixy, SMoA(5).Ix/SMoA(5).Iy, SMoA(5).Imin, ...
            SMoA(5).Imax, SMoA(5).Imax/SMoA(5).Imin, SMoA(5).theta, ...
            extra.diaphyseal_bending];
  endif

  if (nargout == 0)
    ## Save to files
    starting = 'Dgeometry-';
    name = filename([1:length(filename) - 4]);
    endfile = ".csv";
    filename = strcat (starting, name, endfile);
    csvwrite (filename, geometry);
    starting = 'Dinertia-';
    filename = strcat (starting, name, endfile);
    csvwrite (filename, inertia);
    starting = 'Dpolyline2D-';
    filename = strcat (starting, name, endfile);
    csvwrite (filename, polygon2D);
    starting = 'Dpolyline3D-';
    filename = strcat (starting, name, endfile);
    csvwrite (filename, polygon3D);
  elseif (nargout == 1)
    varargout{1} = GEOM;
  elseif (nargout == 2)
    varargout{1} = GEOM;
    varargout{2} = SMoA;
  elseif (nargout == 3)
    varargout{1} = GEOM;
    varargout{2} = SMoA;
    varargout{3} = bone;
  elseif (nargout == 4)
    varargout{1} = GEOM;
    varargout{2} = SMoA;
    varargout{3} = bone;
    varargout{4} = extra;
  elseif (nargout == 5)
    varargout{1} = GEOM;
    varargout{2} = SMoA;
    varargout{3} = bone;
    varargout{4} = extra;
    varargout{5} = data;
  endif

endfunction

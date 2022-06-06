## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn{Function} longbone_CustomGeometry (@var{folder}, @var{filename}, @var{bones}, @var{points})
## @deftypefnx{Function} longbone_CustomGeometry (@var{folder}, @var{filename}, @var{points})
## @deftypefnx{Function} longbone_CustomGeometry (@var{filename}, @var{bones}, @var{points})
## @deftypefnx{Function} longbone_CustomGeometry (@var{filename}, @var{points})
##
## This function analyzes the cross-sectional geometry of a humerus, ulna, femur
## or tibia along their longitudinal axis at custom intervals specified in
## @var{points}, which can be either a string or a numerical vector.  The string
## must specify a Meshlab Point file containing an arbitrary number of points on
## the bone's surface, each one of which specify a cross-sectional plane along
## the bones longitudinal axis.  Alternatively, the numerical vector specifies
## an arbitrary number of cross-sectional plane as a percentage of the bone's
## maximum length in the range [0.1 0.9].  Values beyond this range are ignored.
##
## The function works in a similar fashion to the 'longbone_Geometry' function
## regarding the rest of the input arguments @var{folder}, @var{filename}, and
## @var{bones}, whereas it does not return any output arguments.  Instead, it
## saves all geometric properties in custom .csv files as desccribed below.
##
## For example, for a triangular mesh named as 'bone_ID.obj', the following
## files are produced:
##
## @multitable {Filename} {Description} @columnfractions .35 .65
## @item Cgeometry-bone_ID.csv @tab containing the properties of area, perimeter,
## centroid, slicing and orientation normals for each cross section stored as a
## row vector in the order (:,[2],[3],[4:6],[7:9],[10:12]).  The first collumn
## contains the corresponding cross-sectioning points.  The first and last rows
## contain the default cross sections at 20% and 80% respectively.
## The intermediate rows follow the order at which the cross-secional points are
## stored in the Meshlab Point file or the values given in the numerical vector.
##
## @item Cinertia-bone_ID.csv @tab containing the properties of Ix, Iy, Ixy,
## Imin,Imax and theta angle for each cross section as a row vector (:,[2:7])
## with the first collumn (:,1) containing the cross-sectioning points as
## described in Cgeometry-bone_ID.csv.
##
## @item Cpolyline2D-bone_ID.csv @tab containing the 2D coordinates for each
## cross section as an Nx2 matrix, where N is the number of points for each
## polygon.  The polygons are ordered as (:,[1:2],...,[end-1:end]) starting from
## from the second row ([2:end],:), while the first row (1,[1,3,...,end-1])
## contains the corresponding cross-sectioning points.  The first and
## penultimate columns contain the default cross sections at 20% and 80%
## respectively.  The intermediate odd collumns follow the order at which the
## cross-secional points are stored in the Meshlab Point file or the values
## given in the numerical vector.
## @end multitable
##
## For each analyzed OBJ file, the initial alignment points are either
## read from the corresponding Meshlab Point file (e.g. 'bone_ID.pp'), if
## present in the same folder with the OBJ, or they are automatically registered
## with the 'longbone_Registration' function and they are saved in newly created
## Meshlab Point files. Each OBJ file must explicitly contain a single bone.
##
## All bone models stored as OBJ files must be strictly triangular meshes and
## their coordinate umits are assumed to be in mm. Only intact femur, tibia,
## humerus, and ulna bones can be processed. Side does not matter.
## @seealso{longbone_Geometry, visualize_CrossSections}
## @end deftypefn

function longbone_CustomGeometry (varargin)
  ## check input variables and parse them accordingly
  if nargin < 2 || nargin > 4
    printf ("Invalid number of input arguments\n");
    return;
  endif
  ## 1st arg is filename
  if (nargin == 2 && ischar (varargin{1}(:)'))
    folder = "";
    filename = varargin{1}(:)';
    ## check filename has a valid .obj extension
    if ! (strcmpi (filename([end - 3:end]), ".obj"))
      printf ("Model must be in OBJ file format\n");
      return;
    endif
    find_bone = true;
  elseif (nargin == 1 && ! ischar (varargin{1}(:)'))
    printf("Filename must be a string\n");
    return;
  endif
  ## 1st and 2nd args are either filename and bone or folder and filename
  if (nargin == 3 && ischar (varargin{1}(:)') && ischar (varargin{2}(:)'))
    bone = varargin{2}(:)';
    ## filename and bone (string)
    if (strcmp (bone, "Humerus") || strcmp (bone, "Ulna") || ...
        strcmp (bone, "Femur") || strcmp (bone, "Tibia"))
      folder = "";
      filename = varargin{1}(:)';
      find_bone = false;
    ## folder and filename
    else
      folder = varargin{1}(:)';
      filename = varargin{2}(:)';
      find_bone = true;
      bone = "All";
    endif
  endif
  ## 1st and 2nd args are filename and bone (numeric)
  if (nargin == 3 && ischar (varargin{1}(:)') && isnumeric (varargin{2}))
    folder = "";
    filename = varargin{1}(:)';
    bonesnum = varargin{2};
  endif
  ## 1st, 2nd and 3rd args are folder, filename, and bone (numeric)
  if (nargin == 4 && ischar (varargin{1}(:)') && ischar (varargin{2}(:)') ...
      && isnumeric (varargin{3}))
    folder = varargin{1}(:)';
    filename = varargin{2}(:)';
    bonesnum = varargin{3};
  endif
  ## check last argument being a string or a numerical vector
  if ischar (varargin{end}(:)')
    filenamePP = varargin{end}(:)';
    if (exist (fullfile (folder, filenamePP)) == 2)
      ## load Meshlab points for custom sectioning points from pp file
      CSP_data = read_MeshlabPoints (fullfile (folder, filenamePP));
      CSP_data(:,1) = [];
      CSP_points = true;
    else
      printf("Meshlab Point file for custom sectioning points does not exist\n");
      return;
    endif
  elseif isnumeric (varargin{end})
    CSP_data = varargin{end};
    CSP_points = false;
    if ! (size (CSP_data, 1) == 1)
      printf("Custom sectioning points vector must be a row vector\n");
      return;
    endif
  else
    printf("Last argument must be string or a numerical vector\n");
  endif
  
  ## check numeric argument for bone selection
  if (exist ("bonesnum") == 1)
    find_bone = true;
    if (length (bonesnum) == 1)
      switch (bonesnum)
        case 1
          bone = "Humerus";
        case 2
          bone = "Ulna";
        case 3
          bone = "Femur";
        case 4
          bone = "Tibia";
        case 5
          bone = "All";
        otherwise
          printf ("Invalid numeric argument for bone selection\n");
          return;
       endswitch
    endif
    if (length (bonesnum) > 1)
      if (any (bonesnum < 1) || any (bonesnum > 5))
        printf ("Invalid numeric argument for bone selection\n");
        return;
      endif
      i = 0;
      if (any (bonesnum==1))
        i++;
        bones(i) = {"Humerus"};
      endif
      if (any (bonesnum==2))
        i++;
        bones(i) = {"Ulna"};
      endif
      if (any (bonesnum==3))
        i++;
        bones(i) = {"Femur"};
      endif
      if (any (bonesnum==4))
        i++;
        bones(i) = {"Tibia"};
      endif
      if (i == 4 || any (bonesnum==5))
        bone = "All";
        bonesnum = 0;     % so that length equals 1
      endif
    endif
  endif
    
  ## check if orientation points are available in corresponding pp file
  ## otherwise automatically register new ones on the bone surface
  ## for Meshlab point files
  filenamePP = filename([1:length(filename) - 4]);
  extension = ".pp";
  filenamePP = strcat (filenamePP, extension);
  if (exist (fullfile (folder, filenamePP)) == 2)
    ## load Meshlab points for mediolateral axis from pp file
    MLA_points = read_MeshlabPoints (fullfile (folder, filenamePP));
    MLA_points(:,1) = [];
    register = false;
  else
    register = true;
  endif
  
  ## load vertices and faces of triangular mesh from obj file
  [v,f] = readObj (fullfile (folder, filename));
  ## keep only vertex coordinates (just in case color information is present)
  v = v(:,[1:3]);
  ## find bone and register points as appropriate
  if (find_bone && ! register)
    bonesel = longbone_Registration (v, f);
  elseif (find_bone && register)
    [bonesel, MLA_points] = longbone_Registration (v, f);
  elseif (! find_bone && register)
    [nobone, MLA_points] = longbone_Registration (v, f);
    bonesel = bone;
  elseif (! find_bone && !register)
    [nobone, MLA_points] = longbone_Registration (v, f);
    bonesel = bone;
  endif
  
  ## checking bone selection
  if (nargin == 2)
    bone = bonesel;
  endif
  if (strcmp (bone, "All"))
    bone = bonesel;
  endif
  if (exist ("bonesnum") == 1)
    if (length (bonesnum) == 1 && ! strcmp (bone, bonesel))
      printf ("Model %s is not a %s\n", filename, bone);
      return;
    endif
    if (length(bonesnum) == 2 && ! (strcmp (bones(1), bonesel) ...
        || strcmp (bones(2), bonesel)))
      printf("Model %s is neither a %s nor a %s\n", filename, bones{1}, bones{2});
      return;
    endif
    if (length(bonesnum) == 2 && ! (strcmp (bones(1), bonesel) ...
        || strcmp (bones(2), bonesel) || strcmp (bones(3), bonesel)))
      printf ("Model %s is not a %s, a %s or a %s\n", ...
              filename, bones{1}, bones{2}, bones{3});
      return;
    else
      bone = bonesel;
    endif
  endif
  ## check if bone is properly determined
  if ! (strcmp (bone, "Humerus") || strcmp (bone, "Ulna") || ...
        strcmp (bone, "Femur") || strcmp (bone, "Tibia"))
    printf("Bone should be a Humerus, an Ulna, a Femur or a Tibia\n");
    return;
  else
    clear CS_Geometry SMoA polyline
  endif
  
  ## find the maximum distance of the bone
  [maxDistance, maxd_V1, maxd_V2] = longbone_maxDistance(v);
  ## calculate the normal vector of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt (sum ((maxd_V2 - maxd_V1) .^ 2));
  ## calculate planes at 20 and 80% of the maximum length of the bone
  point_1 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.20;
  point_2 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.35;
  point_3 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.50;
  point_4 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.65;
  point_5 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.80;
  ## calculate the sectioning points for each plane
  plane_1 = meshSection (v, f, point_1, normal);
  plane_2 = meshSection (v, f, point_2, normal);
  plane_3 = meshSection (v, f, point_3, normal);
  plane_4 = meshSection (v, f, point_4, normal);
  plane_5 = meshSection (v, f, point_5, normal);
  ## for each cross section calculate the 3D coordinates of its centroid and area
  CS_Geometry(1) = simple_polygon3D (plane_1, normal);
  CS_Geometry(2) = simple_polygon3D (plane_2, normal);
  CS_Geometry(3) = simple_polygon3D (plane_3, normal);
  CS_Geometry(4) = simple_polygon3D (plane_4, normal);
  CS_Geometry(5) = simple_polygon3D (plane_5, normal);
  ## store the centroids of the initial cross-sectional areas
  Centroid_1 = CS_Geometry(1).Centroid;
  Centroid_2 = CS_Geometry(2).Centroid;
  Centroid_3 = CS_Geometry(3).Centroid;
  Centroid_4 = CS_Geometry(4).Centroid;
  Centroid_5 = CS_Geometry(5).Centroid;
  ## check the centroids' locations from proximal (20%) to distal (80%) and 
  ## if required reverse their order to comply with this standard.
  proximal = distancePoints (Centroid_1, MLA_points(1,:));
  distal = distancePoints (Centroid_5, MLA_points(1,:));
  prox_Edge = maxd_V2;
  dist_Edge = maxd_V1;
  if (strcmp (bone, "Humerus") || strcmp (bone, "Femur")) && (proximal < distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
    prox_Edge = maxd_V1;
    dist_Edge = maxd_V2;
  elseif (strcmp (bone, "Tibia") || strcmp (bone, "Ulna")) && (proximal > distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
    prox_Edge = maxd_V1;
    dist_Edge = maxd_V2;
  endif
  ## check for custom sectioning points from pp file
  if CSP_points
    ## calculate the projection of the user defined points to the longitudinal
    ## axis defined by maxd_V1 and maxd_V2 and keep the point that lie within
    ## the range of 10% to 90% of maximum length
    long_axis = prox_Edge - dist_Edge;
    p = 1;
    ## first centroid at 20%
    Centroid_N(p,:) = [0.2, Centroid_1];
    for i = 1: size (CSP_data, 1)
      dist_CSP = CSP_data(i,:) - dist_Edge;
      point_proj = dist_Edge + (dot (dist_CSP, long_axis) / ...
                                dot (long_axis, long_axis)) * long_axis;
      ## check range
      dist = distancePoints (dist_Edge, point_proj);
      if (dist < maxDistance * 0.9 && dist > maxDistance * 0.1)
        p += 1;
        ratio = dist / maxDistance;
        ratio = round (ratio * 1000) / 1000;
        Centroid_N(p,:) = [ratio, point_proj];
      endif
    endfor
    ## last centroid at 80%
    p += 1;
    Centroid_N(p,:) = [0.8, Centroid_5];
  else
    ## keep the point that lie within the range of 10% to 90% of maximum length
    p = 1;
    ## first centroid at 20%
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
    ## last centroid at 80%
    p += 1;
    Centroid_N(p,:) = [0.8, Centroid_5];
  endif
    
  ## calculate normals of anatomical planes
  ## for humerus, femur, and tibia
  if ! strcmp (bone, "Ulna")
    ## calculate the mediolateral axis vector
    MLA_vector = MLA_points(1,:) - MLA_points(2,:);
    ## calculate longitudinal axis vector
    PDA_vector = Centroid_1 - Centroid_5;
    ## calculate the unit normal of the transverse plane facing upwards
    TransPlane_normal = PDA_vector ./ sqrt (sum (PDA_vector .^ 2));
    ## calculate the vector of the coronal plane
    CorPlane_normal = cross(PDA_vector, MLA_vector);
    ## normalize the vector of the coronal plane
    CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
  else
    ## no optimization of user defined points for ulna
    ## calculate the projection of user defined point to the longitudinal axis
    C5_MLA = MLA_points(1,:) - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    point_proj = Centroid_1 + (dot (C5_MLA, C5_C1) / dot (C5_C1, C5_C1)) * C5_C1;
    ## recalculate the coronal plane unit vector by normalizing the vector from
    ## midpoint_proj to MLA_midpoint
    CorPlane_normal = MLA_points(1,:) - point_proj;
    CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
    ## calculate the unit normal of the transverse plane facing upwards
    TransPlane_normal = C5_C1 ./ sqrt (sum (C5_C1 .^ 2));
    ## calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    MLA_opt_point_A = MLA_points(1,:);
    MLA_opt_point_B = MLA_points(2,:);
  endif
  
  ## find if coronal plane normal points in the right direction (towards the
  ## front) by checking the dot product of the normal with the vector between 
  ## the first point of MLA vector and the nearest centroid. If not, reverse
  ## the normal so that it points to the front
  
  ## for humerus
  if strcmp (bone, "Humerus")
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if sum (MLP1_C5 .* CorPlane_normal) > 0     ## MLA_points should be in front
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## find the points that lie on or below the transverse plane at C5
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find (TP_dotP <= 0),:);
    ## find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 1 and 5 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + ...
                    (dot (C1_midpoint, C1_C5) / dot (C1_C5, C1_C5)) * C1_C5;
    ## recalculate the coronal plane unit vector by normalizing the vector from
    ## midpoint_proj to MLA_midpoint
    CorPlane_normal = MLA_midpoint - midpoint_proj;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    ## calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## find the points that lie in front of the coronal plane at midpoint_proj
    CP_dotP = sum (CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_front = v_below(find (CP_dotP > 0),:);
    ## find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_below_front - MLA_midpoint), 2);
    v_side1 = v_below_front(find (ML_dotP <= 0),:);
    v_side2 = v_below_front(find (ML_dotP >= 0),:);
    ## initialize 5 iterations for finding the optimal points
    for i=1:5
      ## calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj), 2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj), 2);
      ## find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane -optimize local extrema
      [d1,i1] = max (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = max (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt( sum (CorPlane_normal .^ 2));
      ## check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if sum (MLP1_C5 .* CorPlane_normal) > 0   ## MLA_points should be in front
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
    
  ## for femur
  elseif strcmp (bone, "Femur")
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if sum (MLP1_C5 .* CorPlane_normal) < 0       ## MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## find the points that lie on or below the transverse plane at C5
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find (TP_dotP <= 0),:);
    ## find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 1 and 5 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + ...
                    (dot (C1_midpoint, C1_C5) / dot (C1_C5, C1_C5)) * C1_C5;
    ## recalculate the coronal plane unit vector by normalizing the vector from
    ## MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    ## calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## find the points that lie behind the coronal plane of the midpoint_proj
    CP_dotP = sum (CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_behind = v_below(find (CP_dotP < 0),:);
    ## find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_below_behind - MLA_midpoint), 2);
    v_side1 = v_below_behind(find (ML_dotP <= 0),:);
    v_side2 = v_below_behind(find (ML_dotP >= 0),:);
    ## initialize 5 iterations for finding the optimal points
    for i=1:5
      ## calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj), 2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj), 2);
      ## find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane -optimize local extrema
      [d1,i1] = min (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = min (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt (sum (CorPlane_normal .^ 2));
      ## check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if sum (MLP1_C5 .* CorPlane_normal) < 0     ## MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  
  ## for tibia
  elseif strcmp (bone, "Tibia")
    MLP1_C1 = Centroid_1 - MLA_points(1,:);
    if sum (MLP1_C1 .* CorPlane_normal) < 0       ## MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    ## find the points that lie on or above the transverse plane at C1
    TP_dotP = sum (TransPlane_normal .* (v - Centroid_1), 2);
    v_above = v(find (TP_dotP >= 0),:);
    ## find the midpoint between the user defined points and project it to the
    ## vertical vector defined by centroids 5 and 1 to find the closest point
    ## along the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C5_midpoint = MLA_midpoint - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    midpoint_proj = Centroid_5 + ...
                    (dot (C5_midpoint, C5_C1) / dot (C5_C1, C5_C1)) * C5_C1;
    ## recalculate the coronal plane unit vector by normalizing the vector from
    ## MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal .^ 2));
    ## calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross (CorPlane_normal, TransPlane_normal);
    ## find the points that lie behind the coronal plane of the midpoint_pro
    CP_dotP = sum (CorPlane_normal .* (v_above - midpoint_proj), 2);
    v_above_behind = v_above(find (CP_dotP < 0),:);
    ## find the points that lie on one or the other side
    ## of the sagital plane of the MLA_midpoint
    ML_dotP = sum (SagPlane_normal .* (v_above_behind - MLA_midpoint), 2);
    v_side1 = v_above_behind(find (ML_dotP <= 0),:);
    v_side2 = v_above_behind(find (ML_dotP >= 0),:);
    ## initialize 5 iterations for finding the optimal points
    for i=1:5
      ## calculate the distance for each point of each side from the coronal
      ## plane passing through the midpoint_proj point
      d_side1 = sum (CorPlane_normal .* (v_side1 - midpoint_proj),2);
      d_side2 = sum (CorPlane_normal .* (v_side2 - midpoint_proj),2);
      ## find new points for optimized mediolateral axis vector that have the
      ## maximum distance from the initial coronal plane -optimize local extrema
      [d1,i1] = min (d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2] = min (d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      ## calculate the optimized normal of the coronal plane
      CorPlane_normal = cross (PDA_vector, MLA_vector);
      ## normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal .^ 2));
      ## check the optimized CorPlane_normal points towards the right direction
      MLP1_C1 = Centroid_1 - MLA_optimal_point1;
      if sum (MLP1_C1 .* CorPlane_normal) < 0     ## MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  endif
  
  ## flush the screen output to display the results during iterations
  page_screen_output(0);
  page_output_immediately(1);
  if ! strcmp (bone, "Ulna")
    if ! register
      ## print initial user defined MLA_points
      printf ("\nUser defined MLA points A and B are:\n");
      printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
              MLA_points(1,:), MLA_points(2,:));
    else
      ## print initial user defined MLA_points
      printf ("\nAutomatically registered MLA points A and B are:\n");
      printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
              MLA_points(1,:), MLA_points(2,:));
    endif
    ## find which optimized MLA points correspond to the user defined pair
    d1 = distancePoints (MLA_points(1,:), MLA_optimal_point1);
    d2 = distancePoints (MLA_points(1,:), MLA_optimal_point2);
    if d1 > d2
      MLA_opt_point_A = MLA_optimal_point2;
      MLA_opt_point_B = MLA_optimal_point1;
    else
      MLA_opt_point_A = MLA_optimal_point1;
      MLA_opt_point_B = MLA_optimal_point2;
    endif
    ## print optimized MLA_points
    printf ("\nOptimized MLA points A and B are:\n");
    printf ("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", ...
            MLA_opt_point_A, MLA_opt_point_B);
  else
    if ! register
      ## print initial user defined MLA_points
      printf ("\nUser defined point for orienting Ulna is:\n");
      printf ("x: %f y: %f z: %f\n", MLA_points(1,:));
    else
      ## print initial user defined MLA_points
      printf ("\nAutomatically registered point for orienting Ulna is:\n");
      printf ("x: %f y: %f z: %f\n", MLA_points(1,:));
    endif
  endif
  ## save user defined and optimized MLA point back to the Meshlab Point file
  MLP = [MLA_points([1:2],:); MLA_opt_point_A; MLA_opt_point_B];
  write_MeshlabPoints (fullfile (folder, filenamePP), filename, MLP);
	
  ## calculate new normals for each final cross section at 20, 35, 50, 65 and 80%
  ## all normals should be pointing upwards, that is from distal towards proximal
  n1 = (Centroid_1 - Centroid_2) ./ sqrt (sum ((Centroid_1 - Centroid_2) .^2 ));
  n2 = (Centroid_1 - Centroid_3) ./ sqrt (sum ((Centroid_1 - Centroid_3) .^2 ));
  n3 = (Centroid_2 - Centroid_4) ./ sqrt (sum ((Centroid_2 - Centroid_4) .^2 ));
  n4 = (Centroid_3 - Centroid_5) ./ sqrt (sum ((Centroid_3 - Centroid_5) .^2 ));
  n5 = (Centroid_4 - Centroid_5) ./ sqrt (sum ((Centroid_4 - Centroid_5) .^2 ));
	
  ## calculate new cross sections for every point in centroid_N by using the
  ## normal of the nearest default centroid (20%, 35%, 50%, 65%, 80%)
  section_1 = meshSection (v, f, Centroid_1, n1);
  [CS_Geometry, SMoA, polyline] = simple_polygon3D (section_1, n1, ...
                                                    CorPlane_normal);
  ## arrange to matrices for csv files
  i = 1;
  geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter, ...
                   CS_Geometry.Centroid, n1, CorPlane_normal];
  inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                  SMoA.Imax, SMoA.theta];
  polygon2D(1,i*2-1) = Centroid_N(i,1);
  polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
  ## for user defined points
  CENTROIDS = [Centroid_1; Centroid_2; Centroid_3; Centroid_4; Centroid_5];
  NORMALS = [n1; n2; n3; n3; n4; n5];
  for i = 2: size (Centroid_N, 1) - 1
    ## find nearest default centroid
    dist2centroid = distancePoints (Centroid_N(i,[2:4]), CENTROIDS);
    index = find (dist2centroid == min (dist2centroid));
    ## calculate cross section with the normal of the nearest centroid
    section = meshSection (v, f, Centroid_N(i,[2:4]), NORMALS(index,:));
    [CS_Geometry, SMoA, polyline] = simple_polygon3D (section, ...
                                    NORMALS(index,:), CorPlane_normal);
    ## arrange to matrices for csv files
    geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter,...
                     CS_Geometry.Centroid, NORMALS(index,:), CorPlane_normal];
    inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                    SMoA.Imax, SMoA.theta];
    polygon2D(1,i*2-1) = Centroid_N(i,1);
    polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;  
  endfor
  section_5 = meshSection (v, f, Centroid_5, n5);
  [CS_Geometry, SMoA, polyline] = simple_polygon3D (section_5, n5, ...
                                                    CorPlane_normal);
  ## arrange to matrices for csv files
  i += 1;
  geometry(i,:) = [Centroid_N(i,1), CS_Geometry.Area, CS_Geometry.Perimeter, ...
                   CS_Geometry.Centroid, n5, CorPlane_normal];
  inertia(i,:) = [Centroid_N(i,1), SMoA.Ix, SMoA.Iy, SMoA.Ixy, SMoA.Imin, ...
                  SMoA.Imax, SMoA.theta];
  polygon2D(1,i*2-1) = Centroid_N(i,1);
  polygon2D([2:length(polyline.poly2D)+1],[i*2-1:i*2]) = polyline.poly2D;
  
  ## save to files
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
endfunction
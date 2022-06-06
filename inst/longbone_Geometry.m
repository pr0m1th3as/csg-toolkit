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
## @deftypefn{Function} [@var{CS_Geometry}, @var{SMoA}, @var{polyline}] = longbone_Geometry (@var{folder}, @var{filename}, @var{bones})
## @deftypefnx{Function} [@var{CS_Geometry}, @var{SMoA}, @var{polyline}] = longbone_CustomGeometry (@var{folder}, @var{filename})
## @deftypefnx{Function} [@var{CS_Geometry}, @var{SMoA}, @var{polyline}] = longbone_CustomGeometry (@var{filename}, @var{bones})
## @deftypefnx{Function} [@var{CS_Geometry}, @var{SMoA}, @var{polyline}] = longbone_CustomGeometry (@var{filename})
##
## This function slices the provided mesh of a humerus, ulna, femur or tibia at
## 20%, 35%, 50%, 65% and 80% along the bone's maximum length and returns the
## geometric properties for each planar cross section.  The mesh is assumed to
## be in mm.
##
## It loads a triangular mesh stored in Wavefront OBJ file and, if available,
## the respective Meshlab Point file that contains the two initial alignment
## points that define the mediolateral axis of the bone's anatomical orientation.
## Both files must be present in the same directory and share a common base name.
## If the respective Meshlab Point file is missing, the 'longbone_Registration'
## function is utilized to automatically register the initial alignment points.
## If no @var{folder} path is provided, then the working directory is assumed.
##
## The function requires at least 1 input argument, in which case it should be
## a string with the filename of the bone's OBJ mesh file ("bone.obj"). When 
## called with two input arguments there are two options:
##
## @itemize
## @item a string with the absolute path of the containing folder and another
## string with the model's filename.
## @item a string with the model's filename and another variable for bone
## selection, which can also be a string or a numerical vector.
## @end itemize
##
## Single bone selection can be done with the strings, "Humerus", 'Ulna",
## "Femur", and "Tibia".  Any other string value in @var{bones} is ignored and
## the bone is automatically determined by the 'longbone_Registration' function.
## When the function is called with three input arguments, @var{bones} must be
## a numerical vector with int values 1:5, which specify multiple bone selection
## as follows:
## @enumerate
## @item Humerus @item Ulna @item Femur @item Tibia @item All 
## @end enumerate
## The order of integers in @var{bones} is irrelevant.  If @var{folder} is an
## empty string, the current working directory is assumed.
##
## The mediolateral axis alignment points for the humerus need to be positioned
## anteriorly on the trochlea and capitulum at the distal end of the humerus.
## For the femur, the points need to be located posteriorly on the medial and
## lateral condyles at the distal end of the femur, whereas for the tibia, the
## points should be positioned posteriorly of the medial and lateral condyles
## at the proximal end of the tibia. The function will automatically optimize
## the local extremal points of the mediolateral axis based on the initial user
## defined points and will use the optimized coronal plane orientation for the
## calculations of the second moments of area. For ulna, only a single point of
## reference is required placed at the anterior tip of the coronoid process.
## Note, this reference point in not optimized.  If no Meshlab Point file is
## available, the reference point for ulna is registered at the anterior tip
## along the Imin axis as calculated from a 2D top view projection of the
## proximal epiphysis.
##
## When providing orientation points through the respective Meshlab point file,
## caution should be taken to that 'longbone_Geometry.m' function will only
## use the first two points available in the file.  Furthermore, the optimized
## points are appended in the original .pp file after the user defined points.
## So, if more points are present in the original .pp file, these are disregarded
## and eventually deleted from the final .pp file after processing.  In cases of
## automatically registering alignment points, these and their optimized
## counterparts are saved in a newly created pp file in the same folder with
## the .obj file following its name convention.  For ulna, only the first point
## available is considered.
##
## Independently of the number of input arguments, @code{longbone_Geometry}
## returns 3 structure arrays as output arguments with the following fields.
##
## @multitable @columnfractions .25 .75
## @item @code{CS_Geometry(n)} @tab is an array structure, where n=[1:5] defines
## the order of section from proximal 20% to distal 80% of the bone's max length
## and contains the following fields:
## @end multitable
##
## @multitable @columnfractions .05 .4 .55
## @item @tab @code{CS_Geometry(n).Area} @tab calculated in mm^2.
## @item @tab @code{CS_Geometry(n).Perimeter} @tab calculated in mm.
## @item @tab @code{CS_Geometry(n).Centroid} @tab x,y,z centroid coordinates of
## the cross section.
## @item @tab @code{CS_Geometry(n).Section_n} @tab x,y,z normal coordinates of
## the cross sectional plane.
## @item @tab @code{CS_Geometry(n).Coronal_n} @tab x,y,z normal coordinates of
## the coronal plane.
## @end multitable
##
## @multitable @columnfractions .25 .75
## @item @code{SMoA(n)} @tab is an array structure containing the fields:
## @end multitable
##
## @multitable @columnfractions .05 .25 .7
## @item @tab @code{SMoA(n).Ix} @tab 2nd moment of area with respect to x axis,
## which is collinear with the intersection of the coronal plane and the
## sectioning plane.  Calculated in mm^4.
## @item @tab @code{SMoA(n).Iy} @tab 2nd moment of area with respect to y axis,
## which is collinear with the intersection of the sagital plane and the
## sectioning plane.  Calculated in mm^4.
## @item @tab @code{SMoA(n).Ixy} @tab product 2nd moment of area in mm^4.
## @item @tab @code{SMoA(n).Imin} @tab minimum 2nd moment of area in mm^4.
## @item @tab @code{SMoA(n).Imax} @tab maximum 2nd moment of area in mm^4.
## @item @tab @code{SMoA(n).theta} @tab angle of rotation of the principal axis
## of 2nd moment of area with respect to x axis, expressed in degrees.
## @end multitable
##
## @multitable @columnfractions .25 .75
## @item @code{polyline(n)} @tab is an array structure containing the fields:
## @end multitable
##
## @multitable @columnfractions .05 .4 .55
## @item @tab @code{polyline(n).poly2D} @tab Nx2 matrix containing the (x,y)
## coordinates of the cross section on the 2D local axes of the sectioning plane
## ordered counter-clockwise.
## @item @tab @code{polyline(n).poly3D} @tab Nx3 matrix containing the original
## 3D coordinates of the cross section ordered counter clockwise.
## @end multitable
##
## @seealso{inspect_CSG, visualize_CrossSections, longbone_CustomGeometry}
## @end deftypefn

function [CS_Geometry, SMoA, polyline] = longbone_Geometry (varargin)
  ## declare empty output variables so that returning does not produce error
  CS_Geometry = [];
  SMoA = [];
  polyline = [];
  
  ## check input variables and parse them accordingly
  if nargin < 1 || nargin > 3
    printf ("Invalid number of input arguments\n");
    return;
  endif
  ## only 1 arg: filename
  if (nargin == 1 && ischar (varargin{1}(:)'))
    folder = "";
    filename = varargin{1}(:)';
    find_bone = true;
  elseif (nargin == 1 && !ischar (varargin{1}(:)'))
    printf("Filename must be a string\n");
    return;
  endif
  ## only 2 args: filename and bone or folder and filename
  if (nargin == 2 && ischar (varargin{1}(:)') && ischar (varargin{2}(:)'))
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
  ## only 2 args: filename and bone (numeric)
  if (nargin == 2 && ischar (varargin{1}(:)') && isnumeric (varargin{2}))
    folder = "";
    filename = varargin{1}(:)';
    bonesnum = varargin{2};
  endif
  ## only 3 args: folder, filename, and bone (numeric)
  if (nargin == 3 && ischar (varargin{1}(:)') && ischar (varargin{2}(:)') ...
      && isnumeric (varargin{3}))
    folder = varargin{1}(:)';
    filename = varargin{2}(:)';
    bonesnum = varargin{3};
  endif
  
  ## check filename has a valid .obj extension
  if ! (strcmpi (filename([end - 3:end]), ".obj"))
    printf ("Model must be in OBJ file format\n");
    return;
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
  if (nargin == 1)
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
  ## calculate the normal vector and the 5 points in R3 that define the 5
  ## slicing planes at 20, 35, 50, 65 and 80% of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt(sum((maxd_V2 - maxd_V1).^2));
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
	
  ## for each cross section calculate its centroid and area
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
  
  ## calculate the proximal distal axis of the bone and use it together with
  ## the mediolateral axis vector to define the coronal plane and calculate
  ## its normal. Additionally, check if the centroids' locations progress from
  ## proximal (20%) to distal (80%) and if required reverse their order to
  ## comply with this standard.
  proximal = distancePoints (Centroid_1, MLA_points(1,:));
  distal = distancePoints (Centroid_5, MLA_points(1,:));
  if (strcmp (bone, "Humerus") || strcmp (bone, "Femur")) && (proximal < distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
  elseif (strcmp (bone, "Tibia") || strcmp (bone, "Ulna")) && (proximal > distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
  endif

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
    ## calculate the coronal plane unit vector by normalizing the vector from
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
	
  ## calculate the sectioning points for each plane
  section_1 = meshSection (v, f, Centroid_1, n1);
  section_2 = meshSection (v, f, Centroid_2, n2);
  section_3 = meshSection (v, f, Centroid_3, n3);
  section_4 = meshSection (v, f, Centroid_4, n4);
  section_5 = meshSection (v, f, Centroid_5, n5);
	
  ## for each cross section calculate its centroid and area
  [CS_Geometry(1), SMoA(1), polyline(1)] = simple_polygon3D (section_1, n1, ...
                                                             CorPlane_normal);
  [CS_Geometry(2), SMoA(2), polyline(2)] = simple_polygon3D (section_2, n2, ...
                                                             CorPlane_normal);
  [CS_Geometry(3), SMoA(3), polyline(3)] = simple_polygon3D (section_3, n3, ...
                                                             CorPlane_normal);
  [CS_Geometry(4), SMoA(4), polyline(4)] = simple_polygon3D (section_4, n4, ...
                                                             CorPlane_normal);
  [CS_Geometry(5), SMoA(5), polyline(5)] = simple_polygon3D (section_5, n5, ...
                                                             CorPlane_normal);

  ## store the normals for each sectioning plane and the coronal plane normal
  CS_Geometry(1).Section_n = n1; CS_Geometry(1).Coronal_n = CorPlane_normal;
  CS_Geometry(2).Section_n = n2; CS_Geometry(2).Coronal_n = CorPlane_normal;
  CS_Geometry(3).Section_n = n3; CS_Geometry(3).Coronal_n = CorPlane_normal;
  CS_Geometry(4).Section_n = n4; CS_Geometry(4).Coronal_n = CorPlane_normal;
  CS_Geometry(5).Section_n = n5; CS_Geometry(5).Coronal_n = CorPlane_normal;
    
  ## print results for centroids and cross sectional areas
  printf ("\n%s in %s has a maximum distance of %f mm\n\n", ...
          bone, filename, maxDistance);
  printf ("Cross section at 20%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(1).Area, CS_Geometry(1).Perimeter, CS_Geometry(1).Centroid);
  printf ("Cross section at 35%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(2).Area, CS_Geometry(2).Perimeter, CS_Geometry(2).Centroid);
  printf ("Cross section at 50%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(3).Area, CS_Geometry(3).Perimeter, CS_Geometry(3).Centroid);
  printf ("Cross section at 65%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(4).Area, CS_Geometry(4).Perimeter, CS_Geometry(4).Centroid);
  printf ("Cross section at 80%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(5).Area, CS_Geometry(5).Perimeter, CS_Geometry(5).Centroid);
endfunction
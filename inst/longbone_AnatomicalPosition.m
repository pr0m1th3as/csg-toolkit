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
## @deftypefn  {csg-toolkit} {[@var{V2}, @var{F2}] =} longbone_AnatomicalPosition (@var{V1}, @var{F1})
## @deftypefnx {csg-toolkit} {[@var{V2}, @var{F2}, @var{VN2}, @var{FN2}] =} @
## longbone_AnatomicalPosition (@var{V1}, @var{F1}, @var{VN1}, @var{FN1})
##
## This function reorients the mesh of a longbone into its anatomical posistion.
##
## @seealso{longbone_Geometry, meshRotation, readObj}
## @end deftypefn

function [varargout] = longbone_AnatomicalPosition (varargin)

  ## Check input output
  if (! (nargin == 2 || nargin == 4))
    print_usage;
  endif
  if (nargin != nargout)
    print_usage;
  endif

  ## Get input arguments
  v = varargin{1};
  f = varargin{2};
  normals = false;
  if (nargin == 4)
    vn = varargin{3};
    fn = varargin{4};
    normals = true;
  endif

  ## Check input arguments
  if (size (v, 2) != 3)
    error
  endif
  if (size (f, 2) != 3)
    error
  endif
  if (normals)
    if (size (vn, 2) != 3)
      error
    endif
    if (size (fn, 2) != 3)
      error
    endif
  endif

  ## Translate barycenter to origin
  for i = 1:10
    v = v - meshBarycenter (v, f);
  endfor

  ## Compute anatomical plane normals
  [TransPlane_normal, ~] = anatomicalNormals (v, f);

  ## Rotate longtitudinal axis to [0,0,1]
  v_rot = meshRotation (v, TransPlane_normal, [0,0,1]);
  if (normals)
    vn_rot = meshRotation (vn, TransPlane_normal, [0,0,1]);
  endif

  ## Compute anatomical plane normals
  [~, CorPlane_normal] = anatomicalNormals (v_rot, f);

  ## Rotate anteroposterior axis to [1,0,0]
  v_rot = meshRotation (v_rot, CorPlane_normal, [1,0,0]);
  if (normals)
    vn_rot = meshRotation (vn_rot, CorPlane_normal, [1,0,0]);
  endif

  ## Return arguments
  varargout{1} = v_rot;
  varargout{2} = f;
  if (normals)
    varargout{3} = vn_rot;
    varargout{4} = fn;
  endif

endfunction

function [TransPlane_normal, CorPlane_normal] = anatomicalNormals (v, f)
  ## Identify longbone and register initial orientation
  [bone, MLA_points] = longbone_Registration (v, f);

  ## Check for valid longbone in 3D mesh
  if (! (any (strcmpi (bone, {"Humerus", "Ulna", "Femur", "Tibia"}))))
    error ("Bone should be either a Humerus, an Ulna, a Femur, or a Tibia.");
  endif

  ## Compute the maximum distance of the bone
  [~, maxd_V1, maxd_V2] = longbone_maxDistance (v);

  ## Calculate the normal vector of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt (sum ((maxd_V2 - maxd_V1) .^ 2));

  ## Calculate planes at 20 and 80% of the maximum length of the bone
  point_1 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.20;
  point_5 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.80;

  ## Calculate the sectioning points for each plane
  plane_1 = meshSection (v, f, point_1, normal);
  plane_5 = meshSection (v, f, point_5, normal);

  ## Calculate centroid and area for each cross section
  CS_Geometry(1) = simple_polygon3D (plane_1, normal);
  CS_Geometry(2) = simple_polygon3D (plane_5, normal);

  ## Store the centroids of the initial cross-sectional areas
  Centroid_1 = CS_Geometry(1).Centroid;
  Centroid_5 = CS_Geometry(2).Centroid;

  ## Check the centroids' locations from proximal (20%) to distal (80%) and
  ## if required reverse their order to comply with this standard.
  proximal = distancePoints (Centroid_1, MLA_points(1,:));
  distal = distancePoints (Centroid_5, MLA_points(1,:));
  if ((strcmpi (bone, "Humerus") || strcmpi (bone, "Femur")) &&
      (proximal < distal))
    Centroid_1 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
  elseif ((strcmpi (bone, "Tibia") || strcmpi (bone, "Ulna")) &&
          (proximal > distal))
    Centroid_1 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
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
  elseif (strcmpi (bone, "Tibia"))
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
endfunction

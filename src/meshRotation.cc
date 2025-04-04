/*
Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <octave/oct.h>

using namespace std;

struct Vector3            //for 3D points
{
  double x, y, z;
};

// function for normalizing a vector
Vector3 normalizeVector3D (Vector3 A)
{
  Vector3 normal;
  double v_len = sqrt (A.x * A.x + A.y * A.y + A.z * A.z);
  if (v_len != 0)
  {
    normal = {A.x / v_len, A.y / v_len, A.z / v_len};
  }
  else
  {
    normal = {0, 0, 0};
  }
  return normal;
}
// function for calculating the dot product between two 3D vectors
double dotProduct3D (Vector3 A, Vector3 B)
{
	double dotprod = A.x * B.x + A.y * B.y + A.z * B.z;
	return dotprod;
}
// function for calculating CrossProduct between two 3D vectors
Vector3 crossProduct3D (Vector3 A, Vector3 B)
{
  double x, y, z;
  x = A.y * B.z - A.z * B.y;
  y = A.z * B.x - A.x * B.z;
  z = A.x * B.y - A.y * B.x;
  Vector3 crossprod = {x, y, z};
  return crossprod;
}
// function for rotating a 3D point cloud about an arbitrary axis by an angle theta
vector<Vector3> rotatePoints3D (vector<Vector3> Cloud3D, Vector3 rot_axis, double theta)
{
  vector<Vector3> rCloud3D;
  double c = cos (theta);
  double c_ = 1 - c;
  double s = sin (theta);
  for(int i = 0; i < Cloud3D.size(); ++i)
  {
    Vector3 tmp1 = crossProduct3D (rot_axis, Cloud3D[i]);
    double tmp2 = dotProduct3D (rot_axis, Cloud3D[i]);
    double x = Cloud3D[i].x * c + tmp1.x * s + rot_axis.x * tmp2 * c_;
    double y = Cloud3D[i].y * c + tmp1.y * s + rot_axis.y * tmp2 * c_;
    double z = Cloud3D[i].z * c + tmp1.z * s + rot_axis.z * tmp2 * c_;
    rCloud3D.push_back ({x, y, z});
  }
  return rCloud3D;
}

DEFUN_DLD (meshRotation, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {csg-toolkit} {[@var{bone}, @var{Points}] =} meshRotation (@var{v}, @var{n1}, @var{n2})\n\
\n\
\n\
Rotate a 3D mesh. \
\n\
\n\
This function rotates an arbitrary number of points in R3 given in @var{v} \
according to the dihedral angle between a source normal @var{n1} and a target \
normal @var{n2}.  The pivot of rotation is the origin of the cartesian \
coordinates and the axis of rotation is the perpendicular to the source and \
target normals. \
\n\
\n\
@seealso{meshArea, meshBarycenter, meshSection} \n\
@end deftypefn")
{
  // check number of input and output arguments
  if (args.length() != 3)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  if (nargout < 1 && nargout > 2)
  {
    cout << "Invalid number of output arguments.\n";
    return octave_value_list();
  }
  // store vertices and normals into vectors
  Matrix V = args(0).array_value();
  Matrix n1 = args(1).array_value();
  Matrix n2 = args(2).array_value();
  octave_idx_type V_rows = args(0).rows();
  vector<Vector3> Cloud3D;
  double tmpx, tmpy, tmpz;
  for (octave_idx_type i = 0; i < V_rows; i++)
  {
    tmpx = V(i,0);
    tmpy = V(i,1);
    tmpz = V(i,2);
    Vector3 temp_points = {tmpx, tmpy, tmpz};
    Cloud3D.push_back(temp_points);
  }
  tmpx = n1(0,0); tmpy = n1(0,1); tmpz = n1(0,2);
  Vector3 source_normal = {tmpx, tmpy, tmpz};
  tmpx = n2(0,0); tmpy = n2(0,1); tmpz = n2(0,2);
  Vector3 target_normal = {tmpx, tmpy, tmpz};
  // find angle of rotation source and target normals
  Vector3 st_norm = crossProduct3D (source_normal, target_normal);
  Vector3 rotAxis = normalizeVector3D (st_norm);
  double theta_nom = dotProduct3D (st_norm, rotAxis);
  double theta_den = dotProduct3D (source_normal, target_normal);
  double theta = atan2 (theta_nom, theta_den);
  // rotate model so that 20% cross-sectional normal aligns with x axis
  // note: the nearest epiphysis lies on the positive side of x axis
  Cloud3D = rotatePoints3D (Cloud3D, rotAxis, theta);
  // return rotated 3D points
  for (octave_idx_type i = 0; i < V_rows; i++)
  {
    V(i,0) = Cloud3D[i].x;
    V(i,1) = Cloud3D[i].y;
    V(i,2) = Cloud3D[i].z;
  }
  // define return value list
  octave_value_list retval;
  retval(0) = V;
  return octave_value(retval);
}

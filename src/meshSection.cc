/*
Copyright (C) 2018-2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>

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
#include <string>
#include <fstream>
#include <vector>
#include <octave/oct.h>
#include <octave/parse.h>

using namespace std;

struct Vector3          //for point and normal coordinates
{
  double x, y, z;
};
struct Faces
{
  int a, b, c;
};
struct V2F             //for vertices of adjacent faces
{
  int V1, V2;
};
// function for calculating DotProduct between two vectors
double dotProduct3D(Vector3 A, Vector3 B)
{
  double dotprod = A.x * B.x + A.y * B.y + A.z * B.z;
  return dotprod;
}
// function for calculating distance between two 3D points
double distancePoints3D(Vector3 A, Vector3 B)
{
  double distance;
  double x, y, z;
  x = (A.x - B.x) * (A.x - B.x);
  y = (A.y - B.y) * (A.y - B.y);
  z = (A.z - B.z) * (A.z - B.z);
  distance = sqrt(x + y + z);
  return distance;
}
// function for extracting unique ordered intersection points given a plane
// specified by its normal and a point lying on that plane
vector<Vector3> sliceMesh(vector<Vector3> V3D, vector<Faces> Face3V, Vector3 point, Vector3 normal)
{
  // initialize vectors and variables
  vector<Vector3> iPoints, uPoints;
  double V1x, V1y, V1z, V2x, V2y, V2z, V3x, V3y, V3z;
  Vector3 V1, V2, V3;
  double dotV1, dotV2, dotV3;
  // for each triplet of vertices find faces intersected by specified plane
  for(int i = 0; i < Face3V.size(); ++i)
  {
    V1x = V3D[Face3V[i].a].x; V1y = V3D[Face3V[i].a].y; V1z = V3D[Face3V[i].a].z;
    V2x = V3D[Face3V[i].b].x; V2y = V3D[Face3V[i].b].y; V2z = V3D[Face3V[i].b].z;
    V3x = V3D[Face3V[i].c].x; V3y = V3D[Face3V[i].c].y; V3z = V3D[Face3V[i].c].z;
    V1 = {V1x - point.x, V1y - point.y, V1z - point.z};
    V2 = {V2x - point.x, V2y - point.y, V2z - point.z};
    V3 = {V3x - point.x, V3y - point.y, V3z - point.z};
    dotV1 = dotProduct3D(V1, normal);
    dotV2 = dotProduct3D(V2, normal);
    dotV3 = dotProduct3D(V3, normal);
    // true if intersecting
    if(!((dotV1 < 0 && dotV2 < 0 && dotV3 < 0) || 
         (dotV1 > 0 && dotV2 > 0 && dotV3 > 0)))
		{
			if(dotV1 * dotV2 < 0)       // vertices 1 and 2 lie on opposite sides
			{
				Vector3 PL = {point.x - V1x, point.y - V1y, point.z - V1z};
				Vector3 L = {V2x - V1x, V2y - V1y, V2z - V1z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv1v2 = {V1x + d * L.x, V1y + d * L.y, V1z + d * L.z};
				iPoints.push_back(CSv1v2);
			}
			if(dotV1 * dotV3 < 0)       // vertices 1 and 3 lie on opposite sides
			{
				Vector3 PL = {point.x - V1x, point.y - V1y, point.z - V1z};
				Vector3 L = {V3x - V1x, V3y - V1y, V3z - V1z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv1v3 = {V1x + d * L.x, V1y + d * L.y, V1z + d * L.z};
				iPoints.push_back(CSv1v3);
			}
			if(dotV2 * dotV3 < 0)       // vertices 2 and 3 lie on opposite sides
			{
				Vector3 PL = {point.x - V2x, point.y - V2y, point.z - V2z};
				Vector3 L = {V3x - V2x, V3y - V2y, V3z - V2z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv2v3 = {V2x + d * L.x, V2y + d * L.y, V2z + d * L.z};
				iPoints.push_back(CSv2v3);
			}
			if(dotV1 == 0)      // vertex 1 lies on the cross-section plane
			{
				Vector3 CSv1;
				CSv1 = {V1x, V1y, V1z};
				iPoints.push_back(CSv1);
			}
			if(dotV2 == 0)      // vertex 2 lies on the cross-section plane
			{
				Vector3 CSv2;
				CSv2 = {V2x, V2y, V2z};
				iPoints.push_back(CSv2);
			}
			if(dotV3 == 0)      // vertex 3 lies on the cross-section plane
			{
				Vector3 CSv3;
				CSv3 = {V3x, V3y, V3z};
				iPoints.push_back(CSv3);
			}
		}
  }
  // check if any faces were intersected
  if(iPoints.size() > 0)
  {
    // remove duplicate intersected points
    for(int i = 0; i < iPoints.size() - 1; ++i)
    {
      for(int j = i + 1; j < iPoints.size(); ++j)
      {
        double dx_ij = abs(iPoints[i].x - iPoints[j].x);
        double dy_ij = abs(iPoints[i].y - iPoints[j].y);
        double dz_ij = abs(iPoints[i].z - iPoints[j].z);
        if(dx_ij < 0.000000001 && dy_ij < 0.000000001 && dz_ij < 0.000000001)
        {
          iPoints.erase(iPoints.begin() + j);
          iPoints.push_back({iPoints[i].x, iPoints[i].y, iPoints[i].z});
        }
      }
    }
    bool duplicate_found = false;
    int i = 0;
    while(!duplicate_found || i < iPoints.size() - 1)
    {
      for(int j = i + 1; j < iPoints.size(); ++j)
      {
        double dx_ij = abs(iPoints[i].x - iPoints[j].x);
        double dy_ij = abs(iPoints[i].y - iPoints[j].y);
        double dz_ij = abs(iPoints[i].z - iPoints[j].z);
        if(dx_ij < 0.000000001 && dy_ij < 0.000000001 && dz_ij < 0.000000001)
        {
          iPoints.erase(iPoints.begin() + j, iPoints.end());
          duplicate_found = true;
        }
      }
      i++;
    }
    uPoints = iPoints;
  }
  else
  {
    cout << "No faces were intersected.\n";
    uPoints.push_back({NAN, NAN, NAN});
  }
  return uPoints;
}


DEFUN_DLD(meshSection, args, nargout, 
          "-*- texinfo -*-\n\
@deftypefn {Function} @var{cross_section} = meshSection (@var{v}, @var{f}, @var{point}, @var{normal})\n\
\n\
\n\
This function loads the vertices @var{v} and faces @var{f} of a triangular 3D\
Mesh along with a slicing plane defined by its @var{normal} and a @var{point}\
that lies on the slicing plane and it returns the intersection points of the\
face edges between the vertices that lie on opposite sides of the sectioning\
plane.  Duplicate points due to adjacent faces are removed and the points are\
sorted according to the οrder of the faces.\n\n\
@var{v} and @var{f} must be (Nx3) matrices containing 3D coordinates and vertex\
indices respectively.  @var{normal} and @var{point} must be defined as row\
vectors containing x, y, z coordinates in R3.  @var{cross_section} is the\
return variable Nx3 in size, where N is the number of unique intersection\
points of the 3D mesh represented by 'v' and 'f' input arguments.\n\
@end deftypefn")
{
  // count the number of input arguments and store their values
  // check for invalid number of input arguments
  if (args.length() != 4)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  // check vertices containing 3D coordinates
  if (!args(0).is_matrix_type() || args(0).columns() != 3)
  {
    cout << "Vertex matrix should be Nx3 containing x,y,z coordinates.\n";
    return octave_value_list();
  }
  // check mesh for being triangular
  if (!args(1).is_matrix_type() || args(1).columns() != 3)
  {
    cout << "Mesh should be triangular.\n";
    return octave_value_list();
  }
  // check point and normal being 1x3 arrays
  if (args(2).rows() != 1 || args(2).columns() != 3)
  {
    cout << "point should be 1x3 array containing x y z coordinates.\n";
    return octave_value_list();
  }
  if (args(3).rows() != 1 || args(3).columns() != 3)
  {
    cout << "normal should be 1x3 array containing x y z coordinates.\n";
    return octave_value_list();
  }
  // store vertices and faces in Mesh data structure
  Matrix V = args(0).array_value();
  Matrix F = args(1).array_value();
  vector<Vector3> V3D;
  vector<Faces> Face3V;
  // find number of vertices and faces
  octave_idx_type V_rows = args(0).rows();
  octave_idx_type F_rows = args(1).rows();
  // for each face store its corresponding vertex coordinates
  double tmpx, tmpy, tmpz;
  int tmpa, tmpb, tmpc;
  for (octave_idx_type i = 0; i < F_rows; i++)
  {
    tmpa = F(i,0) - 1;
    tmpb = F(i,1) - 1;
    tmpc = F(i,2) - 1;
    Faces temp_face = {tmpa, tmpb, tmpc};
    Face3V.push_back(temp_face);
  }
  for (octave_idx_type i = 0; i < V_rows; i++)
  {
    tmpx = V(i, 0);
    tmpy = V(i, 1);
    tmpz = V(i, 2);
    Vector3 temp_V3D = {tmpx, tmpy, tmpz};
    V3D.push_back(temp_V3D);
  }
  // store point and normal from input arguments
  Matrix P = args(2).array_value();
  Matrix N = args(3).array_value();
  Vector3 point = {P(0), P(1), P(2)};
  Vector3 normal = {N(0), N(1), N(2)};
  // calculate cross-sectioning points
  vector<Vector3> section_points = sliceMesh(V3D, Face3V, point, normal);
  // cast the cross-sectioning points into an octave type Matrix
  Matrix cross_section (section_points.size(), 3);
  for (octave_idx_type i = 0; i < section_points.size(); i++)
  {
    cross_section(i,0) = section_points[i].x;
    cross_section(i,1) = section_points[i].y;
    cross_section(i,2) = section_points[i].z;
  }
  // define return value list and return the cross-sectioning points
  octave_value_list retval;
  retval(0) = cross_section;
  return retval;
}
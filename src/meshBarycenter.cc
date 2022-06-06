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

struct VCoord
{
	    double x, y, z;
};
struct Area
{
      double a;
};  


DEFUN_DLD (meshBarycenter, args, nargout, 
          "-*- texinfo -*-\n\
@deftypefn{Function} @var{centroid} = meshBarycenter (@var{v},@var{f})\n\
\n\
\n\
This function computes the 3D coordinates of the polygon barycenter of a\
triangular 3D Mesh based on its vertices and faces provided as input arguments.\
\n\
The function will only take two input arguments. The first argument should\
be an Nx3 matrix containing the 3-dimensional coordinates of each vertex and\
the second argument should be an Nx3 matrix with each row containing the three\
vertices that form each face of the triangular mesh. The face matrix should\
contain explicitly non-zero integers referring to the existing vertices present\
in the first input argument.\
@end deftypefn")
{
  // count the number of input arguments and store their values
  // into the appropriate variables
  // check for invalid number of input arguments
  if (args.length() != 2)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  // check for both arguments being real matrices
  if (!args(0).is_matrix_type() || !args(1).is_matrix_type())
  {
    cout << "Both input arguments should be real matrices.\n";
    return octave_value_list();
  }
  // store vertices and faces in vectors
  Matrix V = args(0).array_value();
  Matrix F = args(1).array_value();
  // find number of vertices and faces
  octave_idx_type V_rows = args(0).rows();
  octave_idx_type F_rows = args(1).rows();
  octave_idx_type V_columns = args(0).columns();
  octave_idx_type F_columns = args(1).columns();
  // ensure that there are at least 3 vertices and one face in the mesh and
  // vertex and face matrices are Nx3 in size
  if (V_rows < 3)
  {
    cout << "There should be at least 3 vertices in the mesh.\n";
    return octave_value_list();
  }
  if (V_columns != 3)
  {
    cout << "Vertex matrix should be Nx3 containing x,y,z coordinates.\n";
    return octave_value_list();
  }
  if (F_rows < 1)
  {
    cout << "There should be at least 1 face in the mesh.\n";
    return octave_value_list();
  }
  if (F_columns != 3)
  {
    cout << "Face matrix should be Nx3 containing three vertices.\n";
    return octave_value_list();
  }
  // define matrices for storing vertex coordinates
  vector<VCoord> vertex_A;
  vector<VCoord> vertex_B;
  vector<VCoord> vertex_C;
  // define matrix for storing barycenter coordinates for each triangular face
  vector<VCoord> face_barycenter;
  // define a vector for storing face area
  vector<Area> face_area;
  //
  //
  // loop through every face of the mesh and create three vectors with the
  // corresponding vertex coordinates of each face
  for (octave_idx_type i = 0; i < F_rows; i++)
  {
    // define variables for the vertices of each face
    int vert_idx_A, vert_idx_B, vert_idx_C;
    // define variables for the coordinates of each vertex
    double tmpx, tmpy, tmpz;
    // store the vertex indexes
    vert_idx_A = F(i,0) - 1;
    vert_idx_B = F(i,1) - 1;
    vert_idx_C = F(i,2) - 1;
    // store the coordinates of the first vertex
    tmpx = V(vert_idx_A,0);
    tmpy = V(vert_idx_A,1);
    tmpz = V(vert_idx_A,2);
    VCoord temp_vertex_A = {tmpx, tmpy, tmpz};
    vertex_A.push_back(temp_vertex_A);
    // store the coordinates of the second vertex
    tmpx = V(vert_idx_B,0);
    tmpy = V(vert_idx_B,1);
    tmpz = V(vert_idx_B,2);
    VCoord temp_vertex_B = {tmpx, tmpy, tmpz};
    vertex_B.push_back(temp_vertex_B);
    // store the coordinates of the third vertex
    tmpx = V(vert_idx_C,0);
    tmpy = V(vert_idx_C,1);
    tmpz = V(vert_idx_C,2);
    VCoord temp_vertex_C = {tmpx, tmpy, tmpz};
    vertex_C.push_back(temp_vertex_C);
  }
  // loop through every face of the mesh, calculate barycenter coordinates
  // of each triangle and its area, and calculate total area of the mesh
  double total_area = 0;
  for(octave_idx_type i = 0; i < F_rows; i++)
  {
    // define variables for the coordinates of each vertex
    double tmpx, tmpy, tmpz;
    // compute barycenter coordinates, i.e. centroid of face triangle
    tmpx = (vertex_A[i].x + vertex_B[i].x + vertex_C[i].x) / 3;
    tmpy = (vertex_A[i].y + vertex_B[i].y + vertex_C[i].y) / 3;
    tmpz = (vertex_A[i].z + vertex_B[i].z + vertex_C[i].z) / 3;
    VCoord centroid = {tmpx, tmpy, tmpz};
    face_barycenter.push_back(centroid);
    // compute area of triangular face
    // calculate vectors AB and AC
    VCoord vec_AB = {vertex_B[i].x - vertex_A[i].x, vertex_B[i].y
                      - vertex_A[i].y, vertex_B[i].z - vertex_A[i].z};
    VCoord vec_AC = {vertex_C[i].x - vertex_A[i].x, vertex_C[i].y
                      - vertex_A[i].y, vertex_C[i].z - vertex_A[i].z};
    // calculate ABxAC cross product
    tmpx = vec_AB.y * vec_AC.z - vec_AB.z * vec_AC.y;
    tmpy = vec_AB.z * vec_AC.x - vec_AB.x * vec_AC.z;
    tmpz = vec_AB.x * vec_AC.y - vec_AB.y * vec_AC.x;
    double temp_a = 0.5 * sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);
    Area temp_area = {temp_a};
    face_area.push_back(temp_area);
    total_area += temp_a;
  }
  // calculate the barycenter of the mesh as the weighted average of the
  // barycenter coordinates (centroids) of the face triangles weighted by
  // their corresponding area over the total area of the mesh
  Matrix mesh_barycenter (1, 3);
  for(octave_idx_type i = 0; i < F_rows; i++)
  {
    double weight = face_area[i].a / total_area;
    mesh_barycenter(0,0) += face_barycenter[i].x * weight;
    mesh_barycenter(0,1) += face_barycenter[i].y * weight;
    mesh_barycenter(0,2) += face_barycenter[i].z * weight;
  }
  // divide by the total number of faces
  mesh_barycenter(0,0) = (mesh_barycenter(0,0) * total_area) / F_rows;
  mesh_barycenter(0,1) = (mesh_barycenter(0,1) * total_area) / F_rows;
  mesh_barycenter(0,2) = (mesh_barycenter(0,2) * total_area) / F_rows;
  // define return value list
  octave_value_list retval;
  if (nargout == 1)
  {
      retval(0) = mesh_barycenter;
  }
  else
  {
      cout << "Mesh barycentric coordinates are: x=" << mesh_barycenter(0,0) 
           << "  y=" << mesh_barycenter(0,1) << "  z=" << mesh_barycenter(0,2) 
           << "\n";
      return octave_value_list();
  }
  return retval;
}
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


DEFUN_DLD (meshArea, args, nargout, 
          "-*- texinfo -*-\n\
@deftypefn{Function} @var{area} = meshArea (@var{v},@var{f})\n\
\n\
\n\
This function computes the surface area of a triangular mesh based on its\
vertices and faces provided as input arguments.\n\
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
  // loop through every face of the mesh and calculate its area
  double total_area = 0;
  for(octave_idx_type i = 0; i < F_rows; i++)
  {
    // calculate vectors AB and AC
    VCoord vec_AB = {vertex_B[i].x - vertex_A[i].x, vertex_B[i].y
                      - vertex_A[i].y, vertex_B[i].z - vertex_A[i].z};
    VCoord vec_AC = {vertex_C[i].x - vertex_A[i].x, vertex_C[i].y
                      - vertex_A[i].y, vertex_C[i].z - vertex_A[i].z};
    // calculate ABxAC cross product
    double tmpx, tmpy, tmpz;
    tmpx = vec_AB.y * vec_AC.z - vec_AB.z * vec_AC.y;
    tmpy = vec_AB.z * vec_AC.x - vec_AB.x * vec_AC.z;
    tmpz = vec_AB.x * vec_AC.y - vec_AB.y * vec_AC.x;
    double temp_a = 0.5 * sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);
    total_area += temp_a;
  }
  // define return value list
  octave_value_list retval;
  if (nargout == 1)
  {
      retval(0) = total_area;
  }
  else
  {
      cout << "Mesh total area is " << total_area << "\n"; 
      return octave_value_list();
  }
  return retval; 
}
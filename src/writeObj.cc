/*
Copyright (C) 2018-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>

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
struct Faces
{
      int a, b, c;
};
struct TCoord
{
      double u, v  ;
};

// Function for asking user to overwrite existing file
string checkFilename (string filename)
{
  bool filename_exists = ifstream(filename.c_str()).good();
  if (filename_exists)
  {
    cout << "Filename already exists.\n";
    cout << "Do you want to replace? (yes or no)\n";
    string yes_or_no;
    getline(cin, yes_or_no);
    while (yes_or_no.compare("yes") && yes_or_no.compare("no"))
    {
      cout << "Please answer yes or no! ";
      getline(cin, yes_or_no);
    }
    if (yes_or_no.compare("yes"))
    {
      string newfilename;
      cout << "Please enter new filename: ";
      getline(cin, newfilename);
      // check for file extension and append .obj if not present
      string extension = ".obj";
      if (newfilename.compare(newfilename.length() - extension.length(),
          extension.length(), extension) != 0)
      {
        newfilename.append(extension);
      }
      filename = newfilename.c_str();
    }
  }
  return filename;
}


DEFUN_DLD (writeObj, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {csg-toolkit} {} writeObj (@var{v}, @var{f}, @var{filename})\n\
 @deftypefnx {csg-toolkit} {} writeObj (@var{v}, @var{f}, @var{vt}, @var{ft}, @var{filename})\n\
 @deftypefnx {csg-toolkit} {} writeObj (@var{v}, @var{f}, @var{vn}, @var{fn}, @var{filename})\n\
 @deftypefnx {csg-toolkit} {} writeObj (@var{v}, @var{f}, @var{vt}, @var{ft}, @var{vn}, @var{fn}, @var{filename})\n\
 @deftypefnx {csg-toolkit} {} writeObj (@dots{}, @qcode{\"info\"})\n\
 @deftypefnx {csg-toolkit} {} writeObj (@dots{}, @qcode{\"overwrite\"})\n\
 @deftypefnx {csg-toolkit} {@var{filename} =} writeObj (@dots{})\n\
\n\
\n\
This function saves a triangular 3D Mesh in a Wavefront OBJ file format \
according to its elements provided as input arguments.\n\
\n\
The function typically takes 3, 5 or 7 input arguments.  The last argument \
should be a char string referring to the filename of the OBJ file, whereas the \
first 2, 4, or 6 input arguments should be 2-dimensional matrices containing \
the elements of the triangular mesh in the following order:\n\
\n\
@var{v} can be either an @math{Nx3} matrix with @var{x},@var{y},@var{z} vertex \
coordinates or an @math{Nx6} matrix containing the vertex coordinates and the \
@var{r},@var{g},@var{b} color information for each vertex in the range \
@math{[0,1]}.\n\
\n\
@var{f} must be an @math{Mx3} matrix with the indices of each corresponding \
vertex, where @math{M} is the number of faces.\n\
\n\
@var{vt} must be an @math{Nx2} matrix, where @math{N} is the number of texture \
coordinates.\n\
\n\
@var{ft} must be an @math{Mx3} matrix with the indices of each corresponding \
texture coordinates.\n\
\n\
@var{vn} must be an @math{Nx3} matrix with @var{x},@var{y},@var{z} normal \
coordinates, where @math{N} is the number of normals.\n\
\n\
@var{fn} must be an @math{Nx3} matrix with the indices of each corresponding \
face normal.\n\
\n\
When 5 input arguments are provided, the function will determine whether there \
is a texture coordinates matrix or a vertex normals matrix by the size of the \
second dimension of the matrix provided as the third input argument. \
\n\
The function can also take one or two additional optional arguments added \
after the aforementioned syntaxes to control its behavior.  \
Use @qcode{\"info\"} as an optional input argument to display information \
about the saved 3D mesh.  Use @qcode{\"overwrite\"} to overwrite existing file \
without asking back the user for confirmation.\n\
\n\
The function can also return an output argument with the filename used for \
saving the mesh.  This can be useful when prompting the user for a new \
filename when file already exists.\n\
\n\
\n\
@seealso{readObj, renameObj} \n\
@end deftypefn")
{
  int nargin = args.length();
  // Add default options
  bool info = false;
  bool over = false;
  // Parse optional arguments "info" and "overwrite"
  while (args(nargin - 1).is_string() && args(nargin - 2).is_string())
  {
    if (args(nargin - 1).string_value() == "info")
    {
      info = true;
    }
    if (args(nargin - 1).string_value() == "overwrite")
    {
      over = true;
    }
    nargin--;
  }
  // count the number of input arguments and store their values
  // into the appropriate variables
  // check for invalid number of input arguments
  if (nargin != 3 && nargin != 5 && nargin != 7)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  ////
  ////
  ////
  //// for three input arguments
  ////
  ////
  ////
  // check for last argument being a string
  if (nargin == 3 && ! args(2).is_string())
  {
    cout << "Third input argument should be a string.\n";
    return octave_value_list();
  }
  // check for first two arguments being real matrices
  if (! args(0).is_matrix_type() || ! args(1).is_matrix_type())
  {
    cout << "The first two arguments should be real matrices.\n";
    return octave_value_list();
  }
  // considering two input arguments for vertices and faces respcectively and
  // and a third argument as string for filename under which the vertices and
  // will be saved
  if (nargin == 3 && args(2).is_string())
  {
    // store vertices, faces and filename
    Matrix V = args(0).array_value();
    Matrix F = args(1).array_value();
    string filename = args(2).string_value();
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
    if (V_columns != 3 && V_columns != 6)
    {
      cout << "Vertex matrix should be Nx3 containing x,y,z coordinates\n";
      cout << "or Nx6 containing vertex coordinates and r,g,b color information.\n";
      return octave_value_list();
    }
    if (F_rows < 1)
    {
      cout << "There should be at least 1 face in the mesh.\n";
      return octave_value_list();
    }
    if (F_columns != 3)
    {
      cout << "Face matrix should be Nx3 indexing three vertices.\n";
      return octave_value_list();
    }
    // check if filename exists
    if (! over)
    {
      filename = checkFilename (filename);
    }
    // open file
    ofstream outputFile(filename.c_str());
    if (! outputFile.is_open())
    {
      cout << "Error opening " << filename.c_str() << "for write\n";
      return octave_value_list();
    }
    else
    {
      // writing header to file
      outputFile << "#\n# OBJ File generated by GNU Octave\n# using 'writeObj' function\n";
      outputFile << "#\n# Object " << filename.c_str() << "\n#\n";
      outputFile << "# Vertices: " << V_rows << "\n";
      outputFile << "# Faces: " << F_rows << "\n#\n#\n\n";
      if (info) { cout << "Writing to file... "; }
      // write vertices to file
      if (V_columns == 3) // coordinates only
      {
        float tmpx, tmpy, tmpz;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << "\n";
        }
      }
      else                // coordinates and color
      {
        float tmpx, tmpy, tmpz, tmpR, tmpG, tmpB;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          tmpR = V(i,3);
          tmpG = V(i,4);
          tmpB = V(i,5);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << " ";
          outputFile << tmpR << " " << tmpG << " " <<  tmpB << "\n";
        }
      }
      // store faces to file
      int tmpa, tmpb, tmpc;
      for (octave_idx_type i = 0; i < F_rows; i++)
      {
        tmpa = F(i,0);
        tmpb = F(i,1);
        tmpc = F(i,2);
        outputFile << "f " << tmpa << " " << tmpb << " " << tmpc << "\n";
      }
      // close file
      outputFile.close();
      if (info) { cout << "done!\n"; }
    }
    if (info)
    {
      cout << "Mesh filename is " << filename.c_str() << "\n";
      cout << "Mesh has " << V_rows << " vertices.\n";
      cout << "Mesh has " << F_rows << " faces.\n";
    }
    // if output argument is requested return the filename
    if (nargout == 1)
    {
      return octave_value_list(filename.c_str());
    }
    return octave_value_list();
  }
  ////
  ////
  ////
  //// for five input arguments
  ////
  ////
  ////
  // check for last argument being a string
  if (nargin == 5 && ! args(4).is_string())
  {
    cout << "Fifth input argument should be a string.\n";
    return octave_value_list();
  }
  // check for first four arguments being real matrices
  if (! args(0).is_matrix_type() || ! args(1).is_matrix_type() ||
      ! args(2).is_matrix_type() || ! args(3).is_matrix_type())
  {
    cout << "The first four arguments should be real matrices.\n";
    return octave_value_list();
  }
  // considering four input arguments being real matrices, the first two are for
  // vertices and faces respcectively, whereas the latter will either be for
  // texture coordinates and texture faces or for vertex normals and face normals.
  // This will depend on the number of coordinates present in third input argument.
  // The last input argument should again be a string
  //
  // check for 5 input arguments with the third argument having two coordinates
  // present, namely u, v, and the last argument being a string
  if (nargin == 5 && args(2).columns() == 2 && args(4).is_string())
  {
    // store vertices, faces and filename
    Matrix V = args(0).array_value();
    Matrix F = args(1).array_value();
    Matrix VT = args(2).array_value();
    Matrix FT = args(3).array_value();
    string filename = args(4).string_value();
    // find number of vertices and faces
    octave_idx_type V_rows = args(0).rows();
    octave_idx_type F_rows = args(1).rows();
    octave_idx_type VT_rows = args(2).rows();
    octave_idx_type FT_rows = args(3).rows();
    octave_idx_type V_columns = args(0).columns();
    octave_idx_type F_columns = args(1).columns();
    octave_idx_type VT_columns = args(2).columns();
    octave_idx_type FT_columns = args(3).columns();
    // ensure that there are at least 3 vertices and one face in the mesh and
    // vertex and face matrices are Nx3.
    if (V_rows < 3)
    {
      cout << "There should be at least 3 vertices in the mesh.\n";
      return octave_value_list();
    }
    if (V_columns != 3 && V_columns != 6)
    {
      cout << "Vertex matrix should be Nx3 containing x,y,z coordinates\n";
      cout << "or Nx6 containing vertex coordinates and r,g,b color information.\n";
      return octave_value_list();
    }
    if (F_rows < 1)
    {
      cout << "There should be at least 1 face in the mesh.\n";
      return octave_value_list();
    }
    if (F_columns != 3)
    {
      cout << "Face matrix should be Nx3 indexing three vertices.\n";
      return octave_value_list();
    }
    // check that texture faces are an Nx3 matrix
    if (FT_columns != 3)
    {
      cout << "Texture face matrix should be Nx3 indexing three vertices.\n";
      return octave_value_list();
    }
    // check that matrices for faces and texture faces have the same length
    // to ensure that the iterator will not break when writing to the file
    if (F_rows != FT_rows)
    {
      cout << "Faces and texture faces should be the same length.\n";
      return octave_value_list();
    }
    // check if filename exists
    if (! over)
    {
      filename = checkFilename (filename);
    }
    // open file
    ofstream outputFile(filename.c_str());
    if (! outputFile.is_open())
    {
      cout << "Error opening " << filename.c_str() << "for write\n";
      return octave_value_list();
    }
    else
    {
      // writing header to file
      outputFile << "#\n# OBJ File generated by GNU Octave\n# using 'writeObj' function\n";
      outputFile << "#\n# Object " << filename.c_str() << "\n#\n";
      outputFile << "# Vertices: " << V_rows << "\n";
      outputFile << "# Faces: " << F_rows << "\n#\n#\n";
      int i = filename.length() - 3;
      string mtlfilename = filename.c_str();
      outputFile << "mtllib ./" << mtlfilename.replace(i,3, "mtl") << "\n\n";
      if (info) { cout << "Writing to file... "; }
      // write vertices to file
      if (V_columns == 3) // coordinates only
      {
        float tmpx, tmpy, tmpz;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << "\n";
        }
      }
      else                // coordinates and color
      {
        float tmpx, tmpy, tmpz, tmpR, tmpG, tmpB;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          tmpR = V(i,3);
          tmpG = V(i,4);
          tmpB = V(i,5);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << " ";
          outputFile << tmpR << " " << tmpG << " " <<  tmpB << "\n";
        }
      }
      // write texture coordinates to file
      float tmpu, tmpv;
      for (octave_idx_type i = 0; i < VT_rows; i++)
      {
        tmpu = VT(i,0);
        tmpv = VT(i,1);
        outputFile << "vt " << tmpu << " " << tmpv << "\n";
      }
      // write faces and texture faces to file
      int tmpa, tmpb, tmpc, tmpd, tmpe, tmpf;
      for (octave_idx_type i = 0; i < F_rows; i++)
      {
        tmpa = F(i,0);
        tmpb = F(i,1);
        tmpc = F(i,2);
        tmpd = FT(i,0);
        tmpe = FT(i,1);
        tmpf = FT(i,2);
        outputFile << "f " << tmpa << "/" << tmpd << " " << tmpb << "/"
                   << tmpe << " " << tmpc << "/" << tmpf << "\n";
      }
      // close file
      outputFile.close();
      if (info) { cout << "done!\n"; }
    }
    if (info)
    {
      cout << "Mesh filename is " << filename.c_str() << "\n";
      cout << "Mesh has " << V_rows << " vertices.\n";
      cout << "Mesh has " << F_rows << " faces.\n";
    }
    // if output argument is requested return the filename
    if (nargout == 1)
    {
      return octave_value_list(filename.c_str());
    }
    return octave_value_list();
  }
  //
  // check for 5 input arguments with the third argument having three coordinates
  // present, namely x, y, z, and the last argument being a string
  if (nargin == 5 && args(2).columns() == 3 && args(4).is_string())
  {
    // store vertices, faces and filename
    Matrix V = args(0).array_value();
    Matrix F = args(1).array_value();
    Matrix VN = args(2).array_value();
    Matrix FN = args(3).array_value();
    string filename = args(4).string_value();
    // find number of vertices and faces
    octave_idx_type V_rows = args(0).rows();
    octave_idx_type F_rows = args(1).rows();
    octave_idx_type VN_rows = args(2).rows();
    octave_idx_type FN_rows = args(3).rows();
    octave_idx_type V_columns = args(0).columns();
    octave_idx_type F_columns = args(1).columns();
    octave_idx_type VN_columns = args(2).columns();
    octave_idx_type FN_columns = args(3).columns();
    // ensure that there are at least 3 vertices and one face in the mesh and
    // vertex normals and face normals matrices are Nx3.
    if (V_rows < 3)
    {
      cout << "There should be at least 3 vertices in the mesh.\n";
      return octave_value_list();
    }
    if (V_columns != 3 && V_columns != 6)
    {
      cout << "Vertex matrix should be Nx3 containing x,y,z coordinates\n";
      cout << "or Nx6 containing vertex coordinates and r,g,b color information.\n";
      return octave_value_list();
    }
    if (F_rows < 1)
    {
      cout << "There should be at least 1 face in the mesh.\n";
      return octave_value_list();
    }
    if (F_columns != 3)
    {
      cout << "Face matrix should be Nx3 indexing three vertices.\n";
      return octave_value_list();
    }
    // check that face normals are a Nx3 matrix
    if (FN_columns != 3)
    {
      cout << "Face normals should be an Nx3 matrix.\n";
      return octave_value_list();
    }
    // check that matrices for faces and face normals have the same length
    // to ensure that the iterator will not break when writing to the file
    if (F_rows != FN_rows)
    {
      cout << "Faces and face normals should be the same length.\n";
      return octave_value_list();
    }
    // check if filename exists
    if (! over)
    {
      filename = checkFilename (filename);
    }
    // open file
    ofstream outputFile(filename.c_str());
    if (! outputFile.is_open())
    {
      cout << "Error opening " << filename.c_str() << "for write\n";
      return octave_value_list();
    }
    else
    {
      // writing header to file
      outputFile << "#\n# OBJ File generated by GNU Octave\n# using 'writeObj' function\n";
      outputFile << "#\n# Object " << filename.c_str() << "\n#\n";
      outputFile << "# Vertices: " << V_rows << "\n";
      outputFile << "# Faces: " << F_rows << "\n#\n#\n\n";
      if (info) { cout << "Writing to file... "; }
      // write vertices to file
      if (V_columns == 3) // coordinates only
      {
        float tmpx, tmpy, tmpz;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << "\n";
        }
      }
      else                // coordinates and color
      {
        float tmpx, tmpy, tmpz, tmpR, tmpG, tmpB;
        for (octave_idx_type i = 0; i < V_rows; i++)
        {
          tmpx = V(i,0);
          tmpy = V(i,1);
          tmpz = V(i,2);
          tmpR = V(i,3);
          tmpG = V(i,4);
          tmpB = V(i,5);
          outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << " ";
          outputFile << tmpR << " " << tmpG << " " <<  tmpB << "\n";
        }
      }
      // write vertex normals to file
      float tmpx, tmpy, tmpz;
      for (octave_idx_type i = 0; i < VN_rows; i++)
      {
        tmpx = VN(i,0);
        tmpy = VN(i,1);
        tmpz = VN(i,2);
        outputFile << "vn " << tmpx << " " << tmpy << " " << tmpz << "\n";
      }
      // write faces and face normals to file
      int tmpa, tmpb, tmpc, tmpd, tmpe, tmpf;
      for (octave_idx_type i = 0; i < F_rows; i++)
      {
        tmpa = F(i,0);
        tmpb = F(i,1);
        tmpc = F(i,2);
        tmpd = FN(i,0);
        tmpe = FN(i,1);
        tmpf = FN(i,2);
        outputFile << "f " << tmpa << "//" << tmpd << " " << tmpb << "//"
                   << tmpe << " " << tmpc << "//" << tmpf << "\n";
      }
      // close file
      outputFile.close();
      if (info) { cout << "done!\n"; }
    }
    if (info)
    {
      cout << "Mesh filename is " << filename.c_str() << "\n";
      cout << "Mesh has " << V_rows << " vertices.\n";
      cout << "Mesh has " << F_rows << " faces.\n";
    }
    // if output argument is requested return the filename
    if (nargout == 1)
    {
      return octave_value_list(filename.c_str());
    }
    return octave_value_list();
  }
  ////
  ////
  ////
  //// for seven input arguments
  ////
  ////
  ////
  // check for last argument being a string
  if (nargin == 7 && ! args(6).is_string())
  {
    cout << "Seventh input argument should be a string.\n";
    return octave_value_list();
  }
  // check for first six arguments being real matrices
  if (! args(0).is_matrix_type() || ! args(1).is_matrix_type() ||
      ! args(2).is_matrix_type() || ! args(3).is_matrix_type() ||
      ! args(4).is_matrix_type() || ! args(5).is_matrix_type())
  {
    cout << "The first six arguments should be real matrices.\n";
    return octave_value_list();
  }
  // considering six input arguments being real matrices, the first two are for
  // vertices and faces respcectively, the next two are texture vertices and
  // texture faces, and the last two are vertex normals and face normals.
  //
  // store vertices, faces and filename
  Matrix V = args(0).array_value();
  Matrix F = args(1).array_value();
  Matrix VT = args(2).array_value();
  Matrix FT = args(3).array_value();
  Matrix VN = args(4).array_value();
  Matrix FN = args(5).array_value();
  string filename = args(6).string_value();
  // find number of vertices and faces
  octave_idx_type V_rows = args(0).rows();
  octave_idx_type F_rows = args(1).rows();
  octave_idx_type VT_rows = args(2).rows();
  octave_idx_type FT_rows = args(3).rows();
  octave_idx_type VN_rows = args(4).rows();
  octave_idx_type FN_rows = args(5).rows();
  octave_idx_type V_columns = args(0).columns();
  octave_idx_type F_columns = args(1).columns();
  octave_idx_type VT_columns = args(2).columns();
  octave_idx_type FT_columns = args(3).columns();
  octave_idx_type VN_columns = args(4).columns();
  octave_idx_type FN_columns = args(5).columns();
  // ensure that there are at least 3 vertices and one face in the mesh and
  // vertex normals and face normals matrices are Nx3.
  if (V_rows < 3)
  {
    std::cout << "There should be at least 3 vertices in the mesh.\n";
    return octave_value_list();
  }
  if (V_columns != 3 && V_columns != 6)
  {
    cout << "Vertex matrix should be Nx3 containing x,y,z coordinates\n";
    cout << "or Nx6 containing vertex coordinates and r,g,b color information.\n";
    return octave_value_list();
  }
  if (F_rows < 1)
  {
    cout << "There should be at least 1 face in the mesh.\n";
    return octave_value_list();
  }
  if (F_columns != 3)
  {
    cout << "Face matrix should be an Nx3 matrix.\n";
    return octave_value_list();
  }
  // check that texture coordinates are an Nx2 matrix
  if (VT_columns != 2)
  {
    cout << "Texture coordinates should be an Nx2 matrix.\n";
    return octave_value_list();
  }
  // check that texture faces are an Nx3 matrix
  if (FT_columns != 3)
  {
    cout << "Texture faces matrix should be an Nx3 matrix.\n";
    return octave_value_list();
  }
  // check that vertex normals are an Nx3 matrix
  if (VN_columns != 3)
  {
    cout << "Vertex normals should be an Nx3 matrix.\n";
    return octave_value_list();
  }
  // check that face normals are an Nx3 matrix
  if (FN_columns != 3)
  {
    cout << "Face normals should be an Nx3 matrix.\n";
    return octave_value_list();
  }
  // check that matrices for faces, texture faces and face normals have
  // the same length to ensure that the iterator will not break when
  // writing to the file
  if (F_rows != FT_rows && F_rows != FN_rows)
  {
    cout << "Faces, texture faces and face normals should be the same length.\n";
    return octave_value_list();
  }
  // check if filename exists
  if (! over)
  {
    filename = checkFilename (filename);
  }
  // open file
  ofstream outputFile(filename.c_str());
  if (! outputFile.is_open())
  {
    cout << "Error opening " << filename.c_str() << "for write\n";
    return octave_value_list();
  }
  else
  {
    // writing header to file
    outputFile << "#\n# OBJ File generated by GNU Octave\n# using 'writeObj' function\n";
    outputFile << "#\n# Object " << filename.c_str() << "\n#\n";
    outputFile << "# Vertices: " << V_rows << "\n";
    outputFile << "# Faces: " << F_rows << "\n#\n#\n";
    int i = filename.length() - 3;
    string mtlfilename = filename.c_str();
    outputFile << "mtllib ./" << mtlfilename.replace(i,3, "mtl") << "\n\n";
    if (info) { cout << "Writing to file... "; }
    // write vertices to file
    if (V_columns == 3) // coordinates only
    {
      float tmpx, tmpy, tmpz;
      for (octave_idx_type i = 0; i < V_rows; i++)
      {
        tmpx = V(i,0);
        tmpy = V(i,1);
        tmpz = V(i,2);
        outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << "\n";
      }
    }
    else                // coordinates and color
    {
      float tmpx, tmpy, tmpz, tmpR, tmpG, tmpB;
      for (octave_idx_type i = 0; i < V_rows; i++)
      {
        tmpx = V(i,0);
        tmpy = V(i,1);
        tmpz = V(i,2);
        tmpR = V(i,3);
        tmpG = V(i,4);
        tmpB = V(i,5);
        outputFile << "v " << tmpx << " " << tmpy << " " <<  tmpz << " ";
        outputFile << tmpR << " " << tmpG << " " <<  tmpB << "\n";
      }
    }
    // write texture coordinates to file
    float tmpu, tmpv;
    for (octave_idx_type i = 0; i < VT_rows; i++)
    {
      tmpu = VT(i,0);
      tmpv = VT(i,1);
      outputFile << "vt " << tmpu << " " << tmpv << "\n";
    }
    // write vertex normals to file
    float tmpx, tmpy, tmpz;
    for (octave_idx_type i = 0; i < VN_rows; i++)
    {
      tmpx = VN(i,0);
      tmpy = VN(i,1);
      tmpz = VN(i,2);
      outputFile << "vn " << tmpx << " " << tmpy << " " << tmpz << "\n";
    }
    // write faces, texture faces and face normals to file
    // write faces and face normals to file
    int tmpa, tmpb, tmpc, tmpd, tmpe, tmpf, tmpg, tmph, tmpi;
    for (octave_idx_type i = 0; i < F_rows; i++)
    {
      tmpa = F(i,0);
      tmpb = F(i,1);
      tmpc = F(i,2);
      tmpd = FT(i,0);
      tmpe = FT(i,1);
      tmpf = FT(i,2);
      tmpg = FN(i,0);
      tmph = FN(i,1);
      tmpi = FN(i,2);
      outputFile << "f " << tmpa << "/" << tmpd << "/" << tmpg << " "
                         << tmpb << "/" << tmpe << "/" << tmph << " "
                         << tmpc << "/" << tmpf << "/" << tmpi << "\n";
    }
    // close file
    outputFile.close();
    if (info) { cout << "done!\n"; }
  }
    if (info)
    {
      cout << "Mesh filename is " << filename.c_str() << "\n";
      cout << "Mesh has " << V_rows << " vertices.\n";
      cout << "Mesh has " << F_rows << " faces.\n";
    }
  // if output argument is requested return the filename
  if (nargout == 1)
  {
    return octave_value_list(filename.c_str());
  }
  return octave_value_list();
}

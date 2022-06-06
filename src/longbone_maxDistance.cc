/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>

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
#include <fstream>
#include <vector>
#include <octave/oct.h>

using namespace std;

struct Vector3           //for 3D points
{
  double x, y, z;
};

struct maxDpoints       //for max distance 3D points
{
	Vector3 V1, V2;
};

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

// function for finding the most distant points of a longbone model
maxDpoints longbone_maxDistance(vector<Vector3> Cloud3D)
{
  // get distance between the first two points of the models
  Vector3 maxD_V1 = Cloud3D[0];
  Vector3 maxD_V2 = Cloud3D[1];
  double maxD_V1_V2 = distancePoints3D(maxD_V1, maxD_V2);
  double maxD_V2_V1 = 0;
  double iter = 0;
  // iterate through the model for more distant points until maximum distance
  // is found
  while(maxD_V1_V2 != maxD_V2_V1 && iter < 5)
  {
    iter++;
    for(int i = 0; i < Cloud3D.size(); ++i)
    {
      double pairD = distancePoints3D(maxD_V1, Cloud3D[i]);
      if(pairD > maxD_V1_V2)
      {
        maxD_V1_V2 = pairD;
        maxD_V2 = Cloud3D[i];
      }
    }
    for(int i = 0; i < Cloud3D.size(); ++i)
    {
      double pairD = distancePoints3D(maxD_V2, Cloud3D[i]);
      if(pairD > maxD_V2_V1)
      {
        maxD_V2_V1 = pairD;
        maxD_V1 = Cloud3D[i];
      }
    }
  }
  // save max distance extreme points to returning structure
  maxDpoints maxD;
  maxD.V1 = maxD_V1;
  maxD.V2 = maxD_V2;
  return maxD;
}

DEFUN_DLD (longbone_maxDistance, args, nargout,
          "-*- texinfo -*-\n\
@deftypefn{Function} [@var{maxDistance}, @var{maxD_v1}, @var{maxD_v1}] = longbone_maxDistance (@var{v})\n\
\n\
\n\
This function calculates the maximum distance of a long bone as represented \
by its mesh vertices.  It requires an Nx3 matrix as input argument, which \
contains the (x,y,z) coordinates in R3.  The function also returns the point \
coordinates that correspond to the maximum distance.\n\
\n\
\n\
The function 'longbone_maxDistance.m' may be called with varying number of \
output arguments.\n\
If one output argument is defined, then maximum distance is returned.\n\
\n\
@var{maxDistance} = longbone_maxDistance(@var{v})\n\
\n\
If two output arguments are defined, then the corresponding points of maximum \
distance are returned.\n\
\n\
[@var{maxD_v1}, @var{maxD_v2}] = longbone_maxDistance(@var{v})\n\
\n\
If 3 output arguments are defined, then the maximum distance along with \
its corresponding points are returned.\n\
\n\
[@var{maxDistance}, @var{maxD_v1}, @var{maxD_v1}] = longbone_maxDistance(@var{v})\n\
@end deftypefn")
{
  // check number of input and output arguments
  if(nargout < 1 || nargout > 3)
  {
    cout << "Invalid number of output arguments.\n";
    return octave_value_list();
  }
  if(args.length() != 1)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  // store vertices in data structure
  Matrix V = args(0).array_value();
	// store 3d point cloud to std::vector
  vector<Vector3> Cloud3D;
  double tmpx, tmpy, tmpz;
  for(octave_idx_type i = 0; i < args(0).rows(); i++)
  {
    tmpx = V(i,0);
    tmpy = V(i,1);
    tmpz = V(i,2);
    Vector3 temp_points = {tmpx, tmpy, tmpz};
    Cloud3D.push_back(temp_points);
  }
  // find points corresponding to the max distance of the bone model
  maxDpoints maxD = longbone_maxDistance(Cloud3D);
  Matrix maxD_v1(1, 3);
  Matrix maxD_v2(1, 3);
  maxD_v1(0,0) = maxD.V1.x;
  maxD_v1(0,1) = maxD.V1.y;
  maxD_v1(0,2) = maxD.V1.z;
  maxD_v2(0,0) = maxD.V2.x;
  maxD_v2(0,1) = maxD.V2.y;
  maxD_v2(0,2) = maxD.V2.z;
  // calculate max distance
  double maxDistance;
  maxDistance = distancePoints3D(maxD.V1, maxD.V2);
  // return values according to number of output variables
  octave_value_list retval;
  if(nargout == 1)
  {
    retval(0) = maxDistance;
  }
  else if(nargout == 2)
  {
    retval(0) = maxD_v1;
    retval(1) = maxD_v2;
  }
  else if(nargout == 3)
  {
    retval(0) = maxDistance;
    retval(1) = maxD_v1;
    retval(2) = maxD_v2;
  }
  return octave_value(retval);
}
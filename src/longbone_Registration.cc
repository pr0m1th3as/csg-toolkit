/*
Copyright (C) 2021-2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>

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
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <bits/stdc++.h>
#include <octave/oct.h>
#include <octave/parse.h>

using namespace std;

struct Vector2           //for 2D points
{
  double x, y;
};
struct Vector3           //for 3D points
{
  double x, y, z;
};
struct Mesh             //for triangular mesh object data
{
  double V1x, V1y, V1z;
  double V2x, V2y, V2z;
  double V3x, V3y, V3z;
};
struct maxDpoints       //for max distance 3D points
{
  Vector3 V1, V2;
};
struct Raster
{
  bool pixel;
  vector<int> pointer;
};
struct Pixel
{
  int row, col;
};
struct Contour
{
  int row, col;
  int direction;
  vector<int> pointer;
};
struct FCCode
{
  int chain;
  int delta_x;
  int delta_y;
  double delta_t;
};
struct EFDcoef
{
  double a, b, c, d;
};
struct Polygon
{
  double Area;
  Vector2 Centroid;
  double Imax, Imin, theta;
  double Ixy, Ix, Iy;
};
// define namespaces for global variables
namespace global
{
  double resolution = 1;
  string bone, side;
  int epiphysis;
}
namespace PCA
{
  vector<double> PC1{0.93531843007037,0.0183028981666609,0.0369610202429386,
                    -0.0883738237697742,-0.0110708382797718,-0.183380947194812,
                  -0.00141348637990715,-0.00399920897049199,-0.280406254855754,
                   0.00595952325266791,-0.0433463182607595,0.00165435555137897,
                   0.00378463038272828,0.0116475692638028,0.000788228919600664,
                  -0.00043713590973217,0.02924125415881,0.000167059799094242,
                  -0.0149938678383692,0.010248599241707,0.00432520436018409};
  vector<double> PC2{0.177588773121458,-0.00121990747029798,-0.287792734466079,
                  0.867960847284974,-0.000700784087249675,-0.00139263620355797,
                    0.00937293078282894,-0.00749106016039342,0.245542759691789,
                    0.00185153345672356,0.135692334950271,-0.199828570784526,
                  -0.0102100454746126,-0.0246754340040846,0.000547128908593401,
                    0.00625205301220659,-0.086159585310214,0.00462930025019041,
                   0.0723734381719101,-0.00545801757731642,0.00547627334741027};
}
// function for adding two 3D vectors
Vector3 addVectors3D(Vector3 A, Vector3 B)
{
  Vector3 C;
  C.x = A.x + B.x;
  C.y = A.y + B.y;
  C.z = A.z + B.z;
  return C;
}
// function for subtracting two 3D vectors
Vector3 subVectors3D(Vector3 A, Vector3 B)
{
  Vector3 C;
  C.x = A.x - B.x;
  C.y = A.y - B.y;
  C.z = A.z - B.z;
  return C;
}
// function for normalizing a vector
Vector3 normalizeVector3D(Vector3 A)
{
  Vector3 normal;
  double v_len = sqrt(A.x * A.x + A.y * A.y + A.z * A.z);
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
// function for calculating the dot product between two vectors of doubles
double dotProduct(vector<double> A, vector<double> B)
{
	double dotprod = 0;
  for(int i = 0; i < A.size(); ++i)
  {
    dotprod += A[i] * B[i];
  }
	return dotprod;
}
// function for calculating the dot product between two 3D vectors
double dotProduct3D(Vector3 A, Vector3 B)
{
	double dotprod = A.x * B.x + A.y * B.y + A.z * B.z;
	return dotprod;
}
// function for calculating CrossProduct between two 3D vectors
Vector3 crossProduct3D(Vector3 A, Vector3 B)
{
  double x, y, z;
  x = A.y * B.z - A.z * B.y;
  y = A.z * B.x - A.x * B.z;
  z = A.x * B.y - A.y * B.x;
  Vector3 crossprod = {x, y, z};
  return crossprod;
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
// function for translating a vector of 3D points from a given centroid to origin
vector<Vector3> translatePoints3D(vector<Vector3> cloud, Vector3 centroid)
{
  for (int idx = 0; idx < cloud.size(); ++idx)
  {
    cloud[idx].x = cloud[idx].x - centroid.x;
    cloud[idx].y = cloud[idx].y - centroid.y;
    cloud[idx].z = cloud[idx].z - centroid.z;
  }
  return cloud;
}
// function for rotating a 3D point cloud about an arbitrary axis by an angle theta
vector<Vector3> rotatePoints3D(vector<Vector3> Cloud3D, Vector3 rot_axis, double theta)
{
  vector<Vector3> rCloud3D;
  double c = cos(theta);
  double c_ = 1 - c;
  double s = sin(theta);
  for(int i = 0; i < Cloud3D.size(); ++i)
  {
    Vector3 tmp1 = crossProduct3D(rot_axis, Cloud3D[i]);
    double tmp2 = dotProduct3D(rot_axis, Cloud3D[i]);
    double x = Cloud3D[i].x * c + tmp1.x * s + rot_axis.x * tmp2 * c_;
    double y = Cloud3D[i].y * c + tmp1.y * s + rot_axis.y * tmp2 * c_;
    double z = Cloud3D[i].z * c + tmp1.z * s + rot_axis.z * tmp2 * c_;
    rCloud3D.push_back({x, y, z});
  }
  return rCloud3D;
}
// function for comparing the first elements between two Pair vectors
bool compareVectorPair(const pair<double, int> A, const pair<double, int> B)
{
  return (A.first < B.first);
}
// function for calculating distance between two 2D points
double distancePoints2D(Vector2 A, Vector2 B)
{
  double distance;
  double x, y;
  x = (A.x - B.x) * (A.x - B.x);
  y = (A.y - B.y) * (A.y - B.y);
  distance = sqrt(x + y);
  return distance;
}
// function for translating a vector of 2D points from a given centroid to origin
vector<Vector2> translatePoints2D(vector<Vector2> polyline2D, Vector2 centroid)
{
  for (int idx = 0; idx < polyline2D.size(); ++idx)
  {
    polyline2D[idx].x = polyline2D[idx].x - centroid.x;
    polyline2D[idx].y = polyline2D[idx].y - centroid.y;
  }
  return polyline2D;
}
// function for rotating a vector of 2D points around the origin by an angle theta
vector<Vector2> rotatePoints2D(vector<Vector2> polyline2D, double theta)
{
  double c = cos(theta);
  double s = sin(theta);
  vector<Vector2> rPoints2D;
  for(int i = 0; i < polyline2D.size(); ++i)
  {
    double x = polyline2D[i].x * c - polyline2D[i].y * s;
    double y = polyline2D[i].y * c + polyline2D[i].x * s;
    rPoints2D.push_back({x, y});
  }
  return rPoints2D;
}
// function for sorting the points of a 2D polygon in counter clockwise order
vector<Vector2> sortPolygon2D(vector<Vector2> polygon)
{
  vector<Vector2> sPoints;
  // find centroid of polygon's vertices
  Vector2 centroid = {0, 0};
  for(int i = 0; i < polygon.size(); ++i)
  {
    centroid.x += polygon[i].x;
    centroid.y += polygon[i].y;
  }
  centroid.x = centroid.x / polygon.size();
  centroid.y = centroid.y / polygon.size();
  // calculale angles of each point with respect to the centroid
  vector<pair<double, int>> angle;
  for(int i = 0; i < polygon.size(); ++i)
  {
    double theta = atan2(polygon[i].y - centroid.y, polygon[i].x - centroid.x);
    angle.push_back({theta, i});
  }
  // sort rows of anglePoints vector from minimum to maximum angle
  sort(angle.begin(), angle.end(), compareVectorPair);
  // push 2D points sorted in counter clockwise order to the returning vector
  for(int i = 0; i < angle.size(); ++i)
  {
    sPoints.push_back({polygon[angle[i].second].x, polygon[angle[i].second].y});
  }
  // iterate through the points and check for nearest neighbor between each
  // pair of consecutive points to compensate for concave polygons that might
  // produce strange polygon boundaries due to the initial arctan sorting
  vector<int> c_index;
  // create a duplicate index vector for the points to be ordered
  for(int i = 0; i < 2; ++i)
  {
    for(int idx = 0; idx < sPoints.size(); ++idx)
    {
      c_index.push_back(idx);
    }
  }
  // for each point
  for(int i = 0; i < sPoints.size(); ++i)
  {
    // measure distance to next point
    Vector2 c_point = sPoints[i];
    Vector2 n_point = sPoints[c_index[i+1]];
    double dist_next = distancePoints2D(c_point, n_point);
    // create a comparison vector with 80% of the points ahead of the next point
    vector<Vector2> compare_vector;
    for(int idx = i + 2; idx < sPoints.size() * 0.8 + i; ++idx)
    {
      compare_vector.push_back(sPoints[c_index[idx]]);
    }
    // find minimum distance from current point to points lying ahead
    double minDist = distancePoints2D(compare_vector[0], c_point);
    int minDist_idx = 0;
    for(int cv_i = 1; cv_i < compare_vector.size(); ++cv_i)
    {
      double dist = distancePoints2D(compare_vector[cv_i], c_point);
      if(dist < minDist)
      {
        minDist_idx = cv_i;
        minDist = dist;
      }
    }
    Vector2 minDistPoint = sPoints[c_index[i+2+minDist_idx]];
    // check if any point ahead is closer than next point and swap
    if(minDist < dist_next)
    {
      sPoints[c_index[i+1]] = minDistPoint;
      sPoints[c_index[i+2+minDist_idx]] = n_point;
    }
  }
  return sPoints;
}
// function for calculating the centroid of a cross section in 3D space
Vector3 centroidPolygon3D(vector<Vector3> cloud, Vector3 normal)
{
  // find most distant point from the first point in vector of 3D points
  Vector3 point0 = cloud[0];
  double maxD = 0;
  int maxDp_i;
  for(int i = 1; i < cloud.size(); ++i)
  {
    double dist = distancePoints3D(point0, cloud[i]);
    if(dist > maxD)
    {
      maxD = dist;
      maxDp_i = i;
    }
  }
  // calculate the normals of the local x and y 2D coordinates
  Vector3 X2D_vector = subVectors3D(cloud[maxDp_i], point0);
  Vector3 Y2D_vector = crossProduct3D(normal, X2D_vector);
  X2D_vector = normalizeVector3D(X2D_vector);
  Y2D_vector = normalizeVector3D(Y2D_vector);
  // transform cross section to 2D coordinates
  vector<Vector2> XYproj;
  for(int i = 0; i < cloud.size(); i++)
  {
    Vector3 new_point = subVectors3D(cloud[i], point0);
    double x2D = dotProduct3D(new_point, X2D_vector);
    double y2D = dotProduct3D(new_point, Y2D_vector);
    XYproj.push_back({x2D, y2D});
  }
  // sort points of 2D polygon in counter clockwise order
  XYproj = sortPolygon2D(XYproj);
  // close the 2D polygon by appending the first point to the end
  XYproj.push_back(XYproj[0]);
  // calculate the centroid of the polygon
  double Cx = 0;
  double Cy = 0;
  double Area = 0;
  for(int i = 0; i < XYproj.size() - 1; ++i)
  {
    double tmp = XYproj[i].x * XYproj[i+1].y - XYproj[i+1].x * XYproj[i].y;
    Cx += (XYproj[i].x + XYproj[i+1].x) * tmp;
    Cy += (XYproj[i].y + XYproj[i+1].y) * tmp;
    Area += tmp;
  }
  Area = 0.5 * Area;
  Cx = (1 / (6 * Area)) * Cx;
  Cy = (1 / (6 * Area)) * Cy;
  // transform 2D centroid to original 3D coordinates
  double x = point0.x + Cx * X2D_vector.x + Cy * Y2D_vector.x;
  double y = point0.y + Cx * X2D_vector.y + Cy * Y2D_vector.y;
  double z = point0.z + Cx * X2D_vector.z + Cy * Y2D_vector.z;
  Vector3 centroid = {x, y, z};
  return centroid;
}
// function for calculating the centroid and 2nd moments of area of a simple 
// polygon. The polygon should be closed, i.e. the first and last points of
// the vector are the same.
Polygon simplePolygon2D(vector<Vector2> poly2D)
{
  // declare return struct variable
  Polygon GEOM;
  // calculate the centroid of the polygon
  double Cx = 0;
  double Cy = 0;
  double Area = 0;
  for(int i = 0; i < poly2D.size() - 1; ++i)
  {
    double tmp = poly2D[i].x * poly2D[i+1].y - poly2D[i+1].x * poly2D[i].y;
    Cx += (poly2D[i].x + poly2D[i+1].x) * tmp;
    Cy += (poly2D[i].y + poly2D[i+1].y) * tmp;
    Area += tmp;
  }
  GEOM.Area = 0.5 * Area;
  Cx = (1 / (6 * Area)) * Cx;
  Cy = (1 / (6 * Area)) * Cy;
  GEOM.Centroid = {Cx, Cy};
  // move polygon's centroid to origin
  for(int i = 0; i < poly2D.size(); ++i)
  {
    poly2D[i].x = poly2D[i].x - Cx;
    poly2D[i].y = poly2D[i].y - Cy;
  }
  // calculate the parameters of second moments of area
  double Ix = 0;
  double Iy = 0;
  double Ixy = 0;
  for(int i = 0; i < poly2D.size() - 1; ++i)
  {
    double tmp = (poly2D[i].x * poly2D[i+1].y) - (poly2D[i+1].x * poly2D[i].y);
    Ix += ((poly2D[i].y * poly2D[i].y) + (poly2D[i].y * poly2D[i+1].y) +
          (poly2D[i+1].y * poly2D[i+1].y)) * tmp; 
    Iy += ((poly2D[i].x * poly2D[i].x) + (poly2D[i].x * poly2D[i+1].x) +
          (poly2D[i+1].x * poly2D[i+1].x)) * tmp;
    Ixy += ((poly2D[i].x * poly2D[i+1].y) + (2 * poly2D[i].x * poly2D[i].y) +
    (2 * poly2D[i+1].x * poly2D[i+1].y) + (poly2D[i+1].x * poly2D[i].y)) * tmp;
  }
  Ix = Ix / 12;
  Iy = Iy / 12;
  Ixy = Ixy / 24;
  double tmp = sqrt((((Ix - Iy) * 0.5) * ((Ix - Iy) * 0.5)) + (Ixy * Ixy));
  GEOM.Imin = ((Ix + Iy) * 0.5) - tmp;
  GEOM.Imax = ((Ix + Iy) * 0.5) + tmp;
  // check for quadrant and correct angle of theta with respect to the
  // positive side of the x axis
  double PI = 3.14159265358979323846;
  if(Ix >= Iy)
  {
    GEOM.theta = 0.5 * atan(-Ixy / ((Ix - Iy) / 2));
  }
  else if((Ix < Iy) && (Ixy < 0))
  {
    GEOM.theta = 0.5 * atan(-Ixy / ((Ix - Iy) / 2)) + (PI * 0.5);
  }
  else
  {
    GEOM.theta = 0.5 * atan(-Ixy / ((Ix - Iy) / 2)) - (PI * 0.5);
  }
  GEOM.Ixy = Ixy;
  GEOM.Ix = Ix;
  GEOM.Iy = Iy;
  return GEOM;
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
// function for extracting intersection points ginen a plane specified
// by a normal and a point
vector<Vector3> sliceMesh(vector<Mesh> Mesh3D, Vector3 point, Vector3 normal)
{
	// initialize vectors and variables
	vector<Vector3> iPoints;
	Vector3 V1, V2, V3;
	double dotV1, dotV2, dotV3;
	// for each triplet of vertices find which face is intersected by specified plane
	for(int i = 0; i < Mesh3D.size(); ++i)
	{
		V1 = {Mesh3D[i].V1x - point.x, Mesh3D[i].V1y - point.y, Mesh3D[i].V1z - point.z};
		V2 = {Mesh3D[i].V2x - point.x, Mesh3D[i].V2y - point.y, Mesh3D[i].V2z - point.z};
		V3 = {Mesh3D[i].V3x - point.x, Mesh3D[i].V3y - point.y, Mesh3D[i].V3z - point.z};
		dotV1 = dotProduct3D(V1, normal);
		dotV2 = dotProduct3D(V2, normal);
		dotV3 = dotProduct3D(V3, normal);
    // true if intersecting
		if(! ((dotV1 < 0 && dotV2 < 0 && dotV3 < 0) || 
          (dotV1 > 0 && dotV2 > 0 && dotV3 > 0)))
		{
			if(dotV1 * dotV2 < 0)       // vertices 1 and 2 lie on opposite sides
			{
				Vector3 PL = {point.x - Mesh3D[i].V1x, point.y - Mesh3D[i].V1y,
                      point.z - Mesh3D[i].V1z};
				Vector3 L = {Mesh3D[i].V2x - Mesh3D[i].V1x, Mesh3D[i].V2y - Mesh3D[i].V1y,
                     Mesh3D[i].V2z - Mesh3D[i].V1z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv1v2 = {Mesh3D[i].V1x + d * L.x, Mesh3D[i].V1y + d * L.y,
                          Mesh3D[i].V1z + d * L.z};
				iPoints.push_back(CSv1v2);
			}
			if(dotV1 * dotV3 < 0)       // vertices 1 and 3 lie on opposite sides
			{
				Vector3 PL = {point.x - Mesh3D[i].V1x, point.y - Mesh3D[i].V1y,
                      point.z - Mesh3D[i].V1z};
				Vector3 L = {Mesh3D[i].V3x - Mesh3D[i].V1x, Mesh3D[i].V3y - Mesh3D[i].V1y,
                     Mesh3D[i].V3z - Mesh3D[i].V1z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv1v3 = {Mesh3D[i].V1x + d * L.x, Mesh3D[i].V1y + d * L.y,
                          Mesh3D[i].V1z + d * L.z};
				iPoints.push_back(CSv1v3);
			}
			if(dotV2 * dotV3 < 0)       // vertices 2 and 3 lie on opposite sides
			{
				Vector3 PL = {point.x - Mesh3D[i].V2x, point.y - Mesh3D[i].V2y,
                      point.z - Mesh3D[i].V2z};
				Vector3 L = {Mesh3D[i].V3x - Mesh3D[i].V2x, Mesh3D[i].V3y - Mesh3D[i].V2y,
                     Mesh3D[i].V3z - Mesh3D[i].V2z};
				double d = dotProduct3D(PL, normal) / dotProduct3D(L, normal);
				Vector3 CSv2v3 = {Mesh3D[i].V2x + d * L.x, Mesh3D[i].V2y + d * L.y,
                          Mesh3D[i].V2z + d * L.z};
				iPoints.push_back(CSv2v3);
			}
			if(dotV1 == 0)      // vertex 1 lies on the cross-section plane
			{
				Vector3 CSv1;
				CSv1 = {Mesh3D[i].V1x, Mesh3D[i].V1y, Mesh3D[i].V1z};
				iPoints.push_back(CSv1);
			}
			if(dotV2 == 0)      // vertex 2 lies on the cross-section plane
			{
				Vector3 CSv2;
				CSv2 = {Mesh3D[i].V2x, Mesh3D[i].V2y, Mesh3D[i].V2z};
				iPoints.push_back(CSv2);
			}
			if(dotV3 == 0)      // vertex 3 lies on the cross-section plane
			{
				Vector3 CSv3;
				CSv3 = {Mesh3D[i].V3x, Mesh3D[i].V3y, Mesh3D[i].V3z};
				iPoints.push_back(CSv3);
			}
		}
	}
  // remove duplicate intersected points
  vector<Vector3> uPoints;
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
	while(!duplicate_found)
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
	return uPoints;
}
// function for rasterizing the 2D points onto a b&w image.
vector<vector<Raster>> rasterizeCloud2D(vector<Vector2> pCloud)
{
  // find max edge limits in x and y axis
  double left_side = 0;
  double right_side = 0;
  double upper_side = 0;
  double lower_side = 0;
  for(vector<Vector2>::iterator it = pCloud.begin(); it != pCloud.end(); ++it)
  {
    if(it->x < left_side)
    {
      left_side = it->x;
    }
    if(it->x > right_side)
    {
      right_side = it->x;
    }
    if(it->y < lower_side)
    {
      lower_side = it->y;
    }
    if(it->y > upper_side)
    {
      upper_side = it->y;
    }
  }
  // find maximum horizontal and vertical dimensions
  double hor_max_dist = right_side - left_side;
  double ver_max_dist = upper_side - lower_side;
  // get number of necessary pixels for x and y axis
  int hor_pixels = int(hor_max_dist / global::resolution) + 8;
  int ver_pixels = int(ver_max_dist / global::resolution) + 8;
  // determine left and top image boundaries
  double left_boundary = left_side - (4 * global::resolution);
  double top_boundary = upper_side + (4 * global::resolution);
  // create raster image of suitable size
  vector<vector<Raster>> rasterMatrix;
  rasterMatrix.resize(ver_pixels, vector<Raster>(hor_pixels, {false, {-1}}));
  // map 2D points to raster image
  for(int it = 0; it < pCloud.size(); ++it)
  {
    // measure distance from left and top boundaries and calculate the
    // respective pixel coordinates from the top left corner of the image
    double dist_2_left = abs(left_boundary) + pCloud[it].x;
    double dist_2_top = top_boundary - pCloud[it].y;
    int row_index = ceil(dist_2_top / global::resolution);
    int col_index = ceil(dist_2_left / global::resolution);
    // if pixel is not already a valid boundry pixel then set it and write the
    // pointer of the 2D point being checked in place of the default -1 value
    if(!rasterMatrix[row_index][col_index].pixel)
    {
        rasterMatrix[row_index][col_index].pixel = true;
        rasterMatrix[row_index][col_index].pointer[0] = {it};
    }
    // if pixel is already true then just push_back the pointer of the 2D point
    else
    {
        rasterMatrix[row_index][col_index].pointer.push_back(it);
    }
  }
  return rasterMatrix;
}
// function for calculating the Moore Neighbor algorithm on a b&w image.
vector<Contour> calculateMooreNeighbor(vector<vector<Raster>> rasterMatrix)
{
  // define required variables
  vector<Contour> boundaryPixels;
  // populate Neighbor vectors with row,column directions
  // directions in rasterMatrix       [row], [column]
  //              // top -              -1,       0
  //              // top -left          -1,      -1
  //              //     -left           0,      -1
  //              /bottom-left           1,      -1
  //              /bottom-               1,       0
  //              /bottom-right          1,       1
  //              //     -right          0,       1
  //              // top -right         -1,       1
  struct MooreNeighbor
  {
      vector<int> row = {-1, -1, 0, 1, 1, 1, 0, -1};
      vector<int> col = {0, -1, -1, -1, 0, 1, 1, 1};
  };
  MooreNeighbor Mp;
  int direction = 4;
  Pixel b;                            //backtrack of current boundary pixel
  // find an active pixel by starting the search from the top middle of image
  int mid = rasterMatrix[0].size() / 2;
  Pixel p = {0, mid};
  bool pixel_found = false;
  while(!pixel_found)
  {
    p.row++;
    if(rasterMatrix[p.row][p.col].pixel)  //check for active pixels
    {
      b.col = p.col;
      b.row = p.row;
      pixel_found = true;
    }
  }
  vector<int> pointers = rasterMatrix[p.row][p.col].pointer;
  Contour valid_pixel = {p.row, p.col, direction, pointers};
  boundaryPixels.push_back(valid_pixel);
  // start loop until the contour is closed or all pixels are searched
  int it = 0;
  bool contour_open = true;
  int max_it = rasterMatrix.size() * rasterMatrix[0].size();
  while(contour_open && it < max_it)
  {
    //inverse direction for the next search
    int dir = boundaryPixels[it].direction + 4;
    // start from initial searching direction based on the previous valid pixel
    Pixel c;                          //current pixel under consideration
    bool dir_found = false;
    while(!dir_found)
    {
      if(rasterMatrix[p.row + Mp.row[(dir+1)%8]][p.col + Mp.col[(dir+1)%8]].pixel)
      //next pixel found!
      {
        //save indices on the validated pixel under consideration
        c.row = p.row + Mp.row[(dir+1)%8];
        c.col = p.col + Mp.col[(dir+1)%8];
        //check found flag to stop rotating aroung the Neighborhood
        dir_found = true;
        b = p;
        p = c;
        // append the pixel's row and column indices, direction, and associated
        // pointers in the returning vector structure
        direction = (dir+1)%8;
        vector<int> pointers = rasterMatrix[c.row][c.col].pointer;
        Contour valid_pixel = {c.row, c.col, direction, pointers};
        boundaryPixels.push_back(valid_pixel);
      }
      dir++;
    }
    // check if the contour has been closed. This test is initialized after
    // the second loop. in order to compare that the closing pixel with the
    // second found along with the direction from
    // initial pixel towards the second one in the contour.
    if(boundaryPixels.size() > 2)
    {
      int row_0 = boundaryPixels[0].row;
      int col_0 = boundaryPixels[0].col;
      
      int row_1 = boundaryPixels[1].row;
      int col_1 = boundaryPixels[1].col;
      int dir_1 = boundaryPixels[1].direction;
      
      int size = boundaryPixels.size() - 1;
      int row_l = boundaryPixels[size].row;
      int col_l = boundaryPixels[size].col;
      int dir_l = boundaryPixels[size].direction;
      
      int row_bl = boundaryPixels[size - 1].row;
      int col_bl = boundaryPixels[size - 1].col;
      int dir_bl = boundaryPixels[size - 1].direction;
      // the contour is closed if current pixel (last found) matches the second
      // pixel in the vector (position and direction) and the previously found
      // pixel's position matches that of the initial pixel
      if(row_0 == row_bl && col_0 == col_bl && row_1 == row_l && 
                            col_1 == col_l && dir_1 == dir_l)
      {
        // in such case copy the direction of the previously found pixel
        // to its instance at the beginning of the vector (initial pixel)
        // and erase the last element of the vector (also second pixel)
        // so that the contour forms a closed polygon vector.
        boundaryPixels[0].direction = dir_bl;
        boundaryPixels.erase(boundaryPixels.end() - 1);
        // assert the closed contour flag to break the loop
        contour_open = false;
      }
    }
    it++;
  }
  return boundaryPixels;
}
// function for extracting 2D polyline from the boundary pixels vector
vector<Vector2> extract2Dpolyline(vector<Contour> boundaryPixels, vector<Vector2> pCloud)
{
  // initialize variables
  vector<Vector2> Polyline;
  int poly_length = boundaryPixels.size();
  // start iterating through the boundaryPixels vector
  for(int pixel_it = 0; pixel_it < poly_length; ++pixel_it)
  {
    // find the corresponding 2D points from the planar cross-section
    int num_pointers = boundaryPixels[pixel_it].pointer.size();
    // initialize 2D centroid for each pixel
    Vector2 p_centroid = {0, 0};
    for(int pointer_it = 0; pointer_it < num_pointers; ++pointer_it)
    {
      p_centroid.x += pCloud[boundaryPixels[pixel_it].pointer[pointer_it]].x;
      p_centroid.y += pCloud[boundaryPixels[pixel_it].pointer[pointer_it]].y;
    }
    p_centroid.x = p_centroid.x / num_pointers;
    p_centroid.y = p_centroid.y / num_pointers;
    Polyline.push_back(p_centroid);
  }
  return Polyline;
}
// function for extracting Freeman chain code from contour vectors.
// Since the chain coding in Moore's Neighbor algorithm does not use the
// standard direction pattern (0 points upwards instead towards the right),
// the present function makes the necessary shift so that the final chain 
// code follows the standard implementation of Freeman chain coding (0 points to
// the positive side of x axis and coding follows a counter clockwise direction)
// Furthermore, the respective displacement along the x and y axes are calculated
// along with the length of each step.
vector<FCCode> extractFreemanChainCode(vector<Contour> boundaryPixels)
{
  vector<FCCode> chaincode;
  // precalulcate displacements
  double straight = global::resolution;
  double diagonal = sqrt(straight * straight + straight * straight);
  // start with the 2nd element of the segment
  for(int it = 1; it < boundaryPixels.size(); ++it)
  {
    int FCCode = (boundaryPixels[it].direction + 2) % 8;
    double displacement;
    int x, y;
    if(FCCode == 0)
    {
      x = 1; y = 0; displacement = straight;
    }
    else if(FCCode == 1)
    {
      x = 1; y = 1; displacement = diagonal;
    }
    else if(FCCode == 2)
    {
      x = 0; y = 1; displacement = straight;
    }
    else if(FCCode == 3)
    {
     x = -1; y = 1; displacement = diagonal;
    }
    else if(FCCode == 4)
    {
      x = -1; y = 0; displacement = straight;
    }
    else if(FCCode == 5)
    {
      x = -1; y = -1; displacement = diagonal;
    }
    else if(FCCode == 6)
    {
      x = 0; y = -1; displacement = straight;
    }
    else
    {
      x = 1; y = -1; displacement = diagonal;
    }
    chaincode.push_back({FCCode, x, y, displacement});
  }    
  return chaincode;
}
// function for calculating Elliptic Fourier Descriptors
vector<EFDcoef> calculaleEFD(vector<FCCode> FreemanChainCode)
{
  vector<EFDcoef> EFDs;
  double PI = 3.14159265358979323846;
  // calculale the DC components as well as total displacement T along the contour
  double A_0 = 0;
  double C_0 = 0;
  double T = 0;
  double zeta_A = 0;
  double zeta_B = 0;
  double deltaA = 0;
  double deltaB = 0;
  double t_c = 0;
  double t_prev = 0;
  for(int p = 0; p < FreemanChainCode.size(); ++p)
  {
    int delta_x = FreemanChainCode[p].delta_x;
    int delta_y = FreemanChainCode[p].delta_y;
    double delta_t = FreemanChainCode[p].delta_t;
    t_c += delta_t; 
    A_0 += (delta_x / (2 * delta_t)) * (t_c * t_c - t_prev * t_prev);
    A_0 += (zeta_A - (delta_x / delta_t) * zeta_B) * (t_c - t_prev);
    C_0 += (delta_y / (2 * delta_t)) * (t_c * t_c - t_prev * t_prev);
    C_0 += (deltaA - (delta_y / delta_t) * deltaB) * (t_c - t_prev);
    t_prev = t_c;
    if(p > 0)
    {
      zeta_A += FreemanChainCode[p - 1].delta_x;
      zeta_B += FreemanChainCode[p - 1].delta_t;
      deltaA += FreemanChainCode[p - 1].delta_y;
      deltaB += FreemanChainCode[p - 1].delta_t;
    }
    T += FreemanChainCode[p].delta_t;
  }
  A_0 = A_0 / T;
  C_0 = C_0 / T;
  EFDs.push_back({A_0, 0, C_0, 0});
  // calculate the EFD coefficients for a number of harmonics up to half the
  // number of links present in the contour
  // for each harmonic calulcate also the respective Fourier power
  vector<EFDcoef> harmonics;
  vector<double> H_power;
  double power = 0;
  for(int n = 1; n < 1 + (FreemanChainCode.size() / 2); ++n)
  {
    // calculate constant terms for each harmonic
    double T_2npi = T / (2 * n * n * PI * PI);
    double pi2n_T = (2 * PI * n) / T;
    double A_n = 0;
    double B_n = 0;
    double C_n = 0;
    double D_n = 0;
    double t_c = 0;         //accumulated length until current point
    double t_prev = 0;      //accumulated length until previous point
    for(int p = 0; p < FreemanChainCode.size(); ++p)
    {
      int delta_x = FreemanChainCode[p].delta_x;
      int delta_y = FreemanChainCode[p].delta_y;
      double delta_t = FreemanChainCode[p].delta_t;
      t_c += delta_t;       //add pth links length to accumulated length
      A_n += (delta_x / delta_t) * (cos(pi2n_T * t_c) - cos(pi2n_T * t_prev));
      B_n += (delta_x / delta_t) * (sin(pi2n_T * t_c) - sin(pi2n_T * t_prev));
      C_n += (delta_y / delta_t) * (cos(pi2n_T * t_c) - cos(pi2n_T * t_prev));
      D_n += (delta_y / delta_t) * (sin(pi2n_T * t_c) - sin(pi2n_T * t_prev));
      t_prev = t_c;         //save current accumulated length to previous
    }
    A_n = A_n * T_2npi;
    B_n = B_n * T_2npi;
    C_n = C_n * T_2npi;
    D_n = D_n * T_2npi;
    harmonics.push_back({A_n, B_n, C_n, D_n});
    power += (A_n*A_n + B_n*B_n + C_n*C_n + D_n*D_n) / 2;
    H_power.push_back(power);
  }
  // for every harmonic check if the average cumulative Fourier power
  // is less than 99.99% and append it in the EFDs vector
  for(int i = 0; i < harmonics.size(); ++i)
  {
    if(H_power[i] < 0.9999 * (H_power[H_power.size() - 1]))
    {
      EFDs.push_back(harmonics[i]);
    }
  }
  return EFDs;
}
// function for normalizing Elliptic Fourier Descriptors
vector<EFDcoef> normalizeEFD(vector<EFDcoef> EFDs)
{
  vector<EFDcoef> nEFDs;
  // calculate required coefficients
  double nom = 2 * (EFDs[1].a * EFDs[1].b + EFDs[1].c * EFDs[1].d);
  double den = EFDs[1].a * EFDs[1].a + EFDs[1].c * EFDs[1].c;
  den += - EFDs[1].b * EFDs[1].b - EFDs[1].d * EFDs[1].d;
  double theta = 0.5 * atan2(nom,den);
  double a = EFDs[1].a * cos(theta) + EFDs[1].b * sin(theta);
  double c = EFDs[1].c * cos(theta) + EFDs[1].d * sin(theta);
  double Epsilon = sqrt(a * a + c * c);
  double phi = atan2(c,a);
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  // normalize all available harmonics
  for(int n = 0; n < EFDs.size(); ++n)
  {
    double n_theta = n * theta;
    double An = (cosphi * EFDs[n].a + sinphi * EFDs[n].c) / Epsilon;
    double Bn = (cosphi * EFDs[n].b + sinphi * EFDs[n].d) / Epsilon;
    double Cn = (-sinphi * EFDs[n].a + cosphi * EFDs[n].c) / Epsilon;
    double Dn = (-sinphi * EFDs[n].b + cosphi * EFDs[n].d) / Epsilon;
    double a_n = An * cos(n_theta) + Bn * sin(n_theta);
    double b_n = - An * sin(n_theta) + Bn * cos(n_theta);
    double c_n = Cn * cos(n_theta) + Dn * sin(n_theta);
    double d_n = - Cn * sin(n_theta) + Dn * cos(n_theta);
    // correct for near 0 imprecision due to floating point arithmetic
    if(abs(a_n) < 0.0000000001) {a_n = 0;}
    if(abs(b_n) < 0.0000000001) {b_n = 0;}
    if(abs(c_n) < 0.0000000001) {c_n = 0;}
    if(abs(d_n) < 0.0000000001) {d_n = 0;}
    EFDcoef n_harmD = {a_n, b_n, c_n, d_n};
    nEFDs.push_back(n_harmD);
  }
  return nEFDs;
}
// function for arranging specific EFDs to vector of doubles
vector<double> arrangeEFD(vector<EFDcoef> EFDs)
{
  vector<double> EFDvec;
  EFDvec.push_back(EFDs[1].d);
  for(int i = 2; i < 7; ++i)
  {
    EFDvec.push_back(EFDs[i].a);
    EFDvec.push_back(EFDs[i].b);
    EFDvec.push_back(EFDs[i].c);
    EFDvec.push_back(EFDs[i].d);
  }
  return EFDvec;
}
// function for classifying bone and proximal distal epiphysis
void classifyBone(vector<EFDcoef> EFD_1, vector<EFDcoef> EFD_2, Polygon GEOM_1, Polygon GEOM_2)
{
  // check for Ulna by comparing the Area of the distal epiphysis
  double Area1 = GEOM_1.Area;
  double Area2 = GEOM_2.Area;
  if(Area1 < 500 || Area2 < 500)
  {
    global::bone = "Ulna";
    // proximal epiphysis is the larger one
    if(Area1 > Area2)
    {
      global::epiphysis = 1;
    }
    else
    {
      global::epiphysis = 2;
    }
    return;
  }
  // check for Femur, Tibia, and Humerus
  // declare values for group centroids of PC scores {Vector2(PC1, PC2)}
  // [0]: F_distal_1, [1]: F_distal_2, [2]: F_proximal
  // [3]: T_proximal, [4]: T_distal, [5]: H_distal, [6]: H_proximal
  vector<Vector2> BoneGroup;
  BoneGroup.push_back({0.621, 0.389});
  BoneGroup.push_back({0.664,-0.034});
  BoneGroup.push_back({0.524, 0.151});
  BoneGroup.push_back({0.711, 0.184});
  BoneGroup.push_back({0.779, 0.147});
  BoneGroup.push_back({0.509, 0.123});
  BoneGroup.push_back({0.865, 0.169});
  // arrange EFDs to vectors
  vector<double> EFDvec1 = arrangeEFD(EFD_1);
  vector<double> EFDvec2 = arrangeEFD(EFD_2);
  // calculate PC scores for each epiphysis (1 and 2)
  Vector2 E1, E2;
  E1.x = dotProduct(EFDvec1, PCA::PC1);
  E1.y = dotProduct(EFDvec1, PCA::PC2);
  E2.x = dotProduct(EFDvec2, PCA::PC1);
  E2.y = dotProduct(EFDvec2, PCA::PC2);
  // calculale ratios of Imax/Imin for each epiphysis
  double Ratio1 = GEOM_1.Imax / GEOM_1.Imin;
  double Ratio2 = GEOM_2.Imax / GEOM_2.Imin;
  // find nearest neighbor for each epiphysis
  int bg1, bg2;
  double minD1 = 1;
  double minD2 = 1;
  for(int i = 0; i < BoneGroup.size(); ++i)
  {
    double tmpD1 = distancePoints2D(E1, BoneGroup[i]);
    double tmpD2 = distancePoints2D(E2, BoneGroup[i]);
    if(tmpD1 < minD1)
    {
      minD1 = tmpD1;
      bg1 = i;
    }
    if(tmpD2 < minD2)
    {
      minD2 = tmpD2;
      bg2 = i;
    }
  }
  // compare paired groups of centroids and select appropriate bone, use
  // Imax/Imin ratios to verify epiphysis and Areas to select epiphysis
  // for tibia
  // searching for distal femur
  if((bg1 == 0 || bg1 == 1) && Ratio2 > Ratio1)
  {
    global::bone = "Femur";
    global::epiphysis = 1;
  }
  else if((bg2 == 0 || bg2 == 1) && Ratio1 > Ratio2)
  {
    global::bone = "Femur";
    global::epiphysis = 2;
  }
  // searching for distal humerus
  else if((bg1 == 2 || bg1 == 5) && (bg2 == 4 || bg2 == 6) && Ratio1 > Ratio2)
  {
    global::bone = "Humerus";
    global::epiphysis = 1;
  }
  else if((bg2 == 2 || bg2 == 5) && (bg1 == 4 || bg1 == 6) && Ratio2 > Ratio1)
  {
    global::bone = "Humerus";
    global::epiphysis = 2;
  }
  // searching for distal humerus extreme cases (distinct from tibia)
  else if(bg1 == 3 && E1.x < 0.63 && (bg2 == 3 || bg2 == 4 || bg2 == 6) &&
          Ratio1 > Ratio2 && (GEOM_1.Area / GEOM_2.Area) < 1.6)
  {
    global::bone = "Humerus";
    global::epiphysis = 1;
  }
  else if(bg2 == 3 && E2.x < 0.63 && (bg1 == 3 || bg1 == 4 || bg1 == 6) &&
          Ratio2 > Ratio1 && (GEOM_2.Area / GEOM_1.Area) < 1.6)
  {
    global::bone = "Humerus";
    global::epiphysis = 2;
  }
  // searching for proximal tibia
  else if((bg1 == 3 || bg1 == 4) && (bg2 == 3 || bg2 == 4 || bg2 == 6) &&
          GEOM_1.Area > GEOM_2.Area)
  {
    global::bone = "Tibia";
    global::epiphysis = 1;
  }
  else if((bg2 == 3 || bg2 == 4) && (bg1 == 3 || bg1 == 4 || bg1 == 6) &&
          GEOM_2.Area > GEOM_1.Area)
  {
    global::bone = "Tibia";
    global::epiphysis = 2;
  }
  else
  {
    global::bone = "Undetermined";
    global::epiphysis = 0;
  }
  return;
}
// function for registering points in 2D polyline shape
vector<Vector2> register2Dpolyline(vector<Vector2> polyline2D, Polygon GEOM)
{
  // declare variables
  vector<Vector2> t_Poly2D;
  Vector2 centroid;
  // align Imax to x axis. For Femur, Tibia, and Humerus, Imax approximates
  // the sagittal axis (anteroposterior). For Ulna, Imax approximates the
  // coronal axis (mediolateral)
  t_Poly2D = translatePoints2D(polyline2D, GEOM.Centroid);
  t_Poly2D = rotatePoints2D(t_Poly2D, -GEOM.theta);
  centroid = {- GEOM.Centroid.x, -GEOM.Centroid.y};
  t_Poly2D = translatePoints2D(t_Poly2D, centroid);
  // find extreme points along x and y axes to determine a new horizontal
  // midline and the correct side (left [-x] or right [+x]) to register
  // the two alignment points. The correct side is the one with the most
  // distant extremum from origin along the x axis
  double max_x = 0;
  double min_x = 0;
  double max_y = 0;
  double min_y = 0;
  for(int i = 0; i < t_Poly2D.size(); ++i)
  {
    if(t_Poly2D[i].x > max_x)
    {
      max_x = t_Poly2D[i].x;
    }
    if(t_Poly2D[i].x < min_x)
    {
      min_x = t_Poly2D[i].x;
    }
    if(t_Poly2D[i].y > max_y)
    {
      max_y = t_Poly2D[i].y;
    }
    if(t_Poly2D[i].y < min_y)
    {
      min_y = t_Poly2D[i].y;
    }
  }
  double mid_y = (max_y + min_y) / 2;
  // search for the registration points
  Vector2 point_A = {0, 0};
  Vector2 point_B = {0, 0};
  // for Femur, Tibia, and Humerus find two distinct registration points
  // one above and one below the horizontal line set by mid_y
  if(global::bone == "Femur" || global::bone == "Tibia" ||
     global::bone == "Humerus")
  {
    if(max_x > - min_x)   // on the positive x-axis
    {
      for(int i = 0; i < t_Poly2D.size(); ++i)
      {
        // check for upper side [+y]
        if(t_Poly2D[i].x > point_A.x && t_Poly2D[i].y > mid_y)
        {
          point_A = t_Poly2D[i];
        }
        // check for lower side [-y]
        if(t_Poly2D[i].x > point_B.x && t_Poly2D[i].y < mid_y)
        {
          point_B = t_Poly2D[i];
        }
      }
    }
    else    // on the negative x-axis
    {
      for(int i = 0; i < t_Poly2D.size(); ++i)
      {
        // check for upper side [+y]
        if(t_Poly2D[i].x < point_A.x && t_Poly2D[i].y > mid_y)
        {
          point_A = t_Poly2D[i];
        }
        // check for lower side [-y]
        if(t_Poly2D[i].x < point_B.x && t_Poly2D[i].y < mid_y)
        {
          point_B = t_Poly2D[i];
        }
      }
    }
  }
  // for Ulna find a single registration point as the most distant point from
  // origin on the y axis and nearest to Imin (min|x|)
  else if(global::bone == "Ulna")
  {
    if(max_y > - min_y)   // on the positive y-axis
    {
      for(int i = 0; i < t_Poly2D.size(); ++i)
      {
        // check for closest to Imin
        double x_diff = abs(t_Poly2D[i].x);
        if(t_Poly2D[i].y > 0 && t_Poly2D[i].x < x_diff)
        {
          x_diff = abs(t_Poly2D[i].x);
          point_A = t_Poly2D[i];
        }
      }
    }
    else    // on the negative y-axis
    {
      for(int i = 0; i < t_Poly2D.size(); ++i)
      {
        // check for closest to Imin
        double x_diff = abs(t_Poly2D[i].x);
        if(t_Poly2D[i].y < 0 && t_Poly2D[i].x < x_diff)
        {
          x_diff = abs(t_Poly2D[i].x);
          point_A = t_Poly2D[i];
        }
      }
    }
    point_B = point_A;
  }
  // reverse rotation to align registration points on the original coordinates
  // of the polyline
  vector<Vector2> regPoints2D;
  regPoints2D.push_back(point_A);
  regPoints2D.push_back(point_B);
  regPoints2D = translatePoints2D(regPoints2D, GEOM.Centroid);
  regPoints2D = rotatePoints2D(regPoints2D, GEOM.theta);
  regPoints2D = translatePoints2D(regPoints2D, centroid);
  return regPoints2D;
}
// function for registering points in 3D epiphyseal point cloud
vector<Vector3> register3DepCloud(vector<Vector3> epCloud, vector<Vector2> regPoints2D)
{
  // declare variables
  vector<Vector3> regPoints3D;
  Vector3 point_A = epCloud[0];
  Vector3 point_B = epCloud[0];
  // find points in 3D cloud, whose y and z coordinates are closest to regPoints2D
  Vector2 yz_tmp = {epCloud[0].y, epCloud[0].z};
  double minDist_A = distancePoints2D(yz_tmp, regPoints2D[0]);
  double minDist_B = distancePoints2D(yz_tmp, regPoints2D[1]);
  for(int i = 1; i < epCloud.size(); ++i)
  {
    yz_tmp = {epCloud[i].y, epCloud[i].z};
    double tmpDist_A = distancePoints2D(yz_tmp, regPoints2D[0]);
    double tmpDist_B = distancePoints2D(yz_tmp, regPoints2D[1]);
    if(tmpDist_A < minDist_A)
    {
      minDist_A = tmpDist_A;
      point_A = epCloud[i];
    }
    if(tmpDist_B < minDist_B)
    {
      minDist_B = tmpDist_B;
      point_B = epCloud[i];
    }
  }
  // save registered 3D points in returning vector
  regPoints3D.push_back(point_A);
  regPoints3D.push_back(point_B);
  return regPoints3D;
}

DEFUN_DLD (longbone_Registration, args, nargout,
          "-*- texinfo -*-\n\
@deftypefn{Function} [@var{bone}, @var{Points}] = longbone_Registration (@var{v},@var{f})\n\
\n\
\n\
This function finds the long bone (i.e. Femur, Ulna, Tibia, or Humerus)\
represented in the input mesh defined by its vertices @var{v} and faces @var{f}\
and it also registers the two points on the bone's surface required by the\
CSG-Toolkit for anatomical alignment of the bone. The function returns the long\
bone's name in @var{bone} and if a second argument is present it returms the\
3D coordinates of the registered points.\n\
@end deftypefn")
{
  // check number of input and output arguments
  if(args.length() != 2)
  {
    cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  if(nargout < 1 && nargout > 2)
  {
    cout << "Invalid number of output arguments.\n";
    return octave_value_list();
  }
  // store vertices and faces in Mesh data structure
  Matrix V = args(0).array_value();
  Matrix F = args(1).array_value();
	vector<Mesh> Mesh3D;
  // find number of vertices and faces
  octave_idx_type V_rows = args(0).rows();
  octave_idx_type F_rows = args(1).rows();
	// for each face store its corresponding vertex coordinates
	double tmpV1x, tmpV1y, tmpV1z;
	double tmpV2x, tmpV2y, tmpV2z;
	double tmpV3x, tmpV3y, tmpV3z;
	for(octave_idx_type i = 0; i < F_rows; i++)
  {
    tmpV1x = V(F(i,0) - 1, 0); tmpV1y = V(F(i,0) - 1, 1); tmpV1z = V(F(i,0) - 1, 2);
    tmpV2x = V(F(i,1) - 1, 0); tmpV2y = V(F(i,1) - 1, 1); tmpV2z = V(F(i,1) - 1, 2);
    tmpV3x = V(F(i,2) - 1, 0); tmpV3y = V(F(i,2) - 1, 1); tmpV3z = V(F(i,2) - 1, 2);
    Mesh temp_mesh = {tmpV1x, tmpV1y, tmpV1z, tmpV2x, tmpV2y, tmpV2z,
                                              tmpV3x, tmpV3y, tmpV3z};
    Mesh3D.push_back(temp_mesh);
  }
	// store 3d point cloud to vector
  vector<Vector3> Cloud3D;
  double tmpx, tmpy, tmpz;
  for(octave_idx_type i = 0; i < V_rows; i++)
  {
    tmpx = V(i,0);
    tmpy = V(i,1);
    tmpz = V(i,2);
    Vector3 temp_points = {tmpx, tmpy, tmpz};
    Cloud3D.push_back(temp_points);
  }
  // find points corresponding to the max distance of the bone model
  maxDpoints maxD = longbone_maxDistance(Cloud3D);
  // calculate the normal vector and the 2 points in R3 that define the slicing
  // planes at 20, 30, 70, and 80% of the maximum length of the bone
  Vector3 maxdV = subVectors3D(maxD.V2, maxD.V1);
  Vector3 maxDn = normalizeVector3D(maxdV);
  Vector3 section_P1 = {maxdV.x * 0.2, maxdV.y * 0.2, maxdV.z * 0.2};
  Vector3 section_P1a = {maxdV.x * 0.3, maxdV.y * 0.3, maxdV.z * 0.3};
  Vector3 section_P2 = {maxdV.x * 0.8, maxdV.y * 0.8, maxdV.z * 0.8};
  Vector3 section_P2a = {maxdV.x * 0.7, maxdV.y * 0.7, maxdV.z * 0.7};
  section_P1 = addVectors3D(maxD.V1, section_P1);
  section_P1a = addVectors3D(maxD.V1, section_P1a);
  section_P2 = addVectors3D(maxD.V1, section_P2);
  section_P2a = addVectors3D(maxD.V1, section_P2a);
  // calculate the sectioning points for each plane
  vector<Vector3> cs_plane_1 = sliceMesh(Mesh3D, section_P1, maxDn);
  vector<Vector3> cs_plane_1a = sliceMesh(Mesh3D, section_P1a, maxDn);
  vector<Vector3> cs_plane_2 = sliceMesh(Mesh3D, section_P2, maxDn);
  vector<Vector3> cs_plane_2a = sliceMesh(Mesh3D, section_P2a, maxDn);
  // calculate each cross-section's centroid
  Vector3 cs_centroid_1 = centroidPolygon3D(cs_plane_1, maxDn);
  Vector3 cs_centroid_1a = centroidPolygon3D(cs_plane_1a, maxDn);
  Vector3 cs_centroid_2 = centroidPolygon3D(cs_plane_2, maxDn);
  Vector3 cs_centroid_2a = centroidPolygon3D(cs_plane_2a, maxDn);
  // calculate normals between adjacent cross-section centroids
  // so that normals point towards the nearest epiphysis
  Vector3 cscdV1a_1 = subVectors3D(cs_centroid_1, cs_centroid_1a);
  Vector3 cscDn1a_1 = normalizeVector3D(cscdV1a_1);
  Vector3 cscdV2a_2 = subVectors3D(cs_centroid_2, cs_centroid_2a);
  Vector3 cscDn2a_2 = normalizeVector3D(cscdV2a_2);
  // recalculate cross sections with the updated normal
  cs_plane_1 = sliceMesh(Mesh3D, section_P1, cscDn1a_1);
  cs_plane_1a = sliceMesh(Mesh3D, section_P1a, cscDn1a_1);
  cs_plane_2 = sliceMesh(Mesh3D, section_P2, cscDn2a_2);
  cs_plane_2a = sliceMesh(Mesh3D, section_P2a, cscDn2a_2);
  // recalculate each cross-section's centroid
  cs_centroid_1 = centroidPolygon3D(cs_plane_1, cscDn1a_1);
  cs_centroid_1a = centroidPolygon3D(cs_plane_1a, cscDn1a_1);
  cs_centroid_2 = centroidPolygon3D(cs_plane_2, cscDn2a_2);
  cs_centroid_2a = centroidPolygon3D(cs_plane_2a, cscDn2a_2);
  // recalculate normals between adjacent cross-section centroids
  cscdV1a_1 = subVectors3D(cs_centroid_1, cs_centroid_1a);
  cscDn1a_1 = normalizeVector3D(cscdV1a_1);
  cscdV2a_2 = subVectors3D(cs_centroid_2, cs_centroid_2a);
  cscDn2a_2 = normalizeVector3D(cscdV2a_2);
  // translate 20% cross-sectional centroid to origin
  vector<Vector3> Cloud3D_1 = translatePoints3D(Cloud3D, cs_centroid_1);
  // find angle of rotation between 20% cross-sectional normal and x axis
  Vector3 x_axis = {1, 0, 0};
  Vector3 cscDn1_x = crossProduct3D(cscDn1a_1, x_axis);
  Vector3 rotAxis1 = normalizeVector3D(cscDn1_x);
  double theta_nom = dotProduct3D(cscDn1_x, rotAxis1);
  double theta_den = dotProduct3D(cscDn1a_1, x_axis);
  double theta_1_x = atan2(theta_nom, theta_den);
  // rotate model so that 20% cross-sectional normal aligns with x axis
  // note: the nearest epiphysis lies on the positive side of x axis
  Cloud3D_1 = rotatePoints3D(Cloud3D_1, rotAxis1, theta_1_x);
  // map 2D projection as the Y-Z plane of 3D space
  Vector3 X2D_vector = {0, 1, 0};
  Vector3 Y2D_vector = {0, 0, 1};
  // project the epiphysis near the 20% centroid onto the Y-Z plane
  // keep points only on the positive side of x axis to be projected
  vector<Vector3> epCloud_1;
  for(int i = 0; i < Cloud3D_1.size(); i++)
  {
    if(Cloud3D_1[i].x > 0)
    {
      epCloud_1.push_back(Cloud3D_1[i]);
    }
  }
  vector<Vector2> ZYproj_1;
  for(int i = 0; i < epCloud_1.size(); i++)
  {
    double x2D = dotProduct3D(epCloud_1[i], X2D_vector);
    double y2D = dotProduct3D(epCloud_1[i], Y2D_vector);
    ZYproj_1.push_back({x2D, y2D});
  }
  // translate the 80% centroid point to the origin
  vector<Vector3> Cloud3D_2 = translatePoints3D(Cloud3D, cs_centroid_2);
  // find angle of rotation between 80% cross-sectional normal and x axis
  Vector3 cscDn2_x = crossProduct3D(cscDn2a_2, x_axis);
  Vector3 rotAxis2 = normalizeVector3D(cscDn2_x);
  theta_nom = dotProduct3D(cscDn2_x, rotAxis2);
  theta_den = dotProduct3D(cscDn2a_2, x_axis);
  double theta_2_x = atan2(theta_nom, theta_den);
  // rotate model so that 80% cross-sectional normal aligns with x axis
  // note: the nearest epiphysis lies on the positive side of x axis
  Cloud3D_2 = rotatePoints3D(Cloud3D_2, rotAxis2, theta_2_x);
  // project the epiphysis near the 20% centroid onto the Y-Z plane
  // keep points only on the positive side of x axis to be projected
  vector<Vector3> epCloud_2;
  for(int i = 0; i < Cloud3D_2.size(); i++)
  {
    if(Cloud3D_2[i].x > 0)
    {
      epCloud_2.push_back(Cloud3D_2[i]);
    }
  }
  vector<Vector2> ZYproj_2;
  for(int i = 0; i < epCloud_2.size(); i++)
  {
    double x2D = dotProduct3D(epCloud_2[i], X2D_vector);
    double y2D = dotProduct3D(epCloud_2[i], Y2D_vector);
    ZYproj_2.push_back({x2D, y2D});
  }
  // process 2D projections to get the respective polylines and normalized EFDs 
  vector<vector<Raster>> rasterMatrix1 = rasterizeCloud2D(ZYproj_1);
  vector<Contour> boundaryPixels1 = calculateMooreNeighbor(rasterMatrix1);
  vector<Vector2> polyline2D_1 = extract2Dpolyline(boundaryPixels1, ZYproj_1);
  vector<FCCode> FreemanCCode1 = extractFreemanChainCode(boundaryPixels1);
  vector<EFDcoef> EFDs_1 = calculaleEFD(FreemanCCode1);
  vector<EFDcoef> nEFD_1 = normalizeEFD(EFDs_1);
  vector<vector<Raster>> rasterMatrix2 = rasterizeCloud2D(ZYproj_2);
  vector<Contour> boundaryPixels2 = calculateMooreNeighbor(rasterMatrix2);
  vector<Vector2> polyline2D_2 = extract2Dpolyline(boundaryPixels2, ZYproj_2);
  vector<FCCode> FreemanCCode2 = extractFreemanChainCode(boundaryPixels2);
  vector<EFDcoef> EFDs_2 = calculaleEFD(FreemanCCode2);
  vector<EFDcoef> nEFD_2 = normalizeEFD(EFDs_2);
  // get geometric propetries of projected polylines
  Polygon GEOM_1 = simplePolygon2D(polyline2D_1);
  Polygon GEOM_2 = simplePolygon2D(polyline2D_2);
  // classify long bone and proximal distal epiphysis
  classifyBone(nEFD_1, nEFD_2, GEOM_1, GEOM_2);
  // if global::epiphysis == 1 then distal Femur, proximal Tibia, distal
  // Humerus, and proximal Ulna epiphyses are stored in designated _1 variables
  //
  // check for proximal-distal epiphysis, register points in 2D and 3D space on
  // the correct epiphysis and apply the appropriate transformation so that
  // registration points correspond to the original model orientation.
  // 1) rotate around rotAxis1 by -theta_1_x so that x axis aligns with cscDn1a_1
  //    and translate the origin at 20% centroid point, i.e. cs_centroid_1
  // 2) rotate around rotAxis2 by -theta_2_x so that x axis aligns with cscDn2a_2
  //    and translate the origin at 80% centroid point, i.e. cs_centroid_2
  vector<Vector2> regPoints2D;
  vector<Vector3> regPoints3D;
  if(global::epiphysis == 1)
  {
    regPoints2D = register2Dpolyline(polyline2D_1, GEOM_1);
    regPoints3D = register3DepCloud(epCloud_1, regPoints2D);
    regPoints3D = rotatePoints3D(regPoints3D, rotAxis1, -theta_1_x);
    cs_centroid_1 = {-cs_centroid_1.x, -cs_centroid_1.y, -cs_centroid_1.z};
    regPoints3D = translatePoints3D(regPoints3D, cs_centroid_1);
  }
  else if(global::epiphysis == 2)
  {
    regPoints2D = register2Dpolyline(polyline2D_2, GEOM_2);
    regPoints3D = register3DepCloud(epCloud_2, regPoints2D);
    regPoints3D = rotatePoints3D(regPoints3D, rotAxis2, -theta_2_x);
    cs_centroid_2 = {-cs_centroid_2.x, -cs_centroid_2.y, -cs_centroid_2.z};
    regPoints3D = translatePoints3D(regPoints3D, cs_centroid_2);
  }
  else
  {
    Vector3 tmpPoint = {0, 0, 0};
    regPoints3D.push_back(tmpPoint);
    regPoints3D.push_back(tmpPoint);
  }
  // return registered 3D points
  Matrix regPoints(2, 3);
  regPoints(0,0) = regPoints3D[0].x;
  regPoints(0,1) = regPoints3D[0].y;
  regPoints(0,2) = regPoints3D[0].z;
  regPoints(1,0) = regPoints3D[1].x;
  regPoints(1,1) = regPoints3D[1].y;
  regPoints(1,2) = regPoints3D[1].z;
  // define return value list
  octave_value_list retval;
  retval(0) = global::bone;
  if(nargout == 2)
  {
    retval(1) = regPoints;
  }
  return octave_value(retval);
}
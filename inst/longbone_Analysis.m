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
## @deftypefn{Script} longbone_Analysis
##
## This script reads the available bone models stored in OBJ file format that
## are present in some folder and utilizes the 'longbone_Geometry' function to
## analyze their geometric properties. The user is prompted to select which
## long bone(s) should be analyzed and the folder which contains the OBJ files.
## All available OBJ files are evaluated but only the explicitly selected bones
## are fully processed and their corresponding geometric properties stored in the
## current working directory as .csv files named following the name convention
## of their initial mesh filename.
##
## For example, for a triangular mesh named as 'bone_ID.obj', the following
## files are produced:
##
## @multitable {Filename} {Description} @columnfractions .35 .65
## @item geometry-bone_ID.csv @tab containing the properties of area, perimeter,
## centroid, slicing and orientation normals for each cross section stored as a
## row vector in the order (:,[1],[2],[3:5],[6:8],[9:11]) from proximal to distal
## cross sections ([1:5],:).
##
## @item inertia-bone_ID.csv @tab containing the properties of Ix, Iy, Ixy,
## Imin, Imax and theta angle for each cross section as a row vector (:,[1:6])
## from proximal to distal cross sections ([1:5],:).
##
## @item polyline2D-bone_ID.csv @tab containing the 2D coordinates for each
## cross section as an Nx2 matrix, where N is the number of points for each
## polygon. The polygons are ordered as (:,[1:2],[3:4],[5:6],[7:8],[9:10]) from
## proximal to distal cross sections.
##
## @item polyline3D-bone_ID.csv @tab containing the 3D coordinates for each
## cross section as an Nx3 matrix, where N is the number of points for each
## polygon. The polygons are ordered as (:,[1:3],[4:6],[7:9],[10:12],[13:15])
## from proximal to distal cross sections.
## @end multitable
##
## For each selected OBJ file, the initial alignment points are either
## read from the corresponding Meshlab Point file (e.g. 'bone_ID.pp'), if
## present in the same folder with the OBJ, or they are automatically registered
## with the 'longbone_Registration' function and they are saved in newly created
## Meshlab Point files. Both sides of the each bone can be analyzed together,
## but each OBJ file must explicitly contain a single bone.
##
## All bone models stored as OBJ files must be strictly triangular meshes and
## their coordinate umits are assumed to be in mm. Only intact femur, tibia,
## humerus, and ulna bones can be processed.
## @seealso{inspect_CSG, visualize_CrossSections, longbone_Geometry}
## @end deftypefn

## define the bone(s) to be analyzed from the available 3D model files
options = {"Humerus", "Ulna", "Femur", "Tibia", "All"};
pstring = "Select desired bones for analysis";
[bones, OK] = listdlg ("ListString", options, "SelectionMode", "multiple", ...
                      "ListSize", [222, 150], "InitialValue", 4, "Name", ...
                      "longbone_Analysis", "PromptString", pstring);
## define folder containing the available .obj files
dialog = "Select folder containing 3D models for analysis";
folder = uigetdir (dialog);

## list the filenames with .obj extension in the working folder
filenames = dir (fullfile (folder, "*.obj"));
## calculate the geometric properties for each 
for i = 1:length (filenames)
  filename = strcat (filenames(i).name);
  [CS_Geometry, SMoA, polyline] = longbone_Geometry (folder, filename, bones);
  ## check for non empty structures
  if (isstruct(CS_Geometry) && isstruct(SMoA) && isstruct(polyline))
    
    geometry(:,1) = [CS_Geometry(1).Area; CS_Geometry(2).Area; ...
                     CS_Geometry(3).Area; CS_Geometry(4).Area; ...
                     CS_Geometry(5).Area];
    geometry(:,2) = [CS_Geometry(1).Perimeter; CS_Geometry(2).Perimeter; ...
                     CS_Geometry(3).Perimeter; CS_Geometry(4).Perimeter; ...
                     CS_Geometry(5).Perimeter];
    geometry(:,[3:5]) = [CS_Geometry(1).Centroid; CS_Geometry(2).Centroid; ...
                         CS_Geometry(3).Centroid; CS_Geometry(4).Centroid; ...
                         CS_Geometry(5).Centroid];
    geometry(:,[6:8]) = [CS_Geometry(1).Section_n; CS_Geometry(2).Section_n; ...
                         CS_Geometry(3).Section_n; CS_Geometry(4).Section_n; ...
                         CS_Geometry(5).Section_n];
    geometry(:,[9:11]) = [CS_Geometry(1).Coronal_n; CS_Geometry(2).Coronal_n;...
                          CS_Geometry(3).Coronal_n; CS_Geometry(4).Coronal_n;...
                          CS_Geometry(5).Coronal_n];

    inertia(1,:) = [SMoA(1).Ix, SMoA(1).Iy, SMoA(1).Ixy, ...
                    SMoA(1).Imin, SMoA(1).Imax, SMoA(1).theta];
    inertia(2,:) = [SMoA(2).Ix, SMoA(2).Iy, SMoA(2).Ixy, ...
                    SMoA(2).Imin, SMoA(2).Imax, SMoA(2).theta];
    inertia(3,:) = [SMoA(3).Ix, SMoA(3).Iy, SMoA(3).Ixy, ...
                    SMoA(3).Imin, SMoA(3).Imax, SMoA(3).theta];
    inertia(4,:) = [SMoA(4).Ix, SMoA(4).Iy, SMoA(4).Ixy, ...
                    SMoA(4).Imin, SMoA(4).Imax, SMoA(4).theta];
    inertia(5,:) = [SMoA(5).Ix, SMoA(5).Iy, SMoA(5).Ixy, ...
                    SMoA(5).Imin, SMoA(5).Imax, SMoA(5).theta];

    polygon2D([1:length(polyline(1).poly2D)],[1:2]) = polyline(1).poly2D;
    polygon2D([1:length(polyline(2).poly2D)],[3:4]) = polyline(2).poly2D;
    polygon2D([1:length(polyline(3).poly2D)],[5:6]) = polyline(3).poly2D;
    polygon2D([1:length(polyline(4).poly2D)],[7:8]) = polyline(4).poly2D;
    polygon2D([1:length(polyline(5).poly2D)],[9:10]) = polyline(5).poly2D;

    polygon3D([1:length(polyline(1).poly3D)],[1:3]) = polyline(1).poly3D;
    polygon3D([1:length(polyline(2).poly3D)],[4:6]) = polyline(2).poly3D;
    polygon3D([1:length(polyline(3).poly3D)],[7:9]) = polyline(3).poly3D;
    polygon3D([1:length(polyline(4).poly3D)],[10:12]) = polyline(4).poly3D;
    polygon3D([1:length(polyline(5).poly3D)],[13:15]) = polyline(5).poly3D;

    starting = "geometry-";
    name = filename([1:length(filename)-4]);
    endfile = ".csv";
    filename = strcat (starting, name, endfile);
    csvwrite (filename, geometry);
    starting = "inertia-";
    filename = strcat( starting, name, endfile);
    csvwrite (filename, inertia);
    starting = "polyline2D-";
    filename = strcat (starting, name, endfile);
    csvwrite (filename, polygon2D);
    starting = "polyline3D-";
    filename = strcat (starting, name, endfile);
    csvwrite (filename, polygon3D);
  endif
  clear geometry; clear inertia; clear polygon2D; clear polygon3D; 
  clear CS_Geometry; clear SMoA; clear polyline;
endfor
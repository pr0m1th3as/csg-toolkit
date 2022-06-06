## Copyright (C) 2020-2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn{Script} inspect_CSG
##
## This script reads all the available "geometry-$$.csv" files present in the
## working folder to establish the number of different bones present and utilizes
## the 'visualize_CrossSections' function to plot the data for visual inspection
## for each individual bone's cross-sectional contours, while asking the user for
## verification. If contours are shown as a simple polygon, then the results are
## valid and the user may confirm by pressing the "yes" button in order to read
## and append the CSG properties into tabular form. If not, the particular sample
## is ignored and the inspection proceeds to the next available sample. If the
## user cancels, the process is terminated and the already appended values are
## saved in csv file named "CSG Data.csv".
##
## For every sample, the script checks for the presence of the associated 3D
## model as an OBJ file and if it is available it measures the bone's maximum
## length and includes it in the data, otherwise a "NaN" value is appended. The
## OBJ files may be in a different folder/path than the current working directory.
## @seealso{visualize_CrossSections, longbone_Analysis}
## @end deftypefn

## list the filenames with geometry-*.csv pattern in the working folder
filenames = dir ("geometry-*.csv");
## define folder containing the .obj files
dialog = "Select folder containing 3D models for longbone analysis";
folder = uigetdir (dialog);

## create the header of the final csv file in a cell array
complete(1,1) = {"sample"}; complete(1,2) = {"Max Distance"};
complete(1,3) = {"angle-20_35"}; complete(1,4) = {"angle-35_50"};
complete(1,5) = {"angle-50_65"}; complete(1,6)={"angle-65_80"};
complete(1,7) = {"diaphyseal bending"};
complete(1,8) = {"Area 20%"}; complete(1,9) = {"Perimeter 20%"};
complete(1,10) = {"Ix 20%"};complete(1,11) = {"Iy 20%"};
complete(1,12) = {"Ixy 20%"}; complete(1,13) = {"Imin 20%"};
complete(1,14) = {"Imax 20%"}; complete(1,15) = {"Theta 20%"};
complete(1,16) = {"Area 35%"}; complete(1,17) = {"Perimeter 35%"};
complete(1,18) = {"Ix 35%"}; complete(1,19) = {"Iy 35%"};
complete(1,20) = {"Ixy 35%"}; complete(1,21) = {"Imin 35%"};
complete(1,22) = {"Imax 35%"}; complete(1,23) = {"Theta 35%"};
complete(1,24) = {"Area 50%"}; complete(1,25) = {"Perimeter 50%"};
complete(1,26) = {"Ix 50%"}; complete(1,27) = {"Iy 50%"};
complete(1,28) = {"Ixy 50%"}; complete(1,29) = {"Imin 50%"};
complete(1,30) = {"Imax 50%"}; complete(1,31) = {"Theta 50%"};
complete(1,32) = {"Area 65%"}; complete(1,33) = {"Perimeter 65%"};
complete(1,34) = {"Ix 65%"}; complete(1,35) = {"Iy 65%"};
complete(1,36) = {"Ixy 65%"}; complete(1,37) = {"Imin 65%"};
complete(1,38) = {"Imax 65%"}; complete(1,39) = {"Theta 65%"};
complete(1,40) = {"Area 80%"}; complete(1,41) = {"Perimeter 80%"};
complete(1,42) = {"Ix 80%"}; complete(1,43) = {"Iy 80%"};
complete(1,44) = {"Ixy 80%"}; complete(1,45) = {"Imin 80%"};
complete(1,46) = {"Imax 80%"}; complete(1,47) = {"Theta 80%"};

id = 0;
## for each geometry-*.csv files
for i = 1:length (filenames)
  ## extract the OBJ model's name and call the visualization function
  bone_id = strcat (filenames(i).name(10:end-4));
  [geometry, inertia] = visualize_CrossSections (bone_id);
  btn = questdlg ("Are CSG properties shown properly?", bone_id, "Yes", "No");
  if (strcmp (btn, "Yes"))
    close all;
    id += 1;
    ## read full data from relevant geometry csv files
    g_filename = strcat ("geometry-", bone_id, ".csv");
    full_geometry = csvread (g_filename);
    ## write sample's id
    complete(i+1,1) = {bone_id};
    ## check if 3D model's OBJ file is present in the working directory and in
    ## such case calculate its max distance and append it in the cell array,
    ## otherwise append a NaN value
    obj_filename = strcat (bone_id, ".obj");
    if (exist (fullfile (folder, obj_filename)) == 2)
      [v,f] = readObj (fullfile (folder, obj_filename));
      ## keep only vertex coordinates (just in case color information is present)
      v = v(:,[1:3]);
      maxD = longbone_maxDistance (v);
      complete(i+1,2) = {maxD};
    else
      complete(i+1,2) = {NaN};
    endif
    ## compute dihedral angles between cross sectional plane vectors
    a1_2 = rad2deg (acos (dot (full_geometry(1,[6:8]), full_geometry(2,[6:8]))));
    a2_3 = rad2deg (acos (dot (full_geometry(2,[6:8]), full_geometry(3,[6:8]))));
    a3_4 = rad2deg (acos (dot (full_geometry(3,[6:8]), full_geometry(4,[6:8]))));
    a4_5 = rad2deg (acos (dot (full_geometry(4,[6:8]), full_geometry(5,[6:8]))));
    asum = a1_2 + a2_3 + a3_4 + a4_5;
    ## store in a cell array
    complete(i+1,[3:7])={a1_2, a2_3, a3_4, a4_5, asum};
    complete(i+1,[8:47]) = {geometry(1).Area, geometry(1).Perimeter, ...
                          inertia(1).Ix, inertia(1).Iy, inertia(1).Ixy, ...
                          inertia(1).Imin, inertia(1).Imax, inertia(1).theta, ...
                          geometry(2).Area, geometry(2).Perimeter, ...
                          inertia(2).Ix, inertia(2).Iy, inertia(2).Ixy, ...
                          inertia(2).Imin, inertia(2).Imax, inertia(2).theta, ...
                          geometry(3).Area, geometry(3).Perimeter, ...
                          inertia(3).Ix, inertia(3).Iy, inertia(3).Ixy, ...
                          inertia(3).Imin, inertia(3).Imax, inertia(3).theta, ...
                          geometry(4).Area, geometry(4).Perimeter, ...
                          inertia(4).Ix, inertia(4).Iy, inertia(4).Ixy, ...
                          inertia(4).Imin, inertia(4).Imax, inertia(4).theta, ...
                          geometry(5).Area, geometry(5).Perimeter, ...
                          inertia(5).Ix, inertia(5).Iy, inertia(5).Ixy, ...
                          inertia(5).Imin, inertia(5).Imax, inertia(5).theta};
  elseif (strcmp (btn, "No"))
    close all;
  else
    close all;
    break;
  endif
endfor
cell2csv ("CSG Data.csv", complete);
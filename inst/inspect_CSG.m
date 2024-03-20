## Copyright (C) 2020-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} inspect_CSG
## @deftypefnx {csg-toolkit} {} inspect_CSG (@qcode{"default"})
## @deftypefnx {csg-toolkit} {} inspect_CSG (@qcode{"custom"})
## @deftypefnx {csg-toolkit} {} inspect_CSG (@qcode{"fragment"})
## @deftypefnx {csg-toolkit} {} inspect_CSG (@dots{}, @var{filename})
## @deftypefnx {csg-toolkit} {} inspect_CSG (@var{filename})
##
## Batch inspect the CSG properties of previously analyzed long bones.
##
## @code{inspect_CSG} reads all available @qcode{$geometry-$$.csv} files in the
## working folder, utilizes @code{visualize_CrossSections} to retrieve the data
## in the corresponding CSV files for each associated 3D model and plot the
## cross sectional contours, and upon user confirmation, it aggregates the
## CSG properties in a tabular format and save it under @var{filename}.
## If @var{filename} is not provided, @qcode{"CSG Data.csv"} is used by default.
##
## The rationale behind the visual inspection of the cross sectional contours is
## to avoid any erroneous values due to badly shaped long bone models that might
## be processed succesfully by the respective @qcode{Geometry} function.  For
## any given sample, if the contours of are shown as a simple polygon, then the
## computed CSG properties can be safely assumed as accurate and the user can
## confirm by pressing the "yes" button in order to append the CSG properties
## into tabular form. If not, press "no" to ignore the particular sample and
## proceed the inspection to the next available sample. If the user cancels, the
## process is terminated and the already appended values are saved in the CSV
## file.
##
## Prior to batch processing, the user is prompted for the folder containing the
## 3D models of the associated bones and, if available, @code{inspect_CSG}
## measures their maximum length (max Distance) and includes it in the
## aggregated results, otherwise it adds NaN, as a missing value in the stored
## CSV file.
##
## By default, @code{inspect_CSG} processes the available @qcode{Dgeometry-},
## @qcode{Dinertia-}, @qcode{Dpolyline2D-}, and @qcode{Dpolyline3D-} CSV files
## generated with the @code{longbone_Geometry} function.  This is identical to
## calling @code{inspect_CSG} with the optional argument @qcode{"default"}.  In
## this batch prossecing mode, 5 plots are rendered, one for each cross section
## at 20%, 35%, 50%, 65%, and 80%, respectively.
##
## @code{inspect_CSG (@qcode{"custom"})} will process all available
## @qcode{Cgeometry-}, @qcode{Cinertia-}, @qcode{Cpolyline2D-}, and
## @qcode{Cpolyline3D-} CSV files generated with the
## @code{longbone_CustomGeometry} function.  In this batch prossecing mode, an
## arbitrary number of plots is rendered, one for each associated cross section
## according to the contents of the respective CSV files.  Note that all files
## included in the same batch MUST contain the same number of cross sections.
## However, it is not mandatory that they conform to the same points along the
## diaphysis.  This is especially important for 3D modes analyzed with custom PP
## point files, where the ratio used for indexing each cross section might vary
## between samples.  The ratio IDs used in the resulting @qcode{"CSG Data.csv"}
## are taken from the ultimate processed sample.  Also note that the first and
## last plots always refer to the cross sections at 20% and 80%, respectively.
##
## @code{inspect_CSG (@qcode{"fragment"})} will process all available
## @qcode{Fgeometry-}, @qcode{Finertia-}, @qcode{Fpolyline2D-}, and
## @qcode{Fpolyline3D-} CSV files generated with the
## @code{longbone_FragmentGeometry} function.  In this batch prossecing mode, an
## arbitrary number of plots is rendered, one for each associated cross section
## according to the contents of the respective CSV files.  Note that all files
## included in the same batch MUST contain the same number of cross sections.
## However, it is not mandatory that they conform to the same points along the
## diaphysis.
##
## @seealso{visualize_CrossSections, longbone_Analysis}
## @end deftypefn

function inspect_CSG (varargin)

  ## Check input
  if (nargin > 2)
    print_usage;
  endif

  ## Check for default or custom geometry according to input argument (if any)
  if (nargin == 0)      # default slicing
    type = "default";
    fout = "CSG Data.csv";
  elseif (nargin == 1)
    if (isempty (strfind (varargin{1}, ".csv")))
      type = varargin{1};
      fout = "CSG Data.csv";
    else
      type = "default";
      fout = varargin{1};
    endif
  else
    type = varargin{1};
    fout = varargin{2};
  endif

  ## List the filenames with $geometry-*.csv pattern in the working folder
  switch (lower (type))
    case "default"
      filenames = dir ("Dgeometry-*.csv");
    case "custom"
      filenames = dir ("Cgeometry-*.csv");
    case "fragment"
      filenames = dir ("Fgeometry-*.csv");
    otherwise
      print_usage ();
  endswitch

  ## Get user selection for the folder containing the 3D models
  dialog = "Select folder containing analyzed 3D models";
  folder = uigetdir (dialog);

  ## Process $geometry-*.csv files
  switch (lower (type))

    case "default"
      ## Add header in the cell array for the CSV file
      complete(1,1) = {"sample"};            complete(1,2) = {"Max Distance"};
      complete(1,3) = {"angle-20_35"};       complete(1,4) = {"angle-35_50"};
      complete(1,5) = {"angle-50_65"};       complete(1,6) = {"angle-65_80"};
      complete(1,7) = {"diaphyseal bending"};
      complete(1,8) = {"Area 20%"};          complete(1,9) = {"Perimeter 20%"};
      complete(1,10) = {"Ix 20%"};           complete(1,11) = {"Iy 20%"};
      complete(1,12) = {"Ixy 20%"};          complete(1,13) = {"Imin 20%"};
      complete(1,14) = {"Imax 20%"};         complete(1,15) = {"Theta 20%"};
      complete(1,16) = {"Area 35%"};         complete(1,17) = {"Perimeter 35%"};
      complete(1,18) = {"Ix 35%"};           complete(1,19) = {"Iy 35%"};
      complete(1,20) = {"Ixy 35%"};          complete(1,21) = {"Imin 35%"};
      complete(1,22) = {"Imax 35%"};         complete(1,23) = {"Theta 35%"};
      complete(1,24) = {"Area 50%"};         complete(1,25) = {"Perimeter 50%"};
      complete(1,26) = {"Ix 50%"};           complete(1,27) = {"Iy 50%"};
      complete(1,28) = {"Ixy 50%"};          complete(1,29) = {"Imin 50%"};
      complete(1,30) = {"Imax 50%"};         complete(1,31) = {"Theta 50%"};
      complete(1,32) = {"Area 65%"};         complete(1,33) = {"Perimeter 65%"};
      complete(1,34) = {"Ix 65%"};           complete(1,35) = {"Iy 65%"};
      complete(1,36) = {"Ixy 65%"};          complete(1,37) = {"Imin 65%"};
      complete(1,38) = {"Imax 65%"};         complete(1,39) = {"Theta 65%"};
      complete(1,40) = {"Area 80%"};         complete(1,41) = {"Perimeter 80%"};
      complete(1,42) = {"Ix 80%"};           complete(1,43) = {"Iy 80%"};
      complete(1,44) = {"Ixy 80%"};          complete(1,45) = {"Imin 80%"};
      complete(1,46) = {"Imax 80%"};         complete(1,47) = {"Theta 80%"};

      ## Process samples
      id = 0;
      for i = 1:length (filenames)
        ## Extract the OBJ model's name and call the visualization function
        bone_id = strcat (filenames(i).name(11:end-4));
        [geometry, inertia] = visualize_CrossSections (bone_id);

        ## Get user confirmation
        btn = questdlg ("Are cross-sections shown properly?", ...
                        bone_id, "Yes", "No");
        if (strcmp (btn, "Yes"))
          close all;
          id += 1;
          ## Write sample's id
          complete(i+1,1) = {bone_id};

          ## Check if 3D model's OBJ file is present in the user selected
          ## directory and calculate its max distance, otherwise force NaN
          obj_filename = strcat (bone_id, ".obj");
          if (exist (fullfile (folder, obj_filename)) == 2)
            [v,f] = readObj (fullfile (folder, obj_filename));
            ## Keep only vertex coordinates
            v = v(:,[1:3]);
            maxD = longbone_maxDistance (v);
            complete(i+1,2) = {maxD};
          else
            complete(i+1,2) = {NaN};
          endif

          ## Compute dihedral angles between cross sectional plane normals
          a1_2 = rad2deg (acos (dot (geometry(1).Section_n, ...
                                     geometry(2).Section_n)));
          a2_3 = rad2deg (acos (dot (geometry(2).Section_n, ...
                                     geometry(3).Section_n)));
          a3_4 = rad2deg (acos (dot (geometry(3).Section_n, ...
                                     geometry(4).Section_n)));
          a4_5 = rad2deg (acos (dot (geometry(4).Section_n, ...
                                     geometry(5).Section_n)));
          asum = a1_2 + a2_3 + a3_4 + a4_5;

          ## Store angles to cell array
          complete(i+1,[3:7])={a1_2, a2_3, a3_4, a4_5, asum};

          ## Store CSG properties
          complete(i+1,[8:47]) = ...
                  {geometry(1).Area, geometry(1).Perimeter, ...
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
      cell2csv (fout, complete);

    case "custom"
      ## Process samples
      id = 0;
      for i = 1:length (filenames)
        ## Extract the OBJ model's name and call the visualization function
        bone_id = strcat (filenames(i).name(11:end-4));
        [geometry, inertia] = visualize_CrossSections (bone_id, type);

        ## Get user confirmation
        btn = questdlg ("Are cross-sections shown properly?", ...
                        bone_id, "Yes", "No");
        if (strcmp (btn, "Yes"))
          close all;
          id += 1;
          ## Write sample's id
          complete(i+1,1) = {bone_id};

          ## Check if 3D model's OBJ file is present in the user selected
          ## directory and calculate its max distance, otherwise force NaN
          obj_filename = strcat (bone_id, ".obj");
          if (exist (fullfile (folder, obj_filename)) == 2)
            [v,f] = readObj (fullfile (folder, obj_filename));
            ## Keep only vertex coordinates
            v = v(:,[1:3]);
            maxD = longbone_maxDistance (v);
            complete(i+1,2) = {maxD};
          else
            complete(i+1,2) = {NaN};
          endif

          ## Store CSG properties
          for cs = 1:numel (geometry)
            complete(i+1,[(cs-1)*8+3:(cs-1)*8+10]) = ...
                    {geometry(cs).Area, geometry(cs).Perimeter, ...
                     inertia(cs).Ix, inertia(cs).Iy, inertia(cs).Ixy, ...
                     inertia(cs).Imin, inertia(cs).Imax, inertia(cs).theta};
          endfor

        elseif (strcmp (btn, "No"))
          close all;
        else
          close all;
          break;
        endif
      endfor

      ## Add header in the cell array for the CSV file
      complete(1,1) = {"sample"};
      complete(1,2) = {"Max Distance"};
      for cs = 1:numel (geometry)
        csection = geometry(cs).CS;
        s1 = sprintf ("Area %0.2f%%", csection);
        s2 = sprintf ("Perimeter %0.2f%%", csection);
        s3 = sprintf ("Ix %0.2f%%", csection);
        s4 = sprintf ("Iy %0.2f%%", csection);
        s5 = sprintf ("Ixy %0.2f%%", csection);
        s6 = sprintf ("Imin %0.2f%%", csection);
        s7 = sprintf ("Imax %0.2f%%", csection);
        s8 = sprintf ("Theta %0.2f%%", csection);
        complete(1,[(cs-1)*8+3:(cs-1)*8+10]) = {s1, s2, s3, s4, s5, s6, s7, s8};
      endfor
      cell2csv (fout, complete);

    case "fragment"
      ## Process samples
      id = 0;
      for i = 1:length (filenames)
        ## Extract the OBJ model's name and call the visualization function
        bone_id = strcat (filenames(i).name(11:end-4));
        [geometry, inertia] = visualize_CrossSections (bone_id, type);
        ## Get user confirmation
        btn = questdlg ("Are cross-sections shown properly?", ...
                        bone_id, "Yes", "No");
        if (strcmp (btn, "Yes"))
          close all;
          id += 1;
          ## Write sample's id
          complete(i+1,1) = {bone_id};
          ## Check if 3D model's OBJ file is present in the user selected
          ## directory and calculate its max distance, otherwise force NaN
          obj_filename = strcat (bone_id, ".obj");
          if (exist (fullfile (folder, obj_filename)) == 2)
            [v,f] = readObj (fullfile (folder, obj_filename));
            ## Keep only vertex coordinates
            v = v(:,[1:3]);
            maxD = longbone_maxDistance (v);
            complete(i+1,2) = {maxD};
          else
            complete(i+1,2) = {NaN};
          endif

          ## Store CSG properties
          for cs = 1:numel (geometry)
            complete(i+1,[(cs-1)*5+3:(cs-1)*5+7]) = ...
                    {geometry(cs).Area, geometry(cs).Perimeter, ...
                     inertia(cs).Ixy, inertia(cs).Imin, inertia(cs).Imax};
          endfor
        elseif (strcmp (btn, "No"))
          close all;
        else
          close all;
          break;
        endif
      endfor
      ## Add header in the cell array for the CSV file
      complete(1,1) = {"sample"};
      complete(1,2) = {"Max Distance"};
      for cs = 1:numel (geometry)
        csection = geometry(cs).CS;
        s1 = sprintf ("Area %d%%", csection);
        s2 = sprintf ("Perimeter %d%%", csection);
        s3 = sprintf ("Ixy %d%%", csection);
        s4 = sprintf ("Imin %d%%", csection);
        s5 = sprintf ("Imax %d%%", csection);
        complete(1,[(cs-1)*5+3:(cs-1)*5+7]) = {s1, s2, s3, s4, s5};
      endfor
      cell2csv (fout, complete);

  endswitch

endfunction

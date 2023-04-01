## Copyright (C) 2018-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} longbone_Analysis
## @deftypefnx {csg-toolkit} {} longbone_Analysis (@qcode{"custom"})
## @deftypefnx {csg-toolkit} {} longbone_Analysis (@var{points})
## @deftypefnx {csg-toolkit} {} longbone_Analysis (@qcode{"fragment"})
##
## A wrapper for batch processing with @qcode{longbone_*Geometry} functions.
##
## @code{longbone_Analysis} reads the available 3D models stored in Wavefront
## OBJ file format present in some folder and calls the @code{longbone_Geometry}
## function to analyze their geometric properties in a serial manner.  The user
## is prompted for the directory to read from and subsequently the type of bones
## that should be analyzed.  All available 3D bone models are evaluated but only
## those that explicitly match the user selection are fully processed and their
## corresponding geometric properties stored in the corresponding CSV files in
## the current working directory.  The generated CSV files follow the naming
## conventions used by the @code{longbone_Geometry} function.
##
## @code{longbone_Analysis (@qcode{"custom"})} will accordingly utilize the
## @code{longbone_CustomGeometry} function and, besides the folder containing
## the 3D models and the user selected types of bones for processing, it will
## also pass the respective Meshlab PickedPoints filename with the custom
## sectioning points.  These files must reside in the same folder with the 3D
## models and their base filename must be appended with @qcode{"-custom"} with
## respect to the OBJ filename.  For example, assuming a 3D model named
## @qcode{"bone_ID.obj"}, the corresponding custom sectioning points file must
## be named @qcode{"bone_ID-custom.pp"}.
##
## @code{longbone_Analysis (@var{points})} will optionally utilize the
## @code{longbone_CustomGeometry} by parsing a preset numerical vector for
## custom sectioning points as ratios along the long bone's maximum distance.
## Similarly to using the @qcode{"custom"} optional argument, the generated CSV
## files will follow the naming conventions used by the
## @code{longbone_CustomGeometry} function.
##
## @code{longbone_Analysis (@qcode{"fragment"})} will utilize the
## @code{longbone_FragmentGeometry} function, which only requires the folder
## containing the 3D models, but ignores their bone type since it is not
## required for the relevant CSG computations.  The user is only prompted for
## this folder, which must also contain the corresponding Meshlab PickedPoints
## files whose base filename must be appended with @qcode{"-fragment"} with
## respect to the OBJ filename.  For example, assuming a 3D model named
## @qcode{"bone_ID.obj"}, the corresponding custom sectioning points file must
## be named @qcode{"bone_ID-fragment.pp"}.
##
## Under all batch processing schemes except for the @qcode{"fragment"} option,
## the initial alignment points for each 3D model are either read from the
## corresponding Meshlab PickedPoints file (e.g. @qcode{"bone_ID.pp"}), if
## present in the same folder with the OBJ file, or they are automatically
## registered with the @code{longbone_Registration} function.  In either case,
## the opitimized alignment points are saved in the corresponding Meshlab
## PickedPoints file with the same base filename.  For example, assuming a 3D
## model named @qcode{"bone_ID.obj"}, the corresponding initial alignment points
## file must be named @qcode{"bone_ID.pp"}.  Note that existing files will be
## overwritten, hence apart from the required initial alignment points,
## other points will be lost.  Type @code{help longbone_Geometry} for more
## information about initial alignment points.
##
## @seealso{inspect_CSG, longbone_CustomGeometry, longbone_FragmentGeometry,
## longbone_Geometry}
## @end deftypefn

function longbone_Analysis (varargin)

  ## Check input
  if (nargin > 1)
    print_usage;
  endif

  ## Check for default or custom geometry according to input argument (if any)
  if (nargin == 0)        # default slicing
    mode = "default";
  else                    # parse input argument
    argin = varargin{1};
    if (strcmpi (argin, "custom"))                  # custom points file
      mode = "custom";
    elseif (isnumeric (argin) && isvector (argin))  # custom sectioning vector
      mode = "points";
    elseif (strcmpi (argin, "fragment"))            # fragment
      mode = "fragment";
    else
      error ("longbone_Analysis: invalid input argument.");
    endif
  endif

  ## Get user selection for the folder containing the 3D models
  dialog = "Select folder containing 3D models for analysis";
  folder = uigetdir (dialog);

  ## List the filenames with .obj extension in the selected folder
  filenames = dir (fullfile (folder, "*.obj"));

  ## For options other than "fragment"
  if (! strcmpi (mode, "fragment"))

    ## Get user selection for the bone(s) to be analyzed
    options = {"Humerus", "Ulna", "Femur", "Tibia", "All"};
    pstring = "Select desired bones for analysis";
    [bones, OK] = listdlg ("ListString", options, ...
                           "SelectionMode", "multiple", ...
                           "ListSize", [222, 150], ...
                           "InitialValue", 4, ...
                           "Name", "longbone_Analysis", ...
                           "PromptString", pstring);

    ## Process bone selection
    if (isempty (bones))
      bone = "All";
    elseif (numel (bones) == 1)
      switch (bones)
        case 1
          bone = {"Humerus"};
        case 2
          bone = {"Ulna"};
        case 3
          bone = {"Femur"};
        case 4
          bone = {"Tibia"};
        case 5
          bone = {"All"};
       endswitch
    elseif (length (bones) > 1)
      i = 0;
      if (any (bones == 1))
        i++;
        bone(i) = {"Humerus"};
      endif
      if (any (bones == 2))
        i++;
        bone(i) = {"Ulna"};
      endif
      if (any (bones == 3))
        i++;
        bone(i) = {"Femur"};
      endif
      if (any (bones == 4))
        i++;
        bone(i) = {"Tibia"};
      endif
      if (i == 4 || any (bones == 5))
        bone = {"All"};
      endif
    endif

  endif


  switch (lower (mode))
    case "default"      # default
      for i = 1:length (filenames)
        filename = strcat (filenames(i).name);
        longbone_Geometry (folder, filename, bone);
      endfor

    case "custom"       # custom points file
      for i = 1:length (filenames)
        ## Parse Meshlab PP filename for each bone
        filename = strcat (filenames(i).name);
        MLPPfile = strcat (filename([1:end-4]), "-custom.pp");
        if (exist (fullfile (folder, MLPPfile)) == 2)
          longbone_CustomGeometry (folder, filename, bone, MLPPfile);
        else
          printf (strcat (["Model %s does not have an associated %s"], ...
                          [" custom points file\n"]), filename, bone);
        endif
      endfor

    case "points"       # custom sectioning vector
      for i = 1:length (filenames)
        filename = strcat (filenames(i).name);
        longbone_CustomGeometry (folder, filename, bone, argin);
      endfor

    case "fragment"     # fragment
      for i = 1:length (filenames)
        ## Parse Meshlab PP filename for each bone
        filename = strcat (filenames(i).name);
        MLPPfile = strcat (filename([1:end-4]), "-fragment.pp");
        if (exist (fullfile (folder, MLPPfile)) == 2)
          longbone_FragmentGeometry (folder, filename, bone, MLPPfile);
        else
          printf (strcat (["Model %s does not have an associated %s"], ...
                          [" fragment points file\n"]), filename, bone);
        endif
      endfor

  endswitch

endfunction

## Copyright (C) 2018-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} longbone_Scaling
## @deftypefnx {csg-toolkit} {} longbone_Scaling (@var{filename})
## @deftypefnx {csg-toolkit} {} longbone_Scaling (@var{filename}, @var{maxD})
## @deftypefnx {csg-toolkit} {@var{scale} =} longbone_Scaling
## @deftypefnx {csg-toolkit} {@var{scale} =} longbone_Scaling (@dots{})
##
## Scale triangular meshes stored in Wavefront OBJ format.
##
## @code{longbone_Scaling} scales a 3D model in @var{filename} to a specified
## maximum distance in @var{maxD} and returns the scaling factor, while it
## stores the coordinates of the points corresponding to the maximum distance in
## a Meshlab PickedPoints file.
##
## The function may be used for scaling a single mesh or batch scaling multiple
## meshes.  When called with a single input argument, the function will scale
## the particular mesh defined in @var{filename} so that the scaled mesh will
## have a user-defined maximum distance.  The function is primarily intended for
## scaling long bone triangular mesh models produced with 3D photogrammetry,
## hence the use of maximum distance as the target reference for scaling.
## Nevertheless, it can be used for scaling any triangular mesh of arbitrary
## shape according to its maximum distance, which can either be parsed as a
## second input argument @var{maxD} or when prompted by the function.
##
## When called with no input arguments, the function will scan the working
## directory for all .obj files and will scale them iteratively in batch mode.
## For each available OBJ, the user will be prompted for the actual maximum
## distance of each corresponding long bone and it will subsequently scale the
## model to the real world dimensions with units in mm.  Scaling is performed
## about the model's barycentric coordinates, which are translated to origin.
##
## Prior to scaling, when prompted for manual measurement value, the user may
## examine the model's maximum distance corresponding points by opening the
## triangular mesh in Meshlab and utilizing the corresponding @qcode{.pp} file,
## which is created automatically by the present function.  After scaling is
## performed for a given triangular mesh, the user may re-examine the scaled
## model by reloading the mesh and its respective .pp file, which is
## automatically updated with the new points of its maximum distance.
##
## @code{@var{scale} = longbone_Scaling (@dots{})} returns a cell array with
## each 3D model's filename and scale related measurements in the following
## order: filename, ratio, old max distance, new max distance.  If no output
## argument is declared, the function will prompt the user to specify a CSV file
## (new or existing, which will be overwritten) to save these measurements.
## The first row, either in the CSV file or the @var{scale} cell array contains
## the repsective column labels @qcode{"filename"}, @qcode{"ratio"},
## @qcode{"oldMaxD"}, and @qcode{"newMaxD"}.
##
## @seealso{longbone_maxDistance}
## @end deftypefn

function [varargout] = longbone_Scaling (varargin)

  ## Check for valid number of input variables
  if (nargin > 2)
    error ("longbone_Scaling: invalid number of input arguments.");
  endif

  ## Check if input argument is provided
  if (nargin == 1)
    ## Must be a char string
    if (! ischar (varargin{1}))
      error ("longbone_Scaling: 3D mesh filename must be a string.");
    endif
    filenames(1).name = varargin{1};
    ask_user = true;

  elseif (nargin == 2)
    ## First input argument must be a string
    if (! ischar (varargin{1}))
      error ("longbone_Scaling: 3D mesh filename must be a string.");
    endif
    ## Second input argument must be a scalar
    if (! isscalar (varargin{2}) || ! isnumeric (varargin{2}))
      error ("longbone_Scaling: maximum distance must be a scalar.");
    endif
    filenames(1).name = varargin{1};
    RmaxD = varargin{2};
    ask_user = false;

  else
    ## List the filenames with .obj extension in the working folder
    filenames = dir ("*.obj");
    ask_user = true;
  endif

  ## Initialize header for the cell array
  scale = {"filename", "ratio", "oldMaxD", "newMaxD"};

  ## Perform scaling for each mesh object present in the working directory
  for i = 1:length (filenames)
    ## Store filename of current mesh
    filenameOBJ = strcat (filenames(i).name);

    ## Read original obj file
    [v, f, vt, ft, vn, fn, filenameMTL] = readObj (filenameOBJ);

    ## Calculate maximum distance and corresponding points
    [maxD, p1, p2] = longbone_maxDistance (v(:,[1:3]));

    ## Save maxD points to Meshlab .pp file using the name convention of the
    ## original mesh obj file
    filenamePP = strcat (filenameOBJ([1:end-4]), ".pp");
    write_MeshlabPoints (filenamePP, filenameOBJ, [p1; p2]);

    ## Ask user for actual maximum distance (real size)
    if ask_user
      realMaxD = inputdlg ("maximum distance in mm", filenameOBJ, [1,20]);
      RmaxD = str2num (cell2mat (realMaxD));
    endif

    ## Calculate scaling ratio
    ratio = RmaxD / maxD;

    ## Find barycentric coordinates and translate mesh' origin
    origin = meshBarycenter (v(:,[1:3]), f);
    v(:,[1:3]) = v(:,[1:3]) - origin;

    ## Scale mesh
    v(:,[1:3]) = v(:,[1:3]) * ratio;

    ## Check if texture coordinates and material file exist in OBJ
    if (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL))
      ## Read original mtl file
      MTL = mtlread (filenameMTL);
      ## Find which material's name corresponds to OBJ name and keep it or
      ## alternatively search for texture altas in non-empty  'map_Kd' field
      for id = 1:length (MTL)
        if (strcmp (getfield (MTL(id), "newmtl"), filenameOBJ([1:end-4])))
          newmtl_id = id;
        endif
        if (isfield (MTL(id), "map_Kd") &&
            ! isempty (getfield (MTL(id), "map_Kd")))
          map_Kd_id = id;
        endif
      endfor
    endif

    ## Save scaled model according to the elements present in the original OBJ
    if (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL) &&
        !isempty (vn) && !isempty (fn))
      filenameOBJ = writeObj (v, f, vt, ft, vn, fn, filenameOBJ);

    elseif (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL) &&
        isempty (vn) && isempty (fn))
      filenameOBJ = writeObj (v, f, vt, ft, filenameOBJ);

    elseif (isempty (vt) && isempty (ft) && isempty (filenameMTL) &&
        !isempty (vn) && !isempty (fn))
      filenameOBJ = writeObj (v, f, vn, fn, filenameOBJ);

    else
      filenameOBJ = writeObj (v, f, filenameOBJ);
    endif

    ## Calculate scaled model's max distance and save its points in .pp file
    [maxDnew, p1, p2] = longbone_maxDistance (v);
    write_MeshlabPoints (filenamePP, filenameOBJ, [p1; p2]);

    ## Check if texture coordinates and material file exist in OBJ (again!!)
    if (! isempty (vt) && ! isempty (ft) && ! isempty (filenameMTL))
      if (exist ("newmtl_id"))
        MTL = MTL(newmtl_id);
        % rename material to the OBJ's updated name (just in case!!)
        MTL.newmtl = filenameOBJ([1:end-4]);
      else
        MTL = MTL(map_Kd_id);
        % rename material to the OBJ's updated name (just in case!!)
        MTL.newmtl = filenameOBJ([1:end-4]);
      endif

      ## Check for texture map image in 'map_Kd' field and if available, then
      ## compare its filename to OBJ filename and if they don't match make a new
      ## copy of the image with the updated OBJ filename. Also update map_Kd
      if (isfield (MTL, "map_Kd") && ! isempty (getfield (MTL, "map_Kd")))
        filenameIMG = getfield (MTL, "map_Kd");
        ## Check in image exists in working directory
        f_len = length(filenameIMG) - 4;
        if (exist (filenameIMG) == 2 &&
            ! strncmp (filenameIMG, filenameOBJ, f_len))
          new_fnIMG = strcat (filenameOBJ([1:end-4]), filenameIMG([end-3:end]));
          MTL.map_Kd = new_fnIMG;
          copyfile (filenameIMG, new_fnIMG);
        endif
      endif

      ## Make filename for material library file consistent with OBJ filename
      filenameMTL = strcat (filenameOBJ([1:end-4]), ".mtl");
      ## Write to mtl file
      mtlwrite (filenameMTL, MTL);
    endif

    ## Save scaling ratio and related measurements for each mesh object
    scale{i+1,1} = filenameOBJ;
    scale{i+1,2} = ratio;
    scale{i+1,3} = maxD;
    scale{i+1,4} = maxDnew;

    ## Flush the screen output to display the results
    ## during iterations through multiple mesh files
    page_screen_output (0);
    page_output_immediately (1);
    ## Display ratio and measurements for the current mesh object
    printf("Model %s was scaled by %f.\n", scale{i+1,1}, scale{i+1,2});
    printf("Old max distance was %f. New max distance is %f.\n\n", ...
            scale{i+1,3}, scale{i+1,4});
  endfor

  ## Check for output arguments
  if (nargout==0)
    ## Ask user for filename to save scaling values to csv file
    filenameCSV = uiputfile ({'*.csv', 'Supported Formats'});
    cell2csv (filenameCSV, scale);

  elseif (nargout==1)
    varargout{1} = scale;
  endif

endfunction

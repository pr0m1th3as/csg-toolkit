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
## @deftypefn{Function} longbone_Scaling
## @deftypefnx{Function} longbone_Scaling (@var{filename})
## @deftypefnx{Function} longbone_Scaling (@var{filename}, @var{maxD})
## @deftypefnx{Function} @var{scale} = longbone_Scaling
## @deftypefnx{Function} @var{scale} = longbone_Scaling (@var{...})
##
## This function reads the available triangular meshes stored in .obj format that
## are present in the working directory and utilizes the 'longbone_maxDistance'
## function to measure their maximum distance used for scaling and its
## corresponding points, which are stored in a Meshlab .pp file for visualization
## purposes.
##
## The function may be used for scaling a single mesh or batch scaling multiple
## meshes. When called with a single input argument, the function will scale the
## particular mesh defined in @var{filename} so that the scaled mesh will have
## a user-defined maximum distance. The function is primarily intended for
## long bone tringular mesh models produced with 3D photogrammetry, hence the
## use of maximum distance as the target reference for scaling. Nevertheless, it
## can be used for scaling any triangular mesh of arbitrary shape according to
## its maximum distance, which can either be parsed as a second input argument
## @var{maxD} or when prompted by the function. If no input argument is provided,
## the function will scan the working directory for all .obj files and will scale
## them iteratively in batch mode. For each file, the user will be prompted for
## the real maximum distance of each corresponding long bone and it will
## subsequently scale the model to the real world dimensions with units in mm.
## Scaling is performed about the model's barycentric coordinates, which are
## translated to origin.
##
## Prior to scaling, when prompted for manual measurement value, the user may
## examine the model's maximum distance corresponding points by opening the
## triangular mesh in Meshlab and utilizing the corresponding .pp file, which is
## created automatically by the present function. After scaling is performed for
## a given triangular mesh, the user may re-examine the scaled model by reloading
## the mesh and its respective .pp file, which is automatically updated with the
## new points of its maximum distance.
##
## The user may also declare an output argument, which will be returned as a cell
## array containing each model's filename and scale related measurements in the
## following order: "filename", "ratio", "oldMaxD", "newMaxD". If no output
## argument is declared, the function will prompt the user to specify a .csv file
## (new or existing, which will be overwritten) that all scaling measurements
## will be saved. The first row in the csv file will contain the column labels
## "filename", "ratio", "oldMaxD", "newMaxD", which are also saved in
## @var{scale(1,[1:4])}, when output variable is declared by calling the function
## as @var{scale} = longbone_Scaling.
## @seealso{longbone_maxDistance}
## @end deftypefn

function [varargout] = longbone_Scaling (varargin)  
  ## check for valid number of input variables
  if nargin ~= 0 && nargin ~= 1 && nargin ~= 2
    printf ("invalid number of input arguments\n");
    return;
  endif
  ## check if input argument is provided. It must be a string
  if nargin == 1
    ## it must be string
    if ! ischar (varargin{1}(:)')
      printf ("Mesh object filename must be a string\n");
      return;
    endif
    filenames(1).name = varargin{1}(:)';
    ask_user = true;
  elseif nargin == 2
    ## first input argument must be a string
    if ! ischar (varargin{1}(:)')
      printf ("Mesh object filename must be a string\n");
      return;
    endif
    ## second input argument must be a scalar
    if ! isscalar (varargin{2}(:))
      printf ("Max Distance must be a scalar\n");
      return;
    endif
    filenames(1).name = varargin{1}(:)';
    RmaxD = varargin{2}(:);
    ask_user = false;
  else
    ## list the filenames with .obj extension in the working folder
    filenames = dir ("*.obj");
    ask_user = true;
  endif
  ## initialize header for the cell array
  scale = {"filename", "ratio", "oldMaxD", "newMaxD"};
  ## perform scaling for each mesh object present in the working directory
  for i = 1:length (filenames)
    ## store filename of current mesh
    filenameOBJ = strcat (filenames(i).name);
    ## read original obj file
    [v, f, vt, ft, vn, fn, filenameMTL] = readObj (filenameOBJ);
    ## calculate maximum distance and corresponding points
    [maxD, p1, p2] = longbone_maxDistance (v(:,[1:3]));
    ## save maxD points to Meshlab .pp file using the name convention of the
    ## original mesh obj file
    filenamePP = strcat (filenameOBJ([1:end-4]), ".pp");
    write_MeshlabPoints (filenamePP, filenameOBJ, [p1; p2]);
    ## ask user for actual maximum distance (real size)
    if ask_user
      realMaxD = inputdlg ("maximum distance in mm", filenameOBJ, [1,20]);
      RmaxD = str2num (cell2mat (realMaxD));
    endif
    ## calculate scaling ratio
    ratio = RmaxD / maxD;
    ## find barycentric coordinates and translate mesh' origin
    origin = meshBarycenter (v(:,[1:3]), f);
    v(:,[1:3]) = v(:,[1:3]) - origin;
    ## scale mesh
    v(:,[1:3]) = v(:,[1:3]) * ratio;
    ## check if texture coordinates and material file exist in OBJ
    if (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL))
      ## read original mtl file
      MTL = mtlread (filenameMTL);
      ## find which material's name corresponds to OBJ name and keep it or
      ## alternatively search for texture altas in non-empty  'map_Kd' field
      for id = 1:length (MTL)
        if strcmp (getfield (MTL(id), "newmtl"), filenameOBJ([1:end-4]))
          newmtl_id = id;
        endif
        if isfield (MTL(id), "map_Kd") && !isempty (getfield (MTL(id), "map_Kd"))
          map_Kd_id = id;
        endif
      endfor
    endif
    ## save scaled model according to the elements present in the original OBJ
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
    ## calculate scaled model's max distance and save its points in .pp file
    [maxDnew, p1, p2] = longbone_maxDistance (v);
    write_MeshlabPoints (filenamePP, filenameOBJ, [p1; p2]);
    ## check if texture coordinates and material file exist in OBJ (again!!)
    if (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL))
      if exist ("newmtl_id")
        MTL = MTL(newmtl_id);
        % rename material to the OBJ's updated name (just in case!!)
        MTL.newmtl = filenameOBJ([1:end-4]);
      else
        MTL = MTL(map_Kd_id);
        % rename material to the OBJ's updated name (just in case!!)
        MTL.newmtl = filenameOBJ([1:end-4]);
      endif
      ## check for texture map image in 'map_Kd' field and if available, then
      ## compare its filename to OBJ filename and if they don't match make a new
      ## copy of the image with the updated OBJ filename. Also update map_Kd
      if isfield (MTL, "map_Kd") && ! isempty (getfield (MTL, "map_Kd"))
        filenameIMG = getfield (MTL, "map_Kd");
        ## check in image exists in working directory 
        f_len = length(filenameIMG - 4);
        if exist (filenameIMG) == 2 && ! strncmp (filenameIMG, filenameOBJ, f_len)
          new_fnIMG = strcat(filenameOBJ([1:end-4]), filenameIMG([end-3:end]));
          MTL.map_Kd = new_fnIMG;
          copyfile (filenameIMG, new_fnIMG);
        endif
      endif
      ## make filename for material library file consistent with OBJ filename
      filenameMTL = strcat (filenameOBJ([1:end-4]), ".mtl");
      ## write to mtl file
      mtlwrite (filenameMTL, MTL);
    endif
    ## save scaling ratio and related measurements for each mesh object in
    ## a data structure
    scale{i+1,1} = filenameOBJ;
    scale{i+1,2} = ratio;
    scale{i+1,3} = maxD;
    scale{i+1,4} = maxDnew;
    ## flush the screen output to display the results during iterations through
    ## multiple mesh files
    page_screen_output (0);
    page_output_immediately (1);
    ## display ratio and measurements for the current mesh object
    printf("Model %s was scaled by %f.\n", scale{i+1,1}, scale{i+1,2});
    printf("Old max distance was %f. New max distance is %f.\n\n", ...
            scale{i+1,3}, scale{i+1,4});
  endfor
  ## check for output arguments
  if (nargout==0)
    ## ask user for filename to save scaling values to csv file
    filenameCSV = uiputfile ({'*.csv', 'Supported Formats'});
    cell2csv (filenameCSV, scale);
  elseif (nargout==1)
    varargout{1} = scale;
  endif
endfunction
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn{Function} renameObj (@var{source}, @var{target})
##
## This function renames the file containing a triangular mesh in OBJ format and
## its associated material library and texture map sidecar files if available.
## Only single model in OBJ and single material in MTL files are supported.
## as @var{scale} = longbone_Scaling.
## @seealso{readObj, writeObj}
## @end deftypefn

function [varargout] = renameObj (varargin)  
  ## check for valid number of input variables
  if nargin ~= 2
    printf("invalid number of input arguments\n");
    return;
  endif
  ## check if both input arguments are strings
  if (! ischar(varargin{1}(:)') || ! ischar(varargin{2}(:)'))
    printf("Both source and target filenames must be strings\n");
    return;
  else
    source = varargin{1}(:)';
    target = varargin{2}(:)';
  endif
  ## check if source file exists
  if (! exist (source, "file"))
    printf("Source file does not exist\n");
    return;
  endif
  ## check if source target files are identical
  if (strcmp (source, target))
    printf("Target and source filenames are the same\n");
    return;
  endif
  ## read original obj file and delete it
  [v, f, vt, ft, vn, fn, filenameMTL] = readObj (source);
  delete (source);
  ## check if texture coordinates and material file exist in OBJ
  mat_lib = false;
  if (! isempty (vt) && ! isempty (ft) && ! isempty (filenameMTL))
    ## check if material library sidecar file exists
    if (! exist (filenameMTL, "file"))
      printf("Material library file does not exist\n");
      return;
    endif
    ## read original mtl file
    MTL = mtlread (filenameMTL);
    ## find which material's name corresponds to OBJ name and keep it or
    ## alternatively search for texture altas in non-empty  'map_Kd' field
    for id = 1:length (MTL)
      if strcmp (getfield (MTL(id), "newmtl"), source([1:end-4]))
        newmtl_id = id;
      endif
      if isfield (MTL(id), "map_Kd") && ! isempty (getfield (MTL(id), "map_Kd"))
        map_Kd_id = id;
      endif
    endfor
    mat_lib = true;
    ## delete original file
    delete (filenameMTL);
  endif
  ## save scaled model according to the elements present in the original OBJ
  if (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL) &&
      !isempty (vn) && !isempty (fn))
    target = writeObj(v, f, vt, ft, vn, fn, target);
  elseif (!isempty (vt) && !isempty (ft) && !isempty (filenameMTL) &&
      isempty (vn) && isempty (fn))
    target = writeObj(v, f, vt, ft, target);
  elseif (isempty (vt) && isempty (ft) && isempty (filenameMTL) &&
      !isempty (vn) && !isempty (fn))
    target = writeObj(v, f, vn, fn, target);
  else
    target = writeObj(v, f, target);
  endif
  ## check if texture coordinates and material file exist in OBJ (again!!)
  if mat_lib
    if exist ("newmtl_id")
      MTL = MTL(newmtl_id);
      ## rename material to the OBJ's updated name (just in case!!)
      MTL.newmtl = target([1:end-4]);
    elseif exist ("map_Kd_id")
      MTL = MTL(map_Kd_id);
      ## rename material to the OBJ's updated name (just in case!!)
      MTL.newmtl = target([1:end-4]);
    endif
    ## check for texture map image in 'map_Kd' field and if available, then
    ## compare its filename to OBJ filename and if they don't match update map_Kd
    if (isfield (MTL, "map_Kd") && ! isempty (getfield (MTL, "map_Kd")))
      filenameIMG = getfield (MTL, "map_Kd");
      ## check in image exists in working directory 
      f_len = length(filenameIMG - 4);
      if exist (filenameIMG) == 2 && ! strncmp (filenameIMG, target, f_len)
        new_fnIMG = strcat (target([1:end-4]), filenameIMG([end-3:end]));
        MTL.map_Kd = new_fnIMG;
        ## check in image exists and rename it if necessary
        if (exist (filenameIMG, "file") && ! strcmp (filenameIMG, new_fnIMG))
          rename(filenameIMG, new_fnIMG);
        endif
      endif
    endif
    ## make filename for material library file consistent with OBJ filename
    filenameMTL = strcat (target([1:end-4]), ".mtl");
    ## write to mtl file
    mtlwrite(filenameMTL, MTL);
  endif
  ## check for output arguments
  if (nargout==1)
    ## return the final name of target file
    varargout{1} = target;
  endif
endfunction
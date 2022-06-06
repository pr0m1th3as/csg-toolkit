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
## @deftypefn{Function} @var{MLP} = read_MeshlabPoints (@var{filename})
## @deftypefnx{Function} [@var{MLP}, @var{point_id}] = read_MeshlabPoints (@var{filename})
##
## This function reads a .pp MeshLab Point file and returns and Nx4 matrix
## containing the name for each point along with the corresponding x y z 
## coordinates in each row. If name is given as a numerical value then it is
## stored as numerical value. If the name of a point is alpharithmetic or 0
## then NaN is stored in the relevant row of the MPL array and the name is
## stored in a separate cell array as a character string.
## 
## If one output argument is provided then the function returns the Nx4 matrix,
## if two output arguments are defined then an Nx3 matrix containing only the
## x, y, z coordinates and a cell array containing the names of each point are
## returned.
##
## @itemize
## @item
## @var{filename} is char string of the filename of Meshlab's Point filename.
##
## @item @var{MLP} contains the Picked Points in an Nx4 array following the format
## [name, x, y, z].
## @end itemize
## 
## If the names of the Picked Points are alphanumeric strings, then they can be
## retrieved in a separate variable:
##
## @itemize
## @item @var{MLP}(1,:) = [x, y, z]
## @item @var{point_id}(1) = "name of first point"
## @end itemize
## @seealso{read_MeshlabPoints}
## @end deftypefn

function [varargout] = read_MeshlabPoints (filename)
  ## declare output variables
  MLP = zeros (1,4);
  name_list = {''};
  ## scan file and identify line containing points
  point_index = 1;
  fid = fopen (filename, 'rt');
  line = fgets (fid);
  while ischar (line)
    if strfind (line, "<point")
      ## get indices of starting places for each axis coordinate
      x_start = strfind (line, "x=\"") + 3;
      y_start = strfind (line, "y=\"") + 3;
      z_start = strfind (line, "z=\"") + 3;
      name_start = strfind (line, "name=\"") + 6;
      quotes = strfind (line, "\"");
      ## get indices of ending places for each axis coordinate
      x_end = min (quotes(find (quotes > x_start))) - 1;
      y_end = min (quotes(find (quotes > y_start))) - 1;
      z_end = min (quotes(find (quotes > z_start))) - 1;
      name_end = min (quotes(find (quotes>name_start))) - 1;
      ## get values
      x = str2num (line(x_start:x_end));
      y = str2num (line(y_start:y_end));
      z = str2num (line(z_start:z_end));
      ## check for numeric values for each point's name field
      if (isnumeric (str2num (line(name_start:name_end))))
        name = str2num (line(name_start:name_end));
        name_list(point_index) = str2num (line(name_start:name_end));
      else
        name = NaN;
        name_list(point_index) = line(name_start:name_end);
      endif
      MLP(point_index,:) = [name x y z];
      point_index += 1;
    endif
    line = fgets (fid);
  end
  fclose (fid);
  ## return variables
  if (nargout==1)
    varargout{1} = MLP;
  elseif (nargout==2)
    varargout{1} = MLP(:,[2:4]);
    varargout{2} = name_list;
  endif
endfunction
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
## @deftypefn  {csg-toolkit} {@var{MLPP} =} read_MeshlabPoints (@var{filename})
## @deftypefnx {csg-toolkit} {[@var{MLPP}, @var{pnames}] =} read_MeshlabPoints (@var{filename})
##
## Read 3D coordinates from a MeshLab PickedPoints file.
##
## This function reads a @qcode{.pp} MeshLab PickedPoints file and returns an
## @math{Nx4} matrix containing the index for each point along with the
## corresponding @var{x},@var{y},@var{z} coordinates in each row.  The point
## name must be numerical, otherwise @qcode{NaN} is returned.  If the point
## names are stored as alpharithmetic strings, then these can be retrieved in
## a second output argument, @var{pnames}, as a cell array of strings.  In such
## case, @var{MLPP} is an @math{Nx3} matrix containing the
## @var{x},@var{y},@var{z} point coordinates, where @math{N} is the number of
## points.
##
## @seealso{read_MeshlabPoints}
## @end deftypefn

function [varargout] = read_MeshlabPoints (filename)

  ## Initialize output variables
  MLPP = zeros (1,4);
  name_list = {''};

  ## Scan file and identify line containing points
  point_index = 1;
  fid = fopen (filename, 'rt');
  line = fgets (fid);

  while (ischar (line))
    if (strfind (line, "<point"))
      ## Get indices of starting places for each axis coordinate
      x_start = strfind (line, "x=\"") + 3;
      y_start = strfind (line, "y=\"") + 3;
      z_start = strfind (line, "z=\"") + 3;
      name_start = strfind (line, "name=\"") + 6;
      quotes = strfind (line, "\"");

      ## Get indices of ending places for each axis coordinate
      x_end = min (quotes(find (quotes > x_start))) - 1;
      y_end = min (quotes(find (quotes > y_start))) - 1;
      z_end = min (quotes(find (quotes > z_start))) - 1;
      name_end = min (quotes(find (quotes>name_start))) - 1;

      ## Get coordinates
      x = str2num (line(x_start:x_end));
      y = str2num (line(y_start:y_end));
      z = str2num (line(z_start:z_end));

      ## Check for numeric values for each point's name field
      if (isnumeric (str2num (line(name_start:name_end))))
        name = str2num (line(name_start:name_end));
        name_list(point_index) = str2num (line(name_start:name_end));
      else
        name = NaN;
        name_list(point_index) = line(name_start:name_end);
      endif
      MLPP(point_index,:) = [name, x, y, z];
      point_index += 1;
    endif
    ## Go to next line
    line = fgets (fid);
  end
  fclose (fid);

  ## Return variables
  if (nargout==1)
    varargout{1} = MLPP;

  elseif (nargout==2)
    varargout{1} = MLPP(:,[2:4]);
    varargout{2} = name_list;
  endif

endfunction

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
## @deftypefn{Function} write_MeshlabPoints (@var{filename}, @var{meshname}, @var{MLP})
## @deftypefnx{Function} write_MeshlabPoints (@var{filename}, @var{meshname}, @var{MLP}, @var{point_id})
##
## This function writes the 3D coordinates of points along with their associated
## names in a .pp MeshLab Point file. The function requires at leat three input
## arguments. If three arguments are provided, the first two should be character
## strings with the filename under which the points will be saved and the filename
## of the associated mesh. The third input argument should be an Nx4 matrix
## containing name only numeric values) and x,y,z coordinates for each point.
##
## If the point names need be alphanumeric, then an additional variable may be
## parsed into the function, which should be a cell array of strings. The point
## names parsed as a fourth variable should correspond to the points provided in
## the second input argument, which should be an Nx3 matrix containing only the
## 3D coordinates of the points. In case of four input arguments where the third
## argument is an Nx4 matrix, the first column (arithmetic names) is ignored and
## the name list in the fourth argument is used.
## @seealso{write_MeshlabPoints}
## @end deftypefn

function write_MeshlabPoints (varargin)
  ## check the number of input variables
  if length (varargin) < 3 || length (varargin) > 4
    error 'wrong number of input arguments';
  endif
  ## check first two arguments are strings
  if ! ischar (varargin{1}(:)') || ! ischar (varargin{2}(:)')
    error 'first two input arguments should be char strings';
  else
    filename = varargin{1}(:)';
    meshname = varargin{2}(:)';
  endif
  ## chech third input argument
  if length (varargin) == 3 && size (varargin{3},2) == 4
    MLP = varargin{3}(:,[2:4]);
    pointindex = varargin{3}(:,1);
  elseif length (varargin) == 3 && size (varargin{3},2) == 3
    pointindex = [1:size(varargin{3},1)];
    MLP = varargin{3};
  elseif length (varargin) == 4 && size (varargin{3},2) == 3
    namelist = varargin{4};
    MLP = varargin{3};
  elseif length (varargin) == 4 && size (varargin{3},2) == 4
    namelist = varargin{4};
    MLP = varargin{3}(:,[2:4]);
  else
    error 'Invalid arguments';
  endif
  
  ## open .pp file in writing mode and append the required headers
  fid = fopen (filename,'wt');
  fprintf (fid, "<!DOCTYPE PickedPoints>\n<PickedPoints>\n <DocumentData>\n");
  ## get time, date and user from the system to use it for DocumentData section
  [a, user] = system ("users");
  user = user(1:end-1);     % remove trailing newline char from user string
  a = clock;
  a(6) = ceil (a(6));
  fprintf (fid, "  <DateTime time=""%02d:%02d:%02d"" date=""%d-%02d-%02d""/>\n",...
              a(4),a(5),a(6),a(1),a(2),a(3));
  fprintf (fid, "  <User name=""%s""/>\n", user);
  fprintf (fid, "  <DataFileName name=""%s""/>\n", meshname);
  fprintf (fid, "  <templateName name=""""/>\n </DocumentData>\n");
  ## check for three input arguments and add the points to the file
  if exist ('pointindex')
    for i = 1:length (pointindex)
      fprintf (fid, " <point active=""1"" name=""%d"" x=""%0.4f"" y=""%0.4f"" z=""%0.4f""/>\n",...
               pointindex(i), MLP(i,1), MLP(i,2), MLP(i,3));
    endfor
  endif
  ## check for three input arguments and add the points to the file
  if exist('namelist')
    for i = 1:length (namelist)
      fprintf (fid, " <point active=""1"" name=""%s"" x=""%0.4f"" y=""%0.4f"" z=""%0.4f""/>\n",...
               namelist{i}, MLP(i,1), MLP(i,2), MLP(i,3));
    endfor
  endif
  fprintf (fid, "</PickedPoints>");
  fclose (fid);
endfunction
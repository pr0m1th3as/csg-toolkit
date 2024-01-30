## Copyright (C) 2018-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {} write_MeshlabPoints (@var{filename}, @var{meshname}, @var{MLPP})
## @deftypefnx {csg-toolkit} {} write_MeshlabPoints (@var{filename}, @var{meshname}, @var{MLPP}, @var{pnames})
##
## Write 3D coordinates in a MeshLab PickedPoints file.
##
## This function writes the 3D coordinates of points along with their associated
## names in a @qcode{.pp} MeshLab PickedPoints file. The function requires at
## least three input arguments.  When three arguments are provided, the first
## two must be character strings with the filename under which the points will
## be saved, @var{filename} and the filename of the associated 3D mesh,
## @var{meshname}.  The third input argument must be an @math{Nx4} matrix with
## each row containing the point index (only numeric values) along with the
## @var{x},@var{y},@var{z} coordinates of each point.
##
## If the point index is required to be alphanumeric, then an additional input
## ardument may be parsed, @var{pnames}, which must be a cell array of
## strings.  @var{pnames} must correspond to the points provided in the second
## input argument, in which case it can be an @math{Nx3} matrix containing only
## the 3D coordinates of the points.  In case of four input arguments, where the
## third argument is an Nx4 matrix, the first column (point index) is ignored
## and the point names in @var{pnames} is used.
##
## @seealso{write_MeshlabPoints}
## @end deftypefn

function write_MeshlabPoints (varargin)

  ## Check the number of input variables
  if (length (varargin) < 3 || length (varargin) > 4)
    error ("write_MeshlabPoints: wrong number of input arguments.");
  endif
  ## Check first two arguments are strings
  if ! ischar (varargin{1}(:)') || ! ischar (varargin{2}(:)')
    error (strcat (["write_MeshlabPoints: first two input arguments"], ...
                   [" must be char strings."]));
  else
    filename = varargin{1}(:)';
    meshname = varargin{2}(:)';
  endif

  ## Check third input argument
  if (length (varargin) == 3 && size (varargin{3},2) == 4)
    MLPP = varargin{3}(:,[2:4]);
    pointindex = varargin{3}(:,1);

  elseif (length (varargin) == 3 && size (varargin{3},2) == 3)
    pointindex = [1:size(varargin{3},1)];
    MLPP = varargin{3};

  elseif (length (varargin) == 4 && size (varargin{3},2) == 3)
    namelist = varargin{4};
    MLPP = varargin{3};

  elseif (length (varargin) == 4 && size (varargin{3},2) == 4)
    namelist = varargin{4};
    MLPP = varargin{3}(:,[2:4]);

  else
    error 'Invalid arguments';
  endif

  ## Open .pp file in writing mode and append the required headers
  fid = fopen (filename,'wt');
  fprintf (fid, "<!DOCTYPE PickedPoints>\n<PickedPoints>\n <DocumentData>\n");

  ## Get time, date and user from the system to use it for DocumentData section
  [a, user] = system ("users");
  user = user(1:end-1);     # remove trailing newline char from user string
  a = clock;
  a(6) = ceil (a(6));
  fprintf (fid, "  <DateTime time=""%02d:%02d:%02d"" date=""%d-%02d-%02d""/>\n",...
              a(4),a(5),a(6),a(1),a(2),a(3));
  fprintf (fid, "  <User name=""%s""/>\n", user);
  fprintf (fid, "  <DataFileName name=""%s""/>\n", meshname);
  fprintf (fid, "  <templateName name=""""/>\n </DocumentData>\n");

  ## Check for three input arguments and add the points to the file
  if exist ('pointindex')
    for i = 1:length (pointindex)
      fprintf (fid, " <point active=""1"" name=""%d"" x=""%0.4f"" y=""%0.4f"" z=""%0.4f""/>\n",...
               pointindex(i), MLPP(i,1), MLPP(i,2), MLPP(i,3));
    endfor
  endif

  ## Check for three input arguments and add the points to the file
  if exist('namelist')
    for i = 1:length (namelist)
      fprintf (fid, " <point active=""1"" name=""%s"" x=""%0.4f"" y=""%0.4f"" z=""%0.4f""/>\n",...
               namelist{i}, MLPP(i,1), MLPP(i,2), MLPP(i,3));
    endfor
  endif

  ## Close file
  fprintf (fid, "</PickedPoints>");
  fclose (fid);

endfunction

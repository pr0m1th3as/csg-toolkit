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
## @deftypefn{Private Function} @var{mtl} = mtlread (@var{filename})
##
## This function reads the basic elements from a Wavefront material file and
## returns their values in a structure. Input argument filename should be
## a char string with the filename of the material library file.
##
## It recognizes only material name statements, i.e. newmtl; material color and
## illumination statements, i.e. Ka, Kd, Ks, Tr, d, illum and Ns;
## as well as texture map statements, i.e. map_Kd, map_Kd and map_Ks. If multiple
## materials are present, then each material is stored with a diferent index
## along with its corresponding elements as fields of the structure. Each field
## is created only if the corresponding element is present in the .mtl file.
## @seealso{mtlwrite}
## @end deftypefn

function mtl = mtlread (filename)
  mtl_index = 0;
  fid = fopen (filename,'rt');
  line = fgets (fid);
  while ischar (line)
    if (line(1) == 'n' && line(2) == 'e' && line(3) == 'w' && line(4) == 'm' ...
                  && line(5) == 't' && line(6) == 'l')
      mtl_index = mtl_index + 1;
      mtl(mtl_index).newmtl = sscanf(line, 'newmtl %s', 1);
    endif
    if (line(1) == 'K' && line(2) == 'a')
      mtl(mtl_index).Ka = sscanf(line, 'Ka %f %f %f', 3);
    endif
    if (line(1) == 'K' && line(2) == 'd')
      mtl(mtl_index).Kd = sscanf(line, 'Kd %f %f %f', 3);
    endif
    if (line(1) == 'K' && line(2) == 's')
      mtl(mtl_index).Ks = sscanf(line, 'Ks %f %f %f', 3);
    endif
    if (line(1) == 'T' && line(2) == 'r')
      mtl(mtl_index).Tr = sscanf(line, 'Tr %f', 1);
    endif
    if (line(1) == 'd')
      mtl(mtl_index).Tr = 1 - sscanf(line, 'd %f', 1);
    endif
    if (line(1) == 'i' && line(2) == 'l' && line(3) == 'l' && line(4) == 'u' ...
                  && line(5) == 'm')
      mtl(mtl_index).illum = sscanf(line, 'illum %d', 1);
    endif
    if (line(1) == 'N' && line(2) == 's')
      mtl(mtl_index).Ns = sscanf(line, 'Ns %f', 1);
    endif
    if (line(1) == 'm' && line(2) == 'a' && line(3) == 'p' && line(4) == '_' ...
                  && line(5) == 'K' && line(6) == 'a')
      mtl(mtl_index).map_Kd = sscanf(line, 'map_Ka %s', 1);
    endif
    if (line(1) == 'm' && line(2) == 'a' && line(3) == 'p' && line(4) == '_' ...
                  && line(5) == 'K' && line(6) == 'd')
      mtl(mtl_index).map_Kd = sscanf(line, 'map_Kd %s', 1);
    endif
    if (line(1) == 'm' && line(2) == 'a' && line(3) == 'p' && line(4) == '_' ...
                  && line(5) == 'K' && line(6) == 's')
      mtl(mtl_index).map_Kd = sscanf(line, 'map_Ks %s', 1);
    endif
    line = fgets (fid);
  end
  fclose (fid);
endfunction
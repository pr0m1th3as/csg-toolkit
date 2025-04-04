## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {csg-toolkit} {@var{models} =} load_classifiers (@var{bone})
##
## This function loads all classifiers related to @var{bone} into a cell array.
## @end deftypefn

function models = load_classifiers (bone)
  models = cell (65,1);
  if (strcmpi (bone, "Femur"))
    for i = 1:65
      name = fullfile (pwd, "private", sprintf ("F%02d.mat", i));
      mdl = loadmodel (name);
      models(i) = {mdl};
    endfor
  elseif (strcmpi (bone, "Humerus"))
    for i = 1:65
      name = fullfile (pwd, "private", sprintf ("H%02d.mat", i));
      mdl = loadmodel (name);
      models(i) = {mdl};
    endfor
  elseif (strcmpi (bone, "Tibia"))
    for i = 1:65
      name = fullfile (pwd, "private", sprintf ("T%02d.mat", i));
      mdl = loadmodel (name);
      models(i) = {mdl};
    endfor
  elseif (strcmpi (bone, "Ulna"))
    for i = 1:65
      name = fullfile (pwd, "private", sprintf ("U%02d.mat", i));
      mdl = loadmodel (name);
      models(i) = {mdl};
    endfor
  endif
endfunction

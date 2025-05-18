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
## @deftypefn  {csg-toolkit} {@var{h} =} plot_CrossSections (@var{EXTRA})
##
## Plot the default cross sectional contours from a long bone sample.
##
## This fuction takes a structure input argument, @var{EXTRA}, which can be the
## fourth output argument of the @code{longbone_Geometry} function and plots all
## five cross sectional polygon in the a single figure. Unlike the
## @code{visualize_CrossSections} function, @code{plot_CrossSections} does not
## print a legend with CSG properties for each cross section.
##
## Use this function to swiftly evaluate that a 3D model processed with the
## @code{longbone_Geometry} function does not contain any artifacts, which may
## lead to degenerate closed polygons representing the computed diaphyseal cross
## sections and as a result to incorrect values of CSG properties.
##
## The input argument @var{EXTRA} can also be any structure with the following
## fields, as long as they contain equivalent values to those returned by
## @code{longbone_Geometry}.
##
## @enumerate
## @item @qcode{filename}
## @item @qcode{poly2D_20}
## @item @qcode{poly2D_35}
## @item @qcode{poly2D_50}
## @item @qcode{poly2D_65}
## @item @qcode{poly2D_80}
## @end enumerate
##
## Use the @code{plot_CrossSections} function to visualize the cross sections of
## the large collection of long bones from European populations provided by the
## @qcode{csg-dataset} package.
##
## @seealso{longbone_Geometry, visualize_CrossSections}
## @end deftypefn

function h = plot_CrossSections (EXTRA)
  ## Open new figure
  figure ("numbertitle", "off", "menubar", "none", "position", [400 0 1600 400]);
  ## Process all samples
  samples = length (EXTRA);
  for i = 1:samples
    ## Name figure after filename
    set (gcf, 'name', EXTRA(i).filename);
    ## Plot cross sections of current sample
    h(1) = subplot (1, 5, 1);
    x = [EXTRA(i).poly2D_20(:,1); EXTRA(i).poly2D_20(1,1)];
    y = [EXTRA(i).poly2D_20(:,2); EXTRA(i).poly2D_20(1,2)];
    plot (x, y);
    axis ("image");
    title ("Cross section at 20%", 'FontSize', 16);
    xlabel ('Frontal axis (mm)','FontSize', 10);
    ylabel ('Sagital axis (mm)','FontSize', 10);
    grid ("on"); xlim ("padded"); ylim ("padded");
    alim = [xlim, ylim];
    h(2) = subplot (1, 5, 2);
    x = [EXTRA(i).poly2D_35(:,1); EXTRA(i).poly2D_35(1,1)];
    y = [EXTRA(i).poly2D_35(:,2); EXTRA(i).poly2D_35(1,2)];
    plot (x, y);
    axis ("image");
    title ("Cross section at 35%", 'FontSize', 16);
    xlabel ('Frontal axis (mm)','FontSize', 10);
    ylabel ('Sagital axis (mm)','FontSize', 10);
    grid ("on"); xlim ("padded"); ylim ("padded");
    alim = [alim, xlim, ylim];
    h(3) = subplot (1, 5, 3);
    x = [EXTRA(i).poly2D_50(:,1); EXTRA(i).poly2D_50(1,1)];
    y = [EXTRA(i).poly2D_50(:,2); EXTRA(i).poly2D_50(1,2)];
    plot (x, y);
    axis ("image");
    title ("Cross section at 50%", 'FontSize', 16);
    xlabel ('Frontal axis (mm)','FontSize', 10);
    ylabel ('Sagital axis (mm)','FontSize', 10);
    grid ("on"); xlim ("padded"); ylim ("padded");
    alim = [alim, xlim, ylim];
    h(4) = subplot (1, 5, 4);
    x = [EXTRA(i).poly2D_65(:,1); EXTRA(i).poly2D_65(1,1)];
    y = [EXTRA(i).poly2D_65(:,2); EXTRA(i).poly2D_65(1,2)];
    plot (x, y);
    axis ("image");
    title ("Cross section at 65%", 'FontSize', 16);
    xlabel ('Frontal axis (mm)','FontSize', 10);
    ylabel ('Sagital axis (mm)','FontSize', 10);
    grid ("on"); xlim ("padded"); ylim ("padded");
    alim = [alim, xlim, ylim];
    h(5) = subplot (1, 5, 5);
    x = [EXTRA(i).poly2D_80(:,1); EXTRA(i).poly2D_80(1,1)];
    y = [EXTRA(i).poly2D_80(:,2); EXTRA(i).poly2D_80(1,2)];
    plot (x, y);
    axis ("image");
    title ("Cross section at 80%", 'FontSize', 16);
    xlabel ('Frontal axis (mm)','FontSize', 10);
    ylabel ('Sagital axis (mm)','FontSize', 10);
    grid ("on"); xlim ("padded"); ylim ("padded");
    alim = [alim, xlim, ylim];
    ## Fix axes to identical limits
    minlim = min (alim);
    maxlim = max (alim);
    xlim (h(1), [minlim, maxlim]); ylim (h(1), [minlim, maxlim]);
    xlim (h(2), [minlim, maxlim]); ylim (h(2), [minlim, maxlim]);
    xlim (h(3), [minlim, maxlim]); ylim (h(3), [minlim, maxlim]);
    xlim (h(4), [minlim, maxlim]); ylim (h(4), [minlim, maxlim]);
    xlim (h(5), [minlim, maxlim]); ylim (h(5), [minlim, maxlim]);
    ## Get user confirmation
    if (samples > 1)
      btn = questdlg ("Do you want to continue?", ...
                      EXTRA(i).filename, "Yes", "No");
      if (strcmp (btn, "Yes"))
        if (i == samples)
          close (gcf);
          h = NaN;
        endif
      elseif (strcmp (btn, "No"))
        break;
      else
        close (gcf);
        h = NaN;
        break;
      endif
    endif
  endfor
endfunction

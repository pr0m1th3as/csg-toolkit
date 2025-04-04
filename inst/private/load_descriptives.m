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
## @deftypefn  {csg-toolkit} {[@var{idx}, @var{mu}, @var{sigma}, @var{Mlo}, @
## @var{Mhi}, @var{Flo}, @var{Fhi}] =} load_descriptives (@var{bone})
##
## This function returns the descriptives related to @var{bone}.
##
## @itemize
## @item @var{idx} contains the indices to @var{DATA} used by the classifiers.
## @item @var{mu} contains the mean values of the indexed measurements.
## @item @var{sigma} contains the standard deviations of the indexed
## measurements.
## @item @var{Mlo} contains the lower threshold (Q1 - 1.5*IQR) values to detect
## outliers for males.
## @item @var{Mhi} contains the upper threshold (Q3 + 1.5*IQR) values to detect
## outliers for males.
## @item @var{Flo} contains the lower threshold (Q1 - 1.5*IQR) values to detect
## outliers for females.
## @item @var{Fhi} contains the upper threshold (Q3 + 1.5*IQR) values to detect
## outliers for females.
## @end itemize
##
## @end deftypefn

function [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives (bone)
  D = load (fullfile (pwd, "private", 'descriptives.mat'));
  if (strcmpi (bone, "Femur"))
    idx = D.Fvars(:,1)';
    mu = D.Fvars(:,2)';
    sigma = D.Fvars(:,3)';
    Mlo = D.Fvars(:,4)';
    Mhi = D.Fvars(:,5)';
    Flo = D.Fvars(:,6)';
    Fhi = D.Fvars(:,7)';
  elseif (strcmpi (bone, "Humerus"))
    idx = D.Hvars(:,1)';
    mu = D.Hvars(:,2)';
    sigma = D.Hvars(:,3)';
    Mlo = D.Hvars(:,4)';
    Mhi = D.Hvars(:,5)';
    Flo = D.Hvars(:,6)';
    Fhi = D.Hvars(:,7)';
  elseif (strcmpi (bone, "Tibia"))
    idx = D.Tvars(:,1)';
    mu = D.Tvars(:,2)';
    sigma = D.Tvars(:,3)';
    Mlo = D.Tvars(:,4)';
    Mhi = D.Tvars(:,5)';
    Flo = D.Tvars(:,6)';
    Fhi = D.Tvars(:,7)';
  elseif (strcmpi (bone, "Ulna"))
    idx = D.Uvars(:,1)';
    mu = D.Uvars(:,2)';
    sigma = D.Uvars(:,3)';
    Mlo = D.Uvars(:,4)';
    Mhi = D.Uvars(:,5)';
    Flo = D.Uvars(:,6)';
    Fhi = D.Uvars(:,7)';
  endif
endfunction

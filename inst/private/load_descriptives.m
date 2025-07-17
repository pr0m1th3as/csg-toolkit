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
## @deftypefn  {csg-toolkit} {@var{varargout} = } load_descriptives (@var{bone}, @var{type})
## @deftypefnx {csg-toolkit} {[@var{idx}, @var{mu}, @var{sigma}, @var{Mlo}, @
## @var{Mhi}, @var{Flo}, @var{Fhi}] =} load_descriptives (@var{bone}, @qcode{'sex'})
## @deftypefnx {csg-toolkit} {[@var{mu}, @var{sigma}, @var{lower}, @var{upper}, @
## @var{idx}] =} load_descriptives (@var{bone}, @qcode{'pair'})
##
## This function returns the descriptives related to @var{bone} according to the
## requested @var{type}, which can be either @qcode{'sex'} or @qcode{'pair'}.
##
## When @var{type} is @qcode{'sex'}, the following output arguments, which are
## utilized by the @qcode{longbone_Sex} function, are returned:
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
## When @var{type} is @qcode{'pair'}, the following output arguments, which are
## utilized by the @qcode{longbone_Pair} function, are returned:
##
## @itemize
## @item @var{mu} contains the mean differences of absolute values between left
## and right sides for each measurement.
## @item @var{sigma} contains the standard deviations of the differences of the
## absolute values between left and right sides for each measuremens.
## @item @var{lower} contains the lower threshold below which a tested pair can
## be considered a definite mismatch.
## @item @var{upper} contains the upper threshold above which a tested pair can
## be considered a definite mismatch.
## @item @var{idx} contains a logical vector of the measuremens to be used in
## pair matching.
## @end itemize
##
## @end deftypefn

function [varargout] = load_descriptives (bone, type)
  if (strcmpi (type, 'sex'))
    D = load (fullfile (pwd, 'private', 'descriptives_sex.mat'));
    if (strcmpi (bone, 'Femur'))
      varargout{1} = D.Fvars(:,1)';
      varargout{2} = D.Fvars(:,2)';
      varargout{3} = D.Fvars(:,3)';
      varargout{4} = D.Fvars(:,4)';
      varargout{5} = D.Fvars(:,5)';
      varargout{6} = D.Fvars(:,6)';
      varargout{7} = D.Fvars(:,7)';
    elseif (strcmpi (bone, 'Humerus'))
      varargout{1} = D.Hvars(:,1)';
      varargout{2} = D.Hvars(:,2)';
      varargout{3} = D.Hvars(:,3)';
      varargout{4} = D.Hvars(:,4)';
      varargout{5} = D.Hvars(:,5)';
      varargout{6} = D.Hvars(:,6)';
      varargout{7} = D.Hvars(:,7)';
    elseif (strcmpi (bone, 'Tibia'))
      varargout{1} = D.Tvars(:,1)';
      varargout{2} = D.Tvars(:,2)';
      varargout{3} = D.Tvars(:,3)';
      varargout{4} = D.Tvars(:,4)';
      varargout{5} = D.Tvars(:,5)';
      varargout{6} = D.Tvars(:,6)';
      varargout{7} = D.Tvars(:,7)';
    elseif (strcmpi (bone, 'Ulna'))
      varargout{1} = D.Uvars(:,1)';
      varargout{2} = D.Uvars(:,2)';
      varargout{3} = D.Uvars(:,3)';
      varargout{4} = D.Uvars(:,4)';
      varargout{5} = D.Uvars(:,5)';
      varargout{6} = D.Uvars(:,6)';
      varargout{7} = D.Uvars(:,7)';
    endif
  elseif (strcmpi (type, 'pair'))
    D = load (fullfile (pwd, 'private', 'descriptives_pair.mat'));
    if (strcmpi (bone, 'Femur'))
      varargout{1} = D.Fvars(1,:);
      varargout{2} = D.Fvars(2,:);
      varargout{3} = D.Fvars(3,:);
      varargout{4} = D.Fvars(4,:);
      varargout{5} = logical (D.Fvars(5,:));
    elseif (strcmpi (bone, 'Humerus'))
      varargout{1} = D.Hvars(1,:);
      varargout{2} = D.Hvars(2,:);
      varargout{3} = D.Hvars(3,:);
      varargout{4} = D.Hvars(4,:);
      varargout{5} = logical (D.Hvars(5,:));
    elseif (strcmpi (bone, 'Tibia'))
      varargout{1} = D.Tvars(1,:);
      varargout{2} = D.Tvars(2,:);
      varargout{3} = D.Tvars(3,:);
      varargout{4} = D.Tvars(4,:);
      varargout{5} = logical (D.Tvars(5,:));
    elseif (strcmpi (bone, 'Ulna'))
      varargout{1} = D.Uvars(1,:);
      varargout{2} = D.Uvars(2,:);
      varargout{3} = D.Uvars(3,:);
      varargout{4} = D.Uvars(4,:);
      varargout{5} = logical (D.Uvars(5,:));
    endif
  endif
endfunction

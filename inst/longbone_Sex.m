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
## @deftypefn  {csg-toolkit} {@var{table} =} longbone_Sex (@var{filename})
## @deftypefnx {csg-toolkit} {@var{table} =} longbone_Sex (@var{filename}, @var{bones})
## @deftypefnx {csg-toolkit} {@var{table} =} longbone_Sex (@var{folder}, @var{filename})
## @deftypefnx {csg-toolkit} {@var{table} =} longbone_Sex (@var{folder}, @var{filename}, @var{bones})
## @deftypefnx {csg-toolkit} {[@var{sex}, @var{prob}, @var{rep}, @var{vars}, @var{outliers}] =} longbone_Sex (@dots{})
##
## This function analyzes the geometry of an intact humerus, ulna, femur, or
## tibia bone using the @code{longbone_Geometry} function and estimates the
## biological sex of the individual the bone belongs to.
##
## @code{@var{table} = longbone_Sex (@var{filename})} will analyze the 3D model
## specified in @var{filename} provided that it conforms to the Wavefront OBJ
## file format and it is a pure triangular mesh.  The function will
## automatically determine the type of long bone in @var{filename} and register
## the required initial alignment points according to the bone type.  The 3D
## model is assumed to be in @math{mm} units and present in the working
## directory.  When called with a single output argument, @code{longbone_Sex}
## returns a @math{1x7} table, @var{table}, with the following variables:
##
## @enumerate
## @item @code{Name}: a cellstr variable with the filename of the analyzed bone.
## @item @code{Bone}: a cellstr variable with the type of bone as identified by
## the @code{longbone_Registration} function.
## @item @code{Sex}: a categorical variable with the estimated sex.
## @item @code{Prob}: a numerical variable with the posterior probability of the
## estimated sex.
## @item @code{Typical}: a logical variable specifying whether the measurements
## of the analyzed bone are typical of training sample.
## @item @code{Used}: a cell variable specifying the variables used for
## assigning sex in a cell array of character vectors.
## @item @code{Outliers}: a cell variable specifying which of the measurements
## used for assigning sex are outliers in a cell array of character vectors.
## @end enumerate
##
## @var{folder}, which must be a character vector, defines the relative or
## absolute path to the directory containing the 3D model in @var{filename}.
## When omitted, the current working directory is assumed.
##
## @var{bones} must be cell array of character vectors specifying one or more
## long bones that should be analyzed.  @code{longbone_Geometry} will only
## analyze the 3D model if it matches one of the bones named in @var{bones},
## unless the @qcode{'Force'} option is also selected, in which case the user
## must specify a single bone such as in @code{@{'Force', 'Humerus'@}}.  Valid
## options are:
##
## @enumerate
## @item @qcode{'Humerus'}
## @item @qcode{'Ulna'}
## @item @qcode{'Femur'}
## @item @qcode{'Tibia'}
## @item @qcode{'All'}
## @item @qcode{'Force'}
## @end enumerate
##
## Note: When @qcode{'Force'} is included, initial alignment points are also
## required and it is the user's responsibility to ensure the appropriate bone
## is being analyzed.  This option is provided so that damaged or heaviliy
## deformed bones, which may turn out as undefined, can be analyzed.  For more
## information look at the documentation of the @code{longbone_Geometry}
## function by typing @code{help longbone_Geometry}.
##
## When called with more than one output argument, such as in @code{[@var{sex},
## @var{prob}, @var{typ}, @var{vars}, @var{outliers}] = longbone_Sex (@dots{})},
## then the following output arguments as returned:
##
## @itemize
## @item @var{sex} is a numeric scalar value specifying the estimated sex, where
## @qcode{1} corresponds to male, @qcode{2} to female, and @qcode{0} to
## undetermined.
## @item @var{prob} is the posterior probability of the estimated sex.  Unless
## sex is assigned, @var{prob} is @qcode{NaN}.
## @item @var{typ} is a logical scalar value specifying whether the bone
## measurements used for sex estimation are typical @qcode{true} of the
## population samples used to train the classifiers or they contain outliers,
## in which case the sex estimation algorithm defaults to using only the
## univariate linear discriminant classifiers for estimating sex and @var{typ}
## is @qcode{false}.
## @item @var{vars} is a cell array of character vectors containing the
## variables that were significant for assigning sex.  These can be any
## combination of the utilized measurements, corresponding to univariate
## classifiers, as well as @qcode{'LDA ALL'}, @qcode{'KNN ALL'}, and
## @qcode{'FCNN Scores'}, corresponding to multivariate classifiers.
## @item @var{outliers} is a cell array of character vectors containing the
## measurements that were used in sex estimation but they are considered
## outliers according to the descriptive data of the training population
## samples.
## @end itemize
##
## @seealso{longbone_CustomGeometry, longbone_Registration}
## @end deftypefn

function [varargout] = longbone_Sex (varargin)

  ## Analyze 3D model
  [~, ~, BONE, ~, DATA] = longbone_Geometry (varargin{:});

  ## Keep filename
  filename = varargin{1};
  if (! (strcmpi (filename([end - 3:end]), '.obj')))
    filename = varargin{2};
  endif

  ## Load measurements
  measurements = longbone_Measurements ();

  ## Select variables for sex estimation
  if (strcmpi (BONE, 'Femur'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Femur');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = [0.95, 0.62, 0.86, 0.90, NaN];
  elseif (strcmpi (BONE, 'Humerus'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Humerus');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = [0.97, 0.69, 0.82, 0.90, NaN];
  elseif (strcmpi (BONE, 'Tibia'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Tibia');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = [0.94, 0.51, 0.52, 0.95, NaN];
  elseif (strcmpi (BONE, 'Ulna'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Ulna');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = [0.97, 0.55, 0.60, 0.96, NaN];
  else
    sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
    vars = {};
    outliers = {};
    if (nargout == 1 || nargout == 0)
      varNames = {'Name', 'Bone', 'Sex', 'Prob', 'Typical', 'Used', 'Outliers'};
      T = table ({filename}, {BONE}, sex, NaN, false, {vars}, {outliers}, ...
                 'VariableNames', varNames);
      varargout{1} = T;
    else
      varargout{1} = double (sex);
      if (nargout > 1)
        varargout{2} = 5;
      endif
      if (nargout > 2)
        varargout{3} = false;
      endif
      if (nargout > 3)
        varargout{4} = vars;
      endif
      if (nargout > 4)
        varargout{5} = outliers;
      endif
    endif
    return;
  endif

  ## Load classifiers and compute prediction scores
  models = load_classifiers (BONE);
  Scores = [];
  [~, score] = predict (models{1}, Z);
  Scores = [Scores, score];
  [~, score] = predict (models{2}, Z);
  Scores = [Scores, score];
  idx = 3;
  for i = 1:31
    [~, score] = predict (models{idx}, Z(i));
    Scores = [Scores, score];
    idx += 1;
    [~, score] = predict (models{idx}, Z(i));
    Scores = [Scores, score];
    idx += 1;
  endfor
  [~, s] = predict (models{65}, Scores);

  ## For samples without any outliers in measurements
  Mrep = all (X > Mlo & X < Mhi);
  Frep = all (X > Flo & X < Fhi);
  if (Mrep || Frep)
    rep = true;
    ## Group scores according to their magnitude
    odd_cols = Scores([1:2:end]);
    evencols = Scores([2:2:end]);
    m10s = any (odd_cols == 1 & evencols != 1);
    f10s = any (evencols == 1 & odd_cols != 1);
    m09s = any (odd_cols >= 0.9 & odd_cols < 1);
    f09s = any (evencols >= 0.9 & evencols < 1);
    m08n = numel (find (odd_cols >= 0.8 & odd_cols < 0.9));
    f08n = numel (find (evencols >= 0.8 & evencols < 0.9));
    odd_LDAs = odd_cols(3:2:end);
    evenLDAs = evencols(3:2:end);
    m09 = odd_cols >= 0.9;
    f09 = evencols >= 0.9;
    m08 = odd_cols >= 0.8;
    f08 = evencols >= 0.8;

    ## Apply hierarchical decision model
    if (m10s && s(1) > s(2))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 1;
      idx = find (odd_cols == 1);
      vars = resolve_idx (idx, measurements);
    elseif (f10s && s(1) < s(2))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 1;
      idx = find (evencols == 1);
      vars = resolve_idx (idx, measurements);
    elseif (strcmpi (BONE, 'Femur') && sum (m09) > sum (f09))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 2;
      idx = find (odd_cols >= 0.9 & odd_cols < 1);
      vars = resolve_idx (idx, measurements);
    elseif (strcmpi (BONE, 'Femur') && sum (m09) < sum (f09))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 2;
      idx = find (evencols >= 0.9 & evencols < 1);
      vars = resolve_idx (idx, measurements);
    elseif (! strcmpi (BONE, 'Femur') && (m09s && ! f09s && s(1) > s(2)))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 2;
      idx = find (odd_cols >= 0.9 & odd_cols < 1);
      vars = resolve_idx (idx, measurements);
    elseif (! strcmpi (BONE, 'Femur') && (! m09s && f09s && s(1) < s(2)))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 2;
      idx = find (evencols >= 0.9 & evencols < 1);
      vars = resolve_idx (idx, measurements);
    else
      s = round (s * 100) / 100;
      gap = s(2) - s(1) < 0.25 & s(2) - s(1) >= 0;
      if (s(1) > s(2) || (gap && (m08n > f08n)))
        sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
        grp = 3;
        vars = {};
        if (gap && (m08n > f08n))
          idx = find (odd_cols >= 0.8 & odd_cols < 0.9);
          vars = [vars, resolve_idx(idx, measurements)];
        endif
        if (s(1) > s(2))
          vars = [vars, {'FCNN Scores'}];
        endif
      elseif (s(2) >= s(1) + 0.25 || (gap && (f08n > m08n)))
        sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
        grp = 3;
        vars = {};
        if (gap && (f08n > m08n))
          idx = find (evencols >= 0.8 & evencols < 0.9);
          vars = [vars, resolve_idx(idx, measurements)];
        endif
        if (s(2) >= s(1) + 0.25)
          vars = [vars, {'FCNN Scores'}];
        endif
      else
        sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
        grp = 5;
        vars = {};
      endif
    endif
    outliers = {};

  else
    rep = false;
    ## For samples with outliers work only with univariate LDA models
    odd_cols = Scores([1:2:end]);
    evencols = Scores([2:2:end]);
    odd_cols = odd_cols(3:2:end);
    evencols = evencols(3:2:end);
    m09 = odd_cols > 0.9;
    f09 = evencols > 0.9;
    if (sum (m09) > sum (f09))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 4;
      vars = measurements(find (m09));
      outliers = measurements(m09  & (X < Mlo | X > Mhi));
    elseif (sum (m09) < sum (f09))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 4;
      vars = measurements(find (f09));
      outliers = measurements(f09 & (X < Flo | X > Fhi));
    else
      sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
      grp = 5;
      vars = {};
      outliers = {};
    endif
  endif

  ## Prepare output
  if (nargout == 1 || nargout == 0)
    varNames = {'Name', 'Bone', 'Sex', 'Prob', 'Typical', 'Used', 'Outliers'};
    T = table ({filename}, {BONE}, sex, prob(grp), rep, {vars}, {outliers}, ...
               'VariableNames', varNames);
    varargout{1} = T;
    return;
  endif
  varargout{1} = double (sex);
  if (nargout > 1)
    varargout{2} = grp;
  endif
  if (nargout > 2)
    varargout{3} = rep;
  endif
  if (nargout > 3)
    varargout{4} = vars;
  endif
  if (nargout > 4)
    varargout{5} = outliers;
  endif

endfunction

function used = resolve_idx (idx, measurements)
  used = {};
  if (idx(1) == 1)
    used = [used, {'LDA All'}];
    idx(1) = [];
  endif
  if (! isempty (idx))
    if (idx(1) == 2)
      used = [used, {'KNN All'}];
      idx(1) = [];
    endif
  endif
  if (! isempty (idx))
    if (idx(1) == 3)
      idx(1) = 4;
    endif
    used = [used, measurements(unique (floor ((idx - 2) / 2)))];
  endif
endfunction

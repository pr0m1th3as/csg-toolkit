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
## @deftypefn  {csg-toolkit} {@var{sorted} =} longbone_Pair (@var{DATA}, @var{BONE})
## @deftypefnx {csg-toolkit} {@var{sorted} =} longbone_Pair (@var{folder}, @var{BONE})
## @deftypefnx {csg-toolkit} {[@var{sorted}, @var{stats}, @var{unsorted}] =} longbone_Pair (@dots{})
##
## This function analyzes the geometry of intact humerus, ulna, femur, or
## tibia bones and matches the bilateral elements that belong to the same
## individual.
##
## @code{@var{sorted} = longbone_Pair (@var{DATA}, @var{BONE})} will analyze the
## measurements in @var{DATA} for a specific long bone specified in @var{BONE}
## and it will sort the skeletal elements by pair-matching bilateral elements
## that belong to the same individual as well as identify single elements, whose
## complement is missing.  @var{DATA} can be an @math{Nx62} numeric matrix or an
## @math{Nx63} cell array, with each row corresponding to a skeletal element.
##
## @itemize
## @item As a numeric matrix, @var{DATA} uses the first column to identify
## the side of the skeletal element, where left side is specified as @qcode{1}
## and right side denoted by @qcode{2}, whereas the remaining 61 columns contain
## the long bone measurements in the same order as returned in the @var{DATA}
## output argument of the @code{longbone_Geometry} function.  The names of these
## measurements can be retrieved by the @code{longbone_Measurements} function.
##
## @item As a cell array, @var{DATA} uses the first column to specify the
## identifier of each skeletal element, the second column holds the numeric
## value corresponding to side, and the remaining columns hold the measurement
## data.
## @end itemize
##
## The second input argument, @var{BONE}, must be a character vector specifying
## the type of bone under analysis and can take one of following options:
##
## @enumerate
## @item @qcode{'Humerus'}
## @item @qcode{'Ulna'}
## @item @qcode{'Femur'}
## @item @qcode{'Tibia'}
## @end enumerate
##
## The matched skeletal elements are returned in the output table @var{sorted},
## which contains three variables (columns): @qcode{'Left side'}, @qcode{'Right
## side'} and @qcode{'Score'}.  The first two variables contain the identifiers
## of the skeletal elements sorted by individual per each row, and the last one
## holds the cummulative @math{z}-score value used for pair-matching (the lower
## the better) bilateral elements of the same individual.  Single skeletal
## elements, which are identified by means of elimination, hence their
## corresponding @math{z}-score value is @qcode{NaN}, are reported first at the
## top rows of the @var{sorted} table.  The table contains @qcode{'RowNames'}
## using the naming convention @qcode{"Individual 1", @dots{}, "Individual k"}
## up to the total number of distinct individuals identified in the @var{DATA}.
## When @var{DATA} is a numeric matrix, the skeletal element identifiers are
## automatically assigned as @qcode{"1", "2", @dots{}, "N"} with each number
## corresponding to the row each element is stored in.
##
## When a second output argument is requested, @code{longbone_Pair} returns an
## additional table, @var{stats}, with various statistics regarding the results
## of the pair-matching process specified in the following table variables:
##
## @enumerate
## @item @qcode{'Samples'} : the number of processed samples, which is
## equivalent to the number of rows in @var{DATA}.
## @item @qcode{'Sorted'} : the total number of skeletal elements, which have
## been sorted to distinct individuals.
## @item @qcode{'Paired'} the number of bilateral paired elements, which have
## been identified.
## @item @qcode{'Single'} : the number of single skeletal elements, which have
## been identified to distinct individuals.
## @item @qcode{'Individuals'} : the total number of distinct individuals
## identified in the skeletal assemblage.
## @item @qcode{'Unsorted'} : the number of skeletal elements, which could not
## be sorted.
## @item @qcode{'Plausible'} : the number of plausible pair-matches among the
## remaining skeletal elements, which were not sorted.
## @end enumerate
##
## A comprehensive list of the plausible pairs of unsorted elements may be
## returned as the third output argument, @var{unsorted}, which is also a table
## similar to the one returned for the sorted skeletal elements.  Each plausible
## pair is accompanied by the respective cummulative @math{z}-score value, which
## may are may not be indicative of which plausible pair is more likely to be a
## true match.  Keep in mind, that skeletal elements in the plausible pairs list
## may well be single elements without their bilateral counterpart being present
## in the data sample.  Also note, that the algorithm utilized by the
## @code{longbone_Pair} function favors precision over recall and thus is
## optimized for maximizing the ratio of true positives over the sum of true
## positives and false positives.  It is highly advised, to isolate the unsorted
## skeletal elements into a new dataset and re-apply the pair-matching analysis
## on this new smaller subgroup.
##
## As an alternative functionality, @code{longbone_Pair} can work directly with
## 3D models.  @code{@var{sorted} = longbone_Pair (@var{folder}, @var{BONE})}
## will analyze all available bones models of type @var{BONE} which are present
## in the directory specified by @var{folder} and return the output as described
## above.  When @var{folder} is specified, @code{longbone_Pair} utilizes the
## @code{longbone_Geometry} function to retrieve the @var{DATA} from the 3D
## models as well as identify their side.  Keep in mind that in certain cases
## the side of ulna bones might be misidentified, thus it is good practice to
## check the side of each element in the returning tables.  The filenames of the
## 3D models are used as identifiers for the skeletal elements being analyzed.
##
## @seealso{longbone_Sex, longbone_Geometry, longbone_Measurements}
## @end deftypefn

function [varargout] = longbone_Pair (DATA, BONE)

  ## Check number of output arguments
  if (nargout > 3)
    error ("longbone_Pair: too many output arguments.");
  endif

  ## Check for empty input
  if (isempty (DATA))
    error ("longbone_Pair: no DATA or FOLDER.");
  endif

  ## Check bone type
  if (! ismember (BONE, {'Femur', 'Humerus', 'Tibia', 'Ulna'}))
    error ("longbone_Pair: BONE must be a valid bone type.");
  endif

  ## Load pair matching descriptives
  [mu, sd, lb, ub, id] = load_descriptives (BONE, 'pair');

  ## Handle data input
  if (ischar (DATA))
    folder = DATA;
    DATA = [];
    NAME = {};
    samples = 0;
    files = dir (fullfile (folder, "*.obj"));
    for i = 1:length (files)
      [~, ~, bone, extra, data] = longbone_Geometry (folder, files(i).name, {BONE});
      if (strcmpi (bone, BONE))
        if (strcmpi (extra.side, 'left'))
          side = 1;
        else
          side = 2;
        endif
        DATA = [DATA; side, data];
        NAME = [NAME; files(i).name(1:end-4)];
        samples += 1;
      endif
    endfor
  elseif (isnumeric (DATA) && size (DATA, 2) == 62)
    samples = rows (DATA);
    NAME = arrayfun ('num2str', [1:samples]', "UniformOutput", false);
  elseif (iscell (DATA) && size (DATA, 2) == 63)
    samples = rows (DATA);
    NAME = DATA(:,1);
    DATA(:,1) = [];
    DATA = cell2mat (DATA);
  else
    error ("longbone_Pair: invalid DATA.");
  endif

  ## Split DATA into left and right side elements
  LSIDE = DATA(DATA(:,1) == 1, [false, id]);
  LNAME = NAME(DATA(:,1) == 1);
  RSIDE = DATA(DATA(:,1) == 2, [false, id]);
  RNAME = NAME(DATA(:,1) == 2);

  ## Handle extreme case (naive user) only one side is present
  if (isempty (LNAME))
    error ("longbone_Pair: only right-hand side skeletal elements available.");
  elseif (isempty (RNAME))
    error ("longbone_Pair: only left-hand side skeletal elements available.");
  endif

  ## Create a permutation matrix with all possible pairs between left and right
  ## side samples.  First and second columns of the permutation matrix contain
  ## the row (sample) indices to the left and right side variables, respectively.
  TestDATA = [];
  Lvec = abs (LSIDE);
  Rvec = abs (RSIDE);
  R = [1:size(RSIDE, 1)]';
  for i = 1:size (LSIDE, 1)
    L = i * ones (size (R));
    D = Lvec(i,:) - Rvec;
    TestDATA = [TestDATA; L, R, D];
  endfor

  ## Find definite mismatches and exclude them from permutation test matrix
  rejected = [];
  for var = 1:sum (id)
    reject_i = find (TestDATA(:,var+2) < lb(id)(var));
    rejected = [rejected; reject_i];
    reject_i = find (TestDATA(:,var+2) > ub(id)(var));
    rejected = [rejected; reject_i];
  endfor
  rejected = unique (rejected);
  ## Keep remaining pairs as plausible matches
	plausible = TestDATA;
	plausible(rejected,:) = [];

  ## If no plausible matches remain, return the sorted list table
  if (isempty (plausible))
    sorted_elements = samples;
    sorted_pairs = 0;
    sorted_by_elimination = sorted_elements;
    individuals = sorted_elements;
    unsorted_elements = 0;
    plausible_pairs = 0;
    L_side = string ([LNAME; num2cell(nan(size (RNAME)))]);
    R_side = string ([num2cell(nan(size (LNAME))); RNAME; ]);
    scores = nan (samples, 1);
    varnames = {'Left side', 'Right side', 'Score'};
    rownames = arrayfun (@(x) sprintf ("inidividual %d", x), [1:samples]', ...
                         'UniformOutput', false);
    varargout{1} = table (L_side, R_side, scores, 'VariableNames', varnames, ...
                         'RowNames', rownames);
    ## Return stats table if requested
    if (nargout > 1)
      varnames = {'Samples', 'Sorted', 'Paired', 'Single', ...
                  'Individuals', 'Unsorted', 'Plausible'};
      varargout{2} = table (samples, sorted_elements, sorted_pairs, ...
                            sorted_by_elimination, individuals, ...
                            unsorted_elements, plausible_pairs, ...
                            'VariableNames', varnames);
    endif
    ## Return an empty matrix for unsorted list
    if (nargout > 2)
      varnames = {'Left side', 'Right side', 'Score'};
      varargout{3} = table ('Size', [0,3], 'VariableTypes', ...
                           {'string', 'string', 'doublenan'}, ...
                            'VariableNames', varnames);
    endif
    return
  endif

  ## Remove plausible matches from initial left and right side lists
  ## of skeletal elements to identify single-sided samples.
  Lindex = [1:length(LNAME)]';
  Rindex = [1:length(RNAME)]';
  Lindex(ismember (Lindex, unique (plausible(:,1)))) = [];
  Rindex(ismember (Rindex, unique (plausible(:,2)))) = [];
  sorted = [Lindex, nan(length (Lindex), 2); ...
            nan(length (Rindex), 1), Rindex, nan(length (Rindex), 1)];
  ## Keep number of single elements sorted by elimination
  sorted_by_elimination = size (sorted, 1);

  ## Compute sum of absolute z-scores for each paired sample
  ## in the remaining plausible matches
  plausible(:,3) = sum (abs ((plausible(:,[3:end]) - mu(id)) ./ sd(id)), 2);
  plausible(:,[4:end]) = [];

  ## Scan through the remaining plausible cases and cluster the
  ## associated pairs into separate subgroups for each side of bones
  cluster = cluster_pairs (plausible);

  ## Compare the scores of every element with each paired association between
  ## both sides and keep the matching pairs as a sorted pair when both sides
  ## exhibit the lowest score on the same matched pair and the score is
  ## progressively below 30
  [sorted, cluster] = compare_scores_with_threshold (sorted, cluster);

  ## Calculate the number of pairs sorted by lowest score below 30
  sorted_by_threshold = size (sorted, 1) - sorted_by_elimination;

  ## Scan through any remaining plausible matches (if applicable),
  ## recluster plausible pairs, and compare scores by difference
  if (! isempty (cluster))
    s = 0;
    for c = 1:length (cluster)
      s += size (cluster(c).left, 1) + size (cluster(c).right, 1);
    endfor
    if (s > 0)
      ## Recluster the associated pairs into separate subgroups
      cluster = recluster_pairs (cluster);

      ## In each subgroup compare the scores of every element with each paired
      ## association between both sides and keep the matching pairs as a sorted
      ## pair when both sides exhibit the lowest score on the same matched pair
      ## and the score is lower than the second smaller one by at least 5 units
      [sorted, cluster] = compare_scores_with_difference (sorted, cluster);

      ## Calculate the number of pairs sorted by lowest score by difference > 5
      sorted_by_difference = size (sorted,1) - ...
                             (sorted_by_elimination + sorted_by_threshold);
      ## Calculate the total number of pairs sorted
      sorted_pairs = sorted_by_threshold + sorted_by_difference;
    else
      sorted_pairs = sorted_by_threshold;
    endif
  else
    sorted_pairs = sorted_by_threshold;
  endif

  ## Any remaining plausible matches are concatenated into unsorted pairs matrix
  if (! isempty (cluster))
    unsorted = [];
    for c = 1:length (cluster)
      unsorted = [unsorted; cluster(c).left; cluster(c).right];
    endfor
    unsorted = unique (unsorted, "rows");
    unsorted_elements = numel (unique (unsorted(:,1))) + ...
                        numel (unique (unsorted(:,2)));
  else
    unsorted = [];
    unsorted_elements = 0;
  endif
  plausible_pairs = size (unsorted, 1);

  ## Calculate the total number of sorted elements
  ## and distinct individuals identified
  sorted_elements = sum (sum (isfinite (sorted(:,[1:2]))));
  individuals = size (sorted, 1);

  ## Return sorted list table
  L_TF = isfinite (sorted(:,1));
  R_TF = isfinite (sorted(:,2));
  Lindex = LNAME(sorted(:,1)(L_TF));
  Rindex = RNAME(sorted(:,2)(R_TF));
  L_side = num2cell (nan (size (L_TF)));
  L_side(L_TF) = Lindex;
  R_side = num2cell (nan (size (R_TF)));
  R_side(R_TF) = Rindex;
  varnames = {'Left side', 'Right side', 'Score'};
  rownames = arrayfun (@(x) sprintf ("inidividual %d", x), [1:individuals]', ...
                       'UniformOutput', false);
  varargout{1} = table (string (L_side), string (R_side), sorted(:,3), ...
                        'VariableNames', varnames, 'RowNames', rownames);
  ## Return stats table if requested
  if (nargout > 1)
    varnames = {'Samples', 'Sorted', 'Paired', 'Single', ...
                'Individuals', 'Unsorted', 'Plausible'};
    varargout{2} = table (samples, sorted_elements, sorted_pairs, ...
                          sorted_by_elimination, individuals, ...
                          unsorted_elements, plausible_pairs, ...
                          'VariableNames', varnames);
  endif
  ## Return unsorted list table if requested
  if (nargout > 2)
    varnames = {'Left side', 'Right side', 'Score'};
    if (isempty (unsorted))
      varargout{3} = table ('Size', [0,3], 'VariableTypes', ...
                            {'string', 'string', 'doublenan'}, ...
                            'VariableNames', varnames);
    else
      L_side = LNAME(unsorted(:,1));
      R_side = RNAME(unsorted(:,2));
      scores = unsorted(:,3);
      varargout{3} = table (string (L_side), string (R_side), ...
                            scores, 'VariableNames', varnames);
    endif
  endif

endfunction

## Scan through the list of plausible matches and cluster the
## associated pairs into separate subgroups for each side of bones
function	cluster = cluster_pairs (plausible);
  savelist = plausible;
  group = 0;
  while (! isempty (plausible))
    group += 1;
    cluster(group).left = [];
    complete = false;
    ## Find a sample with minimum occurence and use it as a seed
    samples = unique (plausible(:,1));
    clear nsamples;
    for s = 1:length (samples)
      idx = find (plausible(:,1) == samples(s));
      nsamples(s,:) = [length(idx), samples(s)];
    endfor
    nsamples = sortrows (nsamples, 1);
    idx = find (plausible(:,1) == nsamples(1,2));
    left_seed(group) = plausible(idx(1),1);
    right_seed(group) = plausible(idx(1),2);
    left_samples = left_seed(group);
    right_samples = [];
    while (! complete)
      ## Find occurences of right samples according to the left samples
      for i = 1:length (left_samples)
        idx = find (plausible(:,1) == left_samples(i));
        right_samples = [right_samples; plausible(idx,2)];
        cluster(group).left = [cluster(group).left; plausible(idx,:)];
        plausible(idx,:) = [];
      endfor
      k = 0;
      ## Find occurences of left samples according to the right samples
      for i = 1:length (right_samples)
        idx = find (plausible(:,2) == right_samples(i));
        if (! isempty (idx))
          left_samples = [left_samples; plausible(idx,1)];
        else
          k += 1;
        endif
      endfor
      if (i == k)
        complete = true;
      endif
    endwhile
  endwhile
  plausible = savelist;
  group = 0;
  while (! isempty (plausible))
    group += 1;
    cluster(group).right = [];
    complete = false;
    right_samples = right_seed(group);
    left_samples = [];
    while (! complete)
      ## Find occurences of left samples according to the right samples
      for i = 1:length (right_samples)
        idx = find (plausible(:,2) == right_samples(i));
        left_samples = [left_samples; plausible(idx,1)];
        cluster(group).right = [cluster(group).right; plausible(idx,:)];
        plausible(idx,:) = [];
      endfor
      k = 0;
      ## Find occurences of right samples according to the left samples
      for i = 1:length(left_samples)
        idx = find (plausible(:,1) == left_samples(i));
        if (! isempty (idx))
          right_samples = [right_samples; plausible(idx,2)];
        else
          k += 1;
        endif
      endfor
      if (i == k)
        complete = true;
      endif
    endwhile
  endwhile
endfunction

## In each cluster subgroup compare the scores of every element with each paired
## association between both sides and keep the matching pairs as a sorted pair
## when both sides exhibit the lowest score on the same matched pair and the
## score is progressively below 30
function [sorted, cluster] = compare_scores_with_threshold (sorted, cluster);
  for range = 20:1:30
    for c = 1:length (cluster)
      ## Make a list of unique elements
      left_list = unique (cluster(c).left(:,1));
      for s = 1:length (left_list)
        idx = find (cluster(c).left(:,1) == left_list(s));
        if (length (idx) > 0)
          left_side = cluster(c).left(idx,:);
          left_side = sortrows(left_side, 3);
          right_sample = left_side(1,2);
          idx = find (cluster(c).right(:,2) == right_sample);
          right_side = cluster(c).right(idx,:);
          right_side = sortrows (right_side, 3);
          score = right_side(1,3);
          ## If elements in first rows (lowest scores) match (have the
          ## same sample indices in both left and right subgroups) and
          ## their removal does NOT eliminate other samples from the
          ## same cluster, append the pair in the sorted list and remove
          ## their plausible paired instances
          if (score < range && left_side(1,1) == right_side(1,1) &&
                               left_side(1,2) == right_side(1,2))
            ## Left cluster group
            idxL1 = find (cluster(c).left(:,1) == left_side(1,1));
            L2side = cluster(c).left(:,2);
            L2_rem = L2side(idxL1);
            L2_rem(L2_rem == left_side(1,2)) = [];
            if (isempty (L2_rem))
              L2 = true;
            else
              L2 = arrayfun (@(x) sum (L2_rem == x) < sum (L2side == x), L2_rem);
            endif
            idxL2 = find (cluster(c).left(:,2) == left_side(1,2));
            L1side = cluster(c).left(:,1);
            L1_rem = L1side(idxL2);
            L1_rem(L1_rem == left_side(1,1)) = [];
            if (isempty (L1_rem))
              L1 = true;
            else
              L1 = arrayfun (@(x) sum (L1_rem == x) < sum (L1side == x), L1_rem);
            endif
            ## Right cluster group
            idxR1 = find (cluster(c).right(:,1) == right_side(1,1));
            R2side = cluster(c).right(:,2);
            R2_rem = R2side(idxR1);
            R2_rem(R2_rem == right_side(1,2)) = [];
            if (isempty (R2_rem))
              R2 = true;
            else
              R2 = arrayfun (@(x) sum (R2_rem == x) < sum (R2side == x), R2_rem);
            endif
            idxR2 = find (cluster(c).right(:,2) == right_side(1,2));
            R1side = cluster(c).right(:,1);
            R1_rem = R1side(idxR2);
            R1_rem(R1_rem == right_side(1,1)) = [];
            if (isempty (R1_rem))
              R1 = true;
            else
              R1 = arrayfun (@(x) sum (R1_rem == x) < sum (R1side == x), R1_rem);
            endif
            if (all (L1) && all (L2) && all (R1) && all (R2))
              ## Append to sorted list
              sorted = [sorted; left_side(1,:)];
              ## Remove sorted elements from the cluster
              idx = find (cluster(c).left(:,1) == left_side(1,1));
              cluster(c).left(idx,:) = [];
              idx = find (cluster(c).left(:,2) == left_side(1,2));
              cluster(c).left(idx,:) = [];
              idx = find (cluster(c).right(:,1) == right_side(1,1));
              cluster(c).right(idx,:) = [];
              idx = find (cluster(c).right(:,2) == right_side(1,2));
              cluster(c).right(idx,:) = [];
            endif
          endif
        endif
      endfor
    endfor
    ## If any subgroup contains an identical single match, consider it a true
    ## match, append it in the sorted list, and remove it from the cluster
    for c = length (cluster):-1:1
      if (size (cluster(c).left, 1) == 1 && size (cluster(c).right, 1) == 1)
        left_side = cluster(c).left(1,:);
        right_side = cluster(c).right(1,:);
        if (left_side(1,1) == right_side(1,1) &&
            left_side(1,2) == right_side(1,2))
          sorted = [sorted; left_side(1,:)];
          cluster(c) = [];
        endif
      elseif (size (cluster(c).left, 1) == 0 && size (cluster(c).right, 1) == 0)
        cluster(c) = [];
      else
        index = 0;
        for s = 1:size (cluster(c).left, 1)
          left_sample = cluster(c).left(s,1);
          right_sample = cluster(c).left(s,2);
          if (sum (cluster(c).left(:,1) == left_sample) == 1 &&
              sum (cluster(c).left(:,2) == right_sample) == 1)
            index += 1;
            sorted = [sorted; cluster(c).left(s,:)];
            remove(index) = s;
            ## Remove from right cluster as well
            idx = find (cluster(c).right(:,1) == left_sample);
            cluster(c).right(idx,:) = [];
            idx = find (cluster(c).right(:,2) == right_sample);
            cluster(c).right(idx,:) = [];
          endif
        endfor
        if (exist ("remove", "var"))
          cluster(c).left(remove,:) = [];
          clear remove;
        endif
      endif
    endfor
  endfor
endfunction


## Scan through the remaining list of plausible matches and re-cluster
## the associated pairs into separate subgroups for each side of bones
function cluster2 = recluster_pairs (cluster)
  Lgroup = 0;
  Rgroup = 0;
  for c = 1:length (cluster)
    plausible = cluster(c).left;
    while (! isempty (plausible))
      Lgroup += 1;
      cluster2(Lgroup).left = [];
      complete = false;
      ## Find a sample with minimum occurence and use it as a seed
      samples = unique (plausible(:,1));
      clear nsamples;
      for s = 1:length (samples)
        idx = find (plausible(:,1) == samples(s));
        nsamples(s,:) = [length(idx), samples(s)];
      endfor
      nsamples = sortrows (nsamples, 1);
      idx = find (plausible(:,1) == nsamples(1,2));
      left_seed(Lgroup) = plausible(idx(1),1);
      right_seed(Lgroup) = plausible(idx(1),2);
      left_samples = left_seed(Lgroup);
      right_samples = [];
      while (! complete)
        ## Find occurences of right samples according to the left samples
        for i = 1:length (left_samples)
          idx = find (plausible(:,1) == left_samples(i));
          right_samples = [right_samples; plausible(idx,2)];
          cluster2(Lgroup).left = [cluster2(Lgroup).left; plausible(idx,:)];
          plausible(idx,:) = [];
        endfor
        k = 0;
        ## Find occurences of left samples according to the right samples
        for i = 1:length (right_samples)
          idx = find (plausible(:,2) == right_samples(i));
          if (! isempty (idx))
            left_samples = [left_samples; plausible(idx,1)];
          else
            k += 1;
          endif
        endfor
        if (i == k)
          complete = true;
        endif
      endwhile
    endwhile
    plausible = cluster(c).right;
    while (! isempty (plausible))
      Rgroup += 1;
      cluster2(Rgroup).right = [];
      complete = false;
      right_samples = right_seed(Rgroup);
      left_samples = [];
      while (! complete)
        ## Find occurences of left samples according to the right samples
        for i = 1:length (right_samples)
          idx = find (plausible(:,2) == right_samples(i));
          left_samples = [left_samples; plausible(idx,1)];
          cluster2(Rgroup).right = [cluster2(Rgroup).right; plausible(idx,:)];
          plausible(idx,:) = [];
        endfor
        k = 0;
        ## Find occurences of right samples according to the left samples
        for i = 1:length (left_samples)
          idx = find (plausible(:,1) == left_samples(i));
          if (! isempty (idx))
            right_samples = [right_samples; plausible(idx,2)];
          else
            k += 1;
          endif
        endfor
        if (i == k)
          complete = true;
        endif
      endwhile
    endwhile
  endfor
endfunction

## In each cluster subgroup compare the scores of every element with each paired
## association between both sides and keep the matching pairs as a sorted pair
## when both sides exhibit the lowest score on the same matched pair and the
## score is lower than the second smaller one by at least 5 units
function [sorted, cluster] = compare_scores_with_difference (sorted, cluster)
  for iter = 1:2
    for c = length (cluster):-1:1
      ## Make a list of unique elements
      left_list = unique (cluster(c).left(:,1));
      for s = 1:length (left_list)
        idxL = find (cluster(c).left(:,1) == left_list(s));
        idxR = find (cluster(c).right(:,2) == left_list(s));
        if (length (idxL) > 0 && length (idxR) > 0)
          left_side = cluster(c).left(idxL,:);
          right_side = cluster(c).right(idxR,:);
          ## If multiple pairs are present on right side only, use it explicitly
          if (size (right_side, 1) > 1 && size (left_side, 1) == 1)
            left_side = sortrows (left_side, 3);
            right_side = sortrows (right_side, 3);
            score_diff = abs (right_side(1,3) - right_side(2,3));
            ## If elements in first rows (lowest scores) match (have the
            ## same sample indices in both left and right subgroups) and
            ## their removal does NOT eliminate other samples from the
            ## same cluster, append the pair in the sorted list and remove
            ## their plausible paired instances
            if (score_diff > 5 && left_side(1,1) == right_side(1,1)
                               && left_side(1,2) == right_side(1,2))
              ## Left cluster group
              idxL1 = find (cluster(c).left(:,1) == left_side(1,1));
              L2side = cluster(c).left(:,2);
              L2_rem = L2side(idxL1);
              L2_rem(L2_rem == left_side(1,2)) = [];
              if (isempty (L2_rem))
                L2 = true;
              else
                L2 = arrayfun (@(x) sum (L2_rem == x) < sum (L2side == x), L2_rem);
              endif
              idxL2 = find (cluster(c).left(:,2) == left_side(1,2));
              L1side = cluster(c).left(:,1);
              L1_rem = L1side(idxL2);
              L1_rem(L1_rem == left_side(1,1)) = [];
              if (isempty (L1_rem))
                L1 = true;
              else
                L1 = arrayfun (@(x) sum (L1_rem == x) < sum (L1side == x), L1_rem);
              endif
              ## Right cluster group
              idxR1 = find (cluster(c).right(:,1) == right_side(1,1));
              R2side = cluster(c).right(:,2);
              R2_rem = R2side(idxR1);
              R2_rem(R2_rem == right_side(1,2)) = [];
              if (isempty (R2_rem))
                R2 = true;
              else
                R2 = arrayfun (@(x) sum (R2_rem == x) < sum (R2side == x), R2_rem);
              endif
              idxR2 = find (cluster(c).right(:,2) == right_side(1,2));
              R1side = cluster(c).right(:,1);
              R1_rem = R1side(idxR2);
              R1_rem(R1_rem == right_side(1,1)) = [];
              if (isempty (R1_rem))
                R1 = true;
              else
                R1 = arrayfun (@(x) sum (R1_rem == x) < sum (R1side == x), R1_rem);
              endif
              if (all (L1) && all (L2) && all (R1) && all (R2))
                ## Append to sorted list
                sorted = [sorted; left_side(1,:)];
                ## Remove sorted elements from the cluster
                idx = find (cluster(c).left(:,1) == left_side(1,1));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).left(:,2) == left_side(1,2));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).right(:,1) == right_side(1,1));
                cluster(c).right(idx,:) = [];
                idx = find (cluster(c).right(:,2) == right_side(1,2));
                cluster(c).right(idx,:) = [];
              endif
            endif
          ## If multiple pairs are present on left side only, use it explicitly
          elseif (size (left_side, 1) > 1 && size (right_side, 1) == 1)
            left_side = sortrows (left_side, 3);
            right_side = sortrows (right_side, 3);
            score_diff = abs (left_side(1,3) - left_side(2,3));
            ## If elements in first rows (lowest scores) match (have the
            ## same sample indices in both left and right subgroups) and
            ## their removal does NOT eliminate other samples from the
            ## same cluster, append the pair in the sorted list and remove
            ## their plausible paired instances
            if (score_diff > 5 && left_side(1,1) == right_side(1,1)
                               && left_side(1,2) == right_side(1,2))
              ## Left cluster group
              idxL1 = find (cluster(c).left(:,1) == left_side(1,1));
              L2side = cluster(c).left(:,2);
              L2_rem = L2side(idxL1);
              L2_rem(L2_rem == left_side(1,2)) = [];
              if (isempty (L2_rem))
                L2 = true;
              else
                L2 = arrayfun (@(x) sum (L2_rem == x) < sum (L2side == x), L2_rem);
              endif
              idxL2 = find (cluster(c).left(:,2) == left_side(1,2));
              L1side = cluster(c).left(:,1);
              L1_rem = L1side(idxL2);
              L1_rem(L1_rem == left_side(1,1)) = [];
              if (isempty (L1_rem))
                L1 = true;
              else
                L1 = arrayfun (@(x) sum (L1_rem == x) < sum (L1side == x), L1_rem);
              endif
              ## Right cluster group
              idxR1 = find (cluster(c).right(:,1) == right_side(1,1));
              R2side = cluster(c).right(:,2);
              R2_rem = R2side(idxR1);
              R2_rem(R2_rem == right_side(1,2)) = [];
              if (isempty (R2_rem))
                R2 = true;
              else
                R2 = arrayfun (@(x) sum (R2_rem == x) < sum (R2side == x), R2_rem);
              endif
              idxR2 = find (cluster(c).right(:,2) == right_side(1,2));
              R1side = cluster(c).right(:,1);
              R1_rem = R1side(idxR2);
              R1_rem(R1_rem == right_side(1,1)) = [];
              if (isempty (R1_rem))
                R1 = true;
              else
                R1 = arrayfun (@(x) sum (R1_rem == x) < sum (R1side == x), R1_rem);
              endif
              if (all (L1) && all (L2) && all (R1) && all (R2))
                ## Append to sorted list
                sorted = [sorted; left_side(1,:)];
                ## Remove sorted elements from the cluster
                idx = find (cluster(c).left(:,1) == left_side(1,1));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).left(:,2) == left_side(1,2));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).right(:,2) == right_side(1,2));
                cluster(c).right(idx,:) = [];
                idx = find (cluster(c).right(:,1) == right_side(1,1));
                cluster(c).right(idx,:) = [];
              endif
            endif
          ## If multiple pairs are present on both sides,
          ## find the minimum difference from either side
          elseif (size (left_side, 1) > 1 && size (right_side, 1) > 1)
            left_side = sortrows (left_side, 3);
            right_side = sortrows (right_side, 3);
            score_L = abs (left_side(1,3) - left_side(2,3));
            score_R = abs (right_side(1,3) - right_side(2,3));
            score_diff = min ([score_L, score_R]);
            ## If elements in first rows (lowest scores) match (have the
            ## same sample indices in both left and right subgroups) and
            ## their removal does NOT eliminate other samples from the
            ## same cluster, append the pair in the sorted list and remove
            ## their plausible paired instances
            if (score_diff > 5 && left_side(1,1) == right_side(1,1)
                               && left_side(1,2) == right_side(1,2))
              ## Left cluster group
              idxL1 = find (cluster(c).left(:,1) == left_side(1,1));
              L2side = cluster(c).left(:,2);
              L2_rem = L2side(idxL1);
              L2_rem(L2_rem == left_side(1,2)) = [];
              if (isempty (L2_rem))
                L2 = true;
              else
                L2 = arrayfun (@(x) sum (L2_rem == x) < sum (L2side == x), L2_rem);
              endif
              idxL2 = find (cluster(c).left(:,2) == left_side(1,2));
              L1side = cluster(c).left(:,1);
              L1_rem = L1side(idxL2);
              L1_rem(L1_rem == left_side(1,1)) = [];
              if (isempty (L1_rem))
                L1 = true;
              else
                L1 = arrayfun (@(x) sum (L1_rem == x) < sum (L1side == x), L1_rem);
              endif
              ## Right cluster group
              idxR1 = find (cluster(c).right(:,1) == right_side(1,1));
              R2side = cluster(c).right(:,2);
              R2_rem = R2side(idxR1);
              R2_rem(R2_rem == right_side(1,2)) = [];
              if (isempty (R2_rem))
                R2 = true;
              else
                R2 = arrayfun (@(x) sum (R2_rem == x) < sum (R2side == x), R2_rem);
              endif
              idxR2 = find (cluster(c).right(:,2) == right_side(1,2));
              R1side = cluster(c).right(:,1);
              R1_rem = R1side(idxR2);
              R1_rem(R1_rem == right_side(1,1)) = [];
              if (isempty (R1_rem))
                R1 = true;
              else
                R1 = arrayfun (@(x) sum (R1_rem == x) < sum (R1side == x), R1_rem);
              endif
              if (all (L1) && all (L2) && all (R1) && all (R2))
                ## Append to sorted list
                sorted = [sorted; left_side(1,:)];
                ## Remove sorted elements from the cluster
                idx = find (cluster(c).left(:,1) == left_side(1,1));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).left(:,2) == left_side(1,2));
                cluster(c).left(idx,:) = [];
                idx = find (cluster(c).right(:,2) == right_side(1,2));
                cluster(c).right(idx,:) = [];
                idx = find (cluster(c).right(:,1) == right_side(1,1));
                cluster(c).right(idx,:) = [];
              endif
            endif
          else
            sorted = [sorted; left_side(1,:)];
            ## Remove sorted elements from the cluster
            idx = find (cluster(c).left(:,1) == left_side(1,1));
            cluster(c).left(idx,:) = [];
            idx = find (cluster(c).left(:,2) == left_side(1,2));
            cluster(c).left(idx,:) = [];
            idx = find (cluster(c).right(:,2) == right_side(1,2));
            cluster(c).right(idx,:) = [];
            idx = find (cluster(c).right(:,1) == right_side(1,1));
            cluster(c).right(idx,:) = [];
          endif
        endif
      endfor
    endfor
  endfor
endfunction


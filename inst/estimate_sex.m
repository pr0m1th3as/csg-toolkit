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
## @deftypefn  {csg-toolkit} {} estimate_sex
## @deftypefnx {csg-toolkit} {} estimate_sex (@var{filename})
## @deftypefnx {csg-toolkit} {} estimate_sex (@dots{}, @qcode{"classifier"}, @var{fname})
## @deftypefnx {csg-toolkit} {} estimate_sex (@dots{}, @qcode{"bone"}, @var{type})
##
## Estimate biological sex from longbone CSG properties.
##
## @code{estimate_sex} opens a graphical user interface to facilitate biological
## sex estimation from the longbones of the upper and lower limbs.  It loads the
## default classifier from the file @qcode{Athens_Collection_sex_classifier.mat}
## which was based on a modern Greek population sample, derived from the Athens
## Collection.  The GUI allows the user to either select a CSV file containing
## the CSG properties or an OBJ file contating the 3D bone model.  The longbones
## currently supported are the femur, tibia, and humerus.
##
## @seealso{inspect_CSG, longbone_Geometry, match_pairs}
## @end deftypefn

function estimate_sex (varargin)

  ## Check input arguments
  if (nargin == 0)
    classifier = load ("Athens_Collection_sex_classifier.mat");
    sample_IDs = "";
    CSGdata = [];
    bone = "";
    args = {};
  endif

  ## Parse optional Name-Value paired arguments
  if (nargin > 1)
    optNames = {"classifier", "bone"};
    dfValues = {{"Athens_Collection_sex_classifier.mat"}, {""}};
    [class_fname, bone, args] = pairedArgs (optNames, dfValues, varargin(:));
    if (strcmpi (class_fname([end-3:end]), ".mat"))
      classifier = load (class_fname);
    else
      error ("estimate_sex: invalid classifier file type.");
    endif
    if (! isfield (classifier, "classifiers") ||
        ! isfield (classifier, "description"))
      error ("estimate_sex: '%s' is not a valid classifier container.", ...
             class_fname);
    endif
    if (! any (strcmpi (bone, {"Femur", "Tibia", "Humerus", "Ulna"})))
      error ("estimate_sex: invalid type of bone.");
    endif
  else
    classifier = load ("Athens_Collection_sex_classifier.mat");
    args = varargin(:);
  endif

  ## Get filename from first input argument
  if (numel (args) > 0)
    if (! ischar (args{1}))
      error ("estimate_sex: 3D model or CSG data filename must be a string.");
    elseif (exist (args{1}, "file") != 2)
      error ("estimate_sex: '%s' does not exist.", args{1});
    endif
    ## 3D model filename
    if (strcmpi (args{1}([end-3:end]), ".obj"))
      [sample_IDs, CSGdata, bone] = get_OBJ_data (args{1});
    ## CSV data filename
    elseif (strcmpi (args{1}([end-3:end]), ".csv"))
      [sample_IDs, CSGdata] = get_CSV_data (args{1});
      bone = "";
    else
      error ("estimate_sex: invalid data file type.");
    endif
  endif

  ## Call main ui for sex estimation based on CSG data
  h = user_interface_CSG (classifier, sample_IDs, CSGdata, bone);

endfunction

## GUI for sex estimation based on CSG data
function h = user_interface_CSG (classifier, sample_IDs, CSGdata, bone)
  ## List bones and sides for the CSG-Toolkit sex classifier
  if (isempty (bone))
    bone = {"Femur", "Tibia", "Humerus"};
  endif
  side = {"Left", "Right"};
  ## Set UI background colour
  bg_color = [0.8, 0.8, 0.8];
  ## Create a main window for user interface and
  ## populate it with necessary user options
  h.f = figure ("Name", ...
                "skeletal sex estimation based on longbone CSG properties", ...
                "numbertitle", "off", "menubar", "none", "color", bg_color, ...
                "Position", [300, 300, 560, 450], "resize", "off");
  ## Save data into figure
  guidata (h.f, CSGdata);
  ## display available classification methods (RBF or LDA)
  class_methods = {"RBF kernel SVM classification"; ...
                   "Linear Discriminat Function Analysis"};
  h.class_method = uibuttongroup ("title", "Available classification methods",...
                                  "titleposition", "centertop", ...
                                  "backgroundcolor", bg_color,...
                                  "position", [0.2, 0.76, 0.6, 0.2]);
  uicontrol (h.class_method, "style", "radiobutton", ...
                             "string", class_methods{1}, "tag", "RBF",...
                             "backgroundcolor", bg_color, ...
                             "position", [50, 40, 300, 30]);
  uicontrol (h.class_method, "style", "radiobutton", ...
                             "string", class_methods{2}, "tag", "LDA",...
                             "backgroundcolor", bg_color, ...
                             "position", [50, 10, 300, 30]);
  set (get (h.class_method, "children")(2), "selected", "on", "value", 1);
  ## selection of skeletal sample, bone and side (if applicable)
  if (! isempty (sample_IDs))
    h.text_S = uicontrol (h.f, "style", "text", "string", "Available samples",...
                               "backgroundcolor", bg_color, ...
                               "position", [30, 300, 200, 30]);
    h.sample = uicontrol (h.f, "style", "popupmenu", "string", sample_IDs, ...
                               "position", [30, 270, 200, 30]);
  else
    h.text_S = uicontrol (h.f, "style", "text", "string", "Select file", ...
                               "backgroundcolor", bg_color, ...
                               "position", [30, 300, 200, 30]);
    h.sample = uicontrol (h.f, "style", "popupmenu", "string", sample_IDs, ...
                               "position", [30, 270, 200, 30]);
    loaddata = uicontrol (h.f, "style", "pushbutton", "tag", "loadbutton", ...
                               "string", "load 3D model or CSG data",...
                               "position", [30, 270, 200, 30], ...
                               "backgroundcolor", bg_color);
  endif
  if (isempty (bone))
    bone = {"Femur", "Tibia", "Humerus"};
    h.text_B = uicontrol (h.f, "style", "text", "string", "Select bone", ...
                               "backgroundcolor", bg_color, ...
                               "position", [280, 300, 100, 30]);
    h.bone_type = uicontrol (h.f, "style", "popupmenu", "string", bone, ...
                                  "position", [280, 270, 100, 30]);
  else
    h.text_B = uicontrol (h.f, "style", "text", "string", "Identified bone", ...
                               "backgroundcolor", bg_color, ...
                               "position", [280, 300, 100, 30]);
    h.bone_type = uicontrol (h.f, "style", "popupmenu", "string", bone, ...
                                  "position", [280, 270, 100, 30]);
  endif
  h.text_S = uicontrol (h.f, "style", "text", "string", "Select side", ...
                             "backgroundcolor", bg_color, ...
                             "position", [430, 300, 100, 30]);
  h.bone_side = uicontrol (h.f, "style", "popupmenu", "string", side, ...
                                "position", [430, 270, 100, 30]);
  ## selection of particular classifiers
  h.class1 = uicontrol (h.f, "style", "checkbox", "string", "Classifier #1:",...
                             "fontangle", "italic", "fontweight", "bold", ...
                             "backgroundcolor", bg_color, ...
                             "position", [30, 210, 120, 20]);
  h.class2 = uicontrol (h.f, "style", "checkbox", "string", "Classifier #2:",...
                             "fontangle", "italic", "fontweight", "bold", ...
                             "backgroundcolor", bg_color,...
                             "position", [30, 160, 120, 20]);
  h.class3 = uicontrol (h.f, "style", "checkbox", "string", "Classifier #3:",...
                             "fontangle", "italic", "fontweight", "bold", ...
                             "backgroundcolor", bg_color,...
                             "position", [30, 110, 120, 20]);
  set (h.class1, "value", 1, "selected", "on");
  ## containers for each classifier's results
  results_string = {"","",""};
  h.result1 = uicontrol(h.f, "style", "text", "string", results_string{1}, ...
                             "horizontalalignment", "left", ...
                             "backgroundcolor", bg_color, ...
                             "position", [160, 210, 400, 20]);
  h.result2 = uicontrol(h.f, "style", "text", "string", results_string{2}, ...
                             "horizontalalignment", "left", ...
                             "backgroundcolor", bg_color, ...
                             "position", [160, 160, 400, 20]);
  h.result3 = uicontrol(h.f, "style", "text", "string", results_string{3}, ...
                             "horizontalalignment", "left", ...
                             "backgroundcolor", bg_color, ...
                             "position", [160, 110, 400, 20]);
  ## push buttons for estimating sex and saving the results in a csv file
  estimate = uicontrol (h.f, "style", "pushbutton", "string", "estimate sex",...
                             "Position", [130, 20, 100, 50], ...
                             "backgroundcolor", bg_color);
  savedata = uicontrol (h.f, "style", "pushbutton", "string", "save results",...
                             "Position", [340, 20, 100, 50], ...
                             "backgroundcolor", bg_color);
  ## set callbacks
  set (loaddata, "callback", {@run_loaddata_CSG, loaddata, h, classifier});
  set (estimate, "callback", {@run_estimate_CSG, estimate, h, classifier});
  set (savedata, "callback", {@run_savedata_CSG, savedata, h, classifier});
endfunction

## Callback function for load 3D model or CSG data from file
function run_loaddata_CSG (loaddata, dummy1, dummy2, h, classifier)
  ## Select file
  filename = uigetfile ({"*.csv", "CSG data"; "*.obj", "3D model"}, ...
                        ["Select a 3D model or a CSG dataset", ...
                         " created created with 'inspect_CSG'"]);
  ## 3D model filename
  if (strcmpi (filename([end-3:end]), ".obj"))
    [sample_IDs, CSGdata, bone] = get_OBJ_data (filename);
    ## Update UI
    h.text_B = uicontrol (h.f, "style", "text", "string", "Identified bone", ...
                               "backgroundcolor", [0.8, 0.8, 0.8], ...
                               "position", [280, 300, 100, 30]);
    h.bone_type = uicontrol (h.f, "style", "popupmenu", "string", bone, ...
                                  "position", [280, 270, 100, 30]);
    h.text_S = uicontrol (h.f, "style", "text", ...
                               "string", "Loaded 3D model",...
                               "backgroundcolor", [0.8, 0.8, 0.8], ...
                               "position", [30, 300, 200, 30]);
    set (h.sample, "string", sample_IDs);
    guidata (h.f, CSGdata);
    delete (loaddata);
  ## CSV data filename
  elseif (strcmpi (filename([end-3:end]), ".csv"))
    [sample_IDs, CSGdata] = get_CSV_data (filename);
    if (! isempty (sample_IDs))
      ## Update UI
      h.text_S = uicontrol (h.f, "style", "text", ...
                                 "string", "Available samples",...
                                 "backgroundcolor", [0.8, 0.8, 0.8], ...
                                 "position", [30, 300, 200, 30]);
      set (h.sample, "string", sample_IDs);
      guidata (h.f, CSGdata);
      delete (loaddata);
    endif
  endif
endfunction

## Callback function for estimating sex from CSG data
function run_estimate_CSG (estimate, dummy1, dummy2, h, classifier)
  ## Update current user selections
  [method, s_idx, sample_ID, bone, side, skelement, DFs] = update_CSG_UI (h);
  ## Get sample from CSGdata
  CSGdata = guidata (h.f);
  sample_data = CSGdata(s_idx,:);
  ## For each DF keep the appropriate variables, normalize them and
  ## estimate sex and posterior probability according to the chosen method
  [sample_sex, sample_pb] = predict_sex_CSG (DFs, method, skelement, ...
                                             classifier, sample_data);
  ## Reset previous results (if any)
  for i = 1:3
    set (h.(sprintf ("result%i", i)), "string", "");
  endfor
  ## Construct and display results on GUI for each selected DF
  r_idx = 0;
  for c_idx = DFs
    r_idx += 1;
    result_str = sprintf ("Sample is %s with posterior probability %0.2f",...
                           sample_sex{r_idx}, sample_pb(r_idx));
    r_sel = cstrcat ("result", sprintf ("%i", c_idx));
    set (h.(r_sel), "string", result_str);
  endfor
endfunction

## Callback function for saving results from CSG data
function run_savedata_CSG (savedata, dummy1, dummy2, h, classifier)
  ## take current user selections and create indexing accordingly
  [method, s_idx, sample_ID, bone, side, skelement, DFs] = update_CSG_UI (h);
  ## get sample from CSGdata
  CSGdata = guidata (h.f);
  sample_data = CSGdata(s_idx,:);
  ## Compute sex estimation
  [sample_sex, sample_pb] = predict_sex_CSG (DFs, method, skelement, ...
                                             classifier, sample_data);
  ## create a header and a standard filename for the results' file
  results_header = {"Sample ID", "bone", "side", "method", "classifier",...
                    "predicted sex", "posterior probability"};
  ## Select file
  filename = uiputfile ({"*.csv", "CSV"}, ...
                        "Select filename for storing classifcation results");
  ## Check if file was actually selected
  if (filename != 0)
    ## Check if file already exists (if not, create one with header only)
    if (! exist (filename))
      cell2csv (filename, results_header);
    else
      ## Check file and rewrite it if it is not an estimate_sex results file
      results = csv2cell (filename);
      if (size (results, 2) != 7)
        results = results_header;
      endif
    endif
    ## Append new entries at the end of CSV file
    r_idx = 0;
    for c_idx = DFs
      r_idx += 1;
      results(end+1,:) = {sample_ID, bone, side, method, c_idx, ...
                          sample_sex{r_idx}, sample_pb(r_idx)};
    endfor
    cell2csv (filename, results);
  endif
endfunction

## Estimate sex
function [sample_sex, sample_pb] = predict_sex_CSG (DFs, method, skelement, ...
                                                    classifier, sample_data)
  if (length (DFs) > 0)
    ## For each selected DF keep the appropriate variables,
    ## normalize them and estimate sex and posterior probability
    ## according to the selected method regarding CSG data
    for i=1:length (DFs)
      % get index of classifier from description table
      temp_dt = classifier.description.(method);
      classifier_idx = cell2mat (temp_dt(find (strcmp (skelement, ...
                                               temp_dt(:,1))), DFs(i) + 1));
      clear temp_dt;
      ## Get normalization coefficients and variable indices
      ## for corresponding selected classifier
      variable_indices = classifier.normalcoef(classifier_idx).var_num;
      norm_mean = classifier.normalcoef(classifier_idx).mean;
      norm_StD = classifier.normalcoef(classifier_idx).StD;
      ## Truncate to required variables and normalize data
      sample_trunc_data = sample_data(variable_indices);
      sample_norm_data = (sample_trunc_data - norm_mean) ./ norm_StD;
      ## Get appropriate classifier
      sel_class = classifier.classifiers(classifier_idx).(method);
      ## Evaluate sample's data with classifier
      df(i) = evaluate_model(sample_norm_data, sel_class);
      ## some cleanup
      clear variable_indices sample_norm_data
      clear sample_trunc_data norm_mean norm_StD sel_class
      ## Evaluate sample's posterior probability of the calculated outcome
      temp_pb = classifier.post_proba(classifier_idx).(method);
      if (isfield (temp_pb, "sectioning_point"))    # for LDA classifier
        sp = temp_pb.sectioning_point;
        cf = temp_pb.centroid_females;
        cm = temp_pb.centroid_males;
        ## Find distance from sectioning point
        dist = df(i) - sp;
        if (cf < cm)      # female centroid smaller than male centroid
          if (dist < 0)   # then female
            sample_sex(i) = {"female"};
            ## Calculate posterior probability
            sample_pb(i) = exp (-dist) / (exp (-dist) + exp (dist));
          else            # then male
            sample_sex(i) = {"male"};
            ## Calculate posterior probability
            sample_pb(i) = exp (dist) / (exp (dist) + exp (-dist));
          endif
        else              # male centroid smaller than female centroid
          if (dist < 0)   # then male
            sample_sex(i) = {"male"};
            ## Calculate posterior probability
            sample_pb(i) = exp (-dist) / (exp (-dist) + exp (dist));
          else            # then female
            sample_sex(i) = {"female"};
            ## Calculate posterior probability
            sample_pb(i) = exp (dist) / (exp (dist) + exp (-dist));
          endif
        endif
        ## some cleanup
        clear sp cf cm
      elseif isfield(temp_pb, "discrete_PDF")       # for RBF classifier
        dPDF = temp_pb.discrete_PDF;
        fg = temp_pb.female_group;
        mg = temp_pb.male_group;
        ## get posterior probability according to
        ## the ranges provided in the lookup table
        sample_pb(i) = dPDF(find (dPDF(:,1) <= abs (df(i)) & ...
                       dPDF(:,2) > abs (df(i))),3);
        ## Check if first group (negative RBF scores) is males or females
        if (mg < fg)        # male group takes negative scores
          if (df(i) < 0)    # then male
            sample_sex(i) = {"male"};
          else              # then female
            sample_sex(i) = {"female"};
          endif
        else                # female group takes negative scores
          if (df(i) < 0)     # then female
            sample_sex(i) = {"female"};
          else              # then male
            sample_sex(i) = {"male"};
          endif
        endif
        ## some cleanup
        clear dPDF fg mg
      endif
      ## some cleanup
      clear temp_pb
    endfor
  endif
endfunction

## Get current user selections from GUI
function [method, sample_idx, sample_ID, ...
          bone, side, skelement, DFs] = update_CSG_UI (h)
  ## Get classification method
  method = get(get(h.class_method, "selectedobject"), "tag");
  ## Get selected sample
  sample_str = get(h.sample, "string");
  sample_idx = get(h.sample, "value");
  if (iscellstr (sample_str))
    sample_ID = sample_str{sample_idx};
  else
    sample_ID = sample_str;
  endif
  ## Get selected bone
  bone_str = get(h.bone_type, "string");
  bone_idx = get(h.bone_type, "value");
  bone = bone_str{bone_idx};
  ## Get selected side
  side_str = get(h.bone_side, "string");
  side_idx = get(h.bone_side, "value");
  side = side_str{side_idx};
  ## Concatenate bone and side for searching the index table
  ## in the classifier's description
  skelement = cstrcat (bone, " ", side);
  ## Get selected classifiers
  DF1 = get(h.class1, "value");
  DF2 = get(h.class2, "value");
  DF3 = get(h.class3, "value");
  DFs = [1, 2, 3];
  DFs = DFs(find (DFs .* [DF1, DF2, DF3]));
endfunction

## Evaluate LDA and RBF kernel based SVM models
function df = evaluate_model (data, classifier)
	## determine the type of classifier by the number of fields
	if (length (fieldnames (classifier)) == 4)
		## RBF
		df = 0;
		for i=1:length (classifier.dual_coef)
			d(i) = norm(classifier.support_vector(i,:) - data);
			t(i) = exp(-classifier.gamma_param * (d(i) ^ 2));
			df += classifier.dual_coef(i) * t(i);
		endfor
		df += classifier.rho_intercept;
	elseif (length (fieldnames (classifier))==1)
		## LDA (constant + coefficients x measurements)
		df = classifier.weights(1) + sum (classifier.weights([2:end]) .* data);
	endif
endfunction

## Retrieve CSG data from CSV file
function [sample_IDs, CSGdata] = get_CSV_data (filename)
  ## Read sample IDs
  all_data = csv2cell (filename);
  ## Check for valid CSV
  if (size(all_data,1) < 2 || size(all_data,2) ~= 47)
    error ("estimate_sex: selected CSV file does not contain CSG properties.");
  endif
  sample_IDs = all_data([2:end],1);
  ## Read necessary CSG data
  CSG = csvread (filename);
  ## Keep maxDistance, Areas, Perimeters, Ix, Iy, Imin, Imax
  CSG(1,:) = [];
  CSG(:,[1,3:7,12,15,20,23,28,31,36,39,44,47]) = [];
  ## Calculate extra variables (Ix/Iy, Imax/Imin, ArPerIndex) used
  ## for sex estimation and rearrange the matrix accordingly
  CSGdata = CSG(:,[1:3]);
  CSGdata = [CSGdata, ((CSG(:,2))*(4*pi))./(CSG(:,3).^2), ...
             CSG(:,[4:5]), CSG(:,4)./CSG(:,5)];
  CSGdata = [CSGdata, CSG(:,[6:7]), CSG(:,7)./CSG(:,6), CSG(:,[8:9])];
  CSGdata = [CSGdata, ((CSG(:,8))*(4*pi))./(CSG(:,9).^2), ...
             CSG(:,[10:11]), CSG(:,10)./CSG(:,11)];
  CSGdata = [CSGdata, CSG(:,[12:13]), CSG(:,13)./CSG(:,12), CSG(:,[14:15])];
  CSGdata = [CSGdata, ((CSG(:,14))*(4*pi))./(CSG(:,15).^2), ...
             CSG(:,[16:17]), CSG(:,16)./CSG(:,17)];
  CSGdata = [CSGdata, CSG(:,[18:19]), CSG(:,19)./CSG(:,18), CSG(:,[20:21])];
  CSGdata = [CSGdata, ((CSG(:,20))*(4*pi))./(CSG(:,21).^2), ...
             CSG(:,[22:23]), CSG(:,22)./CSG(:,23)];
  CSGdata = [CSGdata, CSG(:,[24:25]), CSG(:,25)./CSG(:,24), CSG(:,[26:27])];
  CSGdata = [CSGdata, ((CSG(:,26))*(4*pi))./(CSG(:,27).^2), ...
             CSG(:,[28:29]), CSG(:,28)./CSG(:,29)];
  CSGdata = [CSGdata, CSG(:,[30:31]), CSG(:,31)./CSG(:,30)];
endfunction

## Retrieve CSG data from OBJ file
function [sample_IDs, CSGdata, bone] = get_OBJ_data (filename)
  ## Get sample ID from filename
  sample_IDs = filename([1:end-4]);
  ## Retrieve CSG properties from 3D model
  [GEOM, SMoA, bone] = longbone_Geometry (filename);
  ## Get Max Distance
  v = readObj (filename);
  CSGdata = longbone_maxDistance (v);
  ## Start appending required values to CSG data vector
  for i = 1:5
    CSGdata = [CSGdata, GEOM(i).Area, GEOM(i).Perimeter];
    ArPerIndex = (GEOM(i).Area * (4 * pi)) ./ (GEOM(i).Perimeter .^ 2);
    CSGdata = [CSGdata, ArPerIndex];
    Ix_Iy = SMoA(i).Ix / SMoA(i).Iy;
    CSGdata = [CSGdata, SMoA(i).Ix, SMoA(i).Iy, Ix_Iy];
    CSGdata = [CSGdata, SMoA(i).Imin, SMoA(i).Imax];
    Imax_Imin = SMoA(i).Imin / SMoA(i).Imax;
    CSGdata = [CSGdata, Imax_Imin];
  endfor
endfunction

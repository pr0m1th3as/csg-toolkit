## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef LongBone

  properties (GetAccess = public)
    Vertices = [];
    Faces = [];
    Colors = [];
    TextureCoords = [];
    TextureFaces = [];
    Normals = [];
    FaceNormals = []
    MaterialLib = [];
    TextureImage = [];
    obj_filename = "";
    mtl_filename = "";
    bone = "";
    side = "";
  endproperties

  methods (Hidden)

    function this = LongBone (varargin)
      if (nargin > 0)
        this.Vertices = varargin{1};
      endif
      if (nargin > 1)
        this.Faces = varargin{2};
      endif
      if (nargin > 2)
        this.Colors = varargin{3};
      endif
      if (nargin > 3)
        this.TextureCoords = varargin{4};
      endif
      if (nargin > 4)
        this.TextureFaces = varargin{5};
      endif
      if (nargin > 5)
        this.Normals = varargin{6};
      endif
      if (nargin > 6)
        this.FaceNormals = varargin{7};
      endif
      if (nargin > 7)
        this.MaterialLib = varargin{8};
      endif
      if (nargin > 8)
        this.TextureImage = varargin{9};
      endif
      if (nargin > 9)
        this.obj_filename = varargin{10};
      endif
      if (nargin > 10)
        this.mtl_filename = varargin{11};
      endif
    endfunction

  endmethods

  methods (Access = public)

    function LBM = identifyBone (this)
      ## Identify bone
      bone = longbone_Registration (this.Vertices, this.Faces);
      if (strcmp (bone, "Undetermined"))
        name = inputname (1);
        warning ("identifyBone: 3D mesh in %s is not a valid long bone.", name);
      endif
      ## Append bone info
      LBM = this;
      LBM.bone = bone;
    endfunction

    function LBM = anatomicalPosition (this)
      LBM = this;
      if (isempty (this.Normals))
        [V, F] = longbone_AnatomicalPosition (this.Vertices, this.Faces);
        LBM.Vertices = V;
        LBM.Faces = F;
      else
        [V, F, VN, FN] = longbone_AnatomicalPosition  (this.Vertices, ...
                         this.Faces, this.Normals, this.FaceNormals);
        LBM.Vertices = V;
        LBM.Faces = F;
        LBM.Normals = VN;
        LBM.FaceNormals = FN;
      endif
    endfunction

    function centroid = barycenter (this)
      centroid = meshBarycenter (this.Vertices, this.Faces);
    endfunction

    function D = maxDistance (this)
      D = longbone_maxDistance (this.Vertices);
    endfunction

    function LBM = scale (this, maxD)
      LBM = this;
      D = maxDistance (this);
      ## Calculate scaling ratio
      ratio = maxD / D;
      ## Scale mesh
      LBM.Vertices = this.Vertices * ratio;
    endfunction

    function LBM = renameMesh (this, filename)
      ## Check input arguments
      if (! ischar (filename))
        error ("renameMesh: FILENAME must be a character vector.");
      endif
      this.obj_filename = filename;
      if (! isempty (this.TextureCoords))
        basename = filename([1:end-4]);
        this.mtl_filename = [basename, ".mtl"];
        this.MaterialLib.newmtl = basename
        IMG = this.MaterialLib.map_Kd;
        if (strcmpi (IMG([end-4:end]), ".png"))
          this.MaterialLib.map_Kd = [basename, ".png"];
        elseif (strcmpi (IMG([end-4:end]), ".tiff"))
          this.MaterialLib.map_Kd = [basename, ".tiff"];
        else
          this.MaterialLib.map_Kd = [basename, ".jpg"];
        endif
      endif
      LBM = this;
    endfunction

    function writeObj (this)
      ## Add default options
      if (isempty (this.Normals))
        normals = false;
      else
        normals = true;
      endif
      if (isempty (this.TextureCoords))
        texture = false;
      else
        texture = true;
      endif
      filename = this.obj_filename;
      ## Write to file
      if (! (normals || texture))
        writeObj (this.Vertices, this.Faces, filename);
      elseif (normals && ! texture)
        writeObj (this.Vertices, this.Faces, ...
                  this.Normals, this.FaceNormals, filename);
      elseif (! normals && texture)
        writeObj (this.Vertices, this.Faces, ...
                  this.TextureCoords, this.TextureFaces, filename);
      else
        writeObj (this.Vertices, this.Faces, this.TextureCoords, ...
                  this.TextureFaces, this.Normals, this.FaceNormals, filename);
      endif
      ## Save material library file and texture image
      if (texture)
        mtlwrite (this.mtl_filename, this.MaterialLib);
        imwrite (this.TextureImage, this.MaterialLib.map_Kd);
      endif
    endfunction

  endmethods

  methods (Static, Hidden)

    function LBM = readObj (obj_name, varargin)
      ## Check input arguments
      if (! ischar (obj_name))
        error ("LongBone.readObj: FILENAME must be a character vector.");
      endif
      getinfo = false;
      texture = false;
      while (numel (varargin) > 0)
        if (! ischar (varargin{1}))
          error ("LongBone.readObj: optional args must be character vectors.");
        endif
        switch (tolower (varargin{1}))
          case "info"
            getinfo = true;
          case "texture"
            texture = true;
          otherwise
            error ("LongBone.readObj: unrecognized optional argument.");
        endswitch
        varargin(1) = [];
      endwhile
      ## Read OBJ file
      if (getinfo)
        [V, F, VT, FT, VN, FN, mtl_name] = readObj (obj_name, "info");
      else
        [V, F, VT, FT, VN, FN, mtl_name] = readObj (obj_name);
      endif
      ## Check for color information in vertices
      if (size (V, 2) > 3)
        C = V(:,[4:end]);
        V = V(:,[1:3]);
      else
        C = [];
      endif
      ## Check if texture coordinates and material file exist in OBJ
      if (! (isempty (VT) || isempty (FT) || isempty (mtl_name)) && texture)
        ## Read original mtl file
        MTL = mtlread (mtl_name);
        ## Find which material's name corresponds to OBJ name and keep it or
        ## alternatively search for texture altas in non-empty 'map_Kd' field
        for id = 1:length (MTL)
          if (strcmp (getfield (MTL(id), "newmtl"), obj_name([1:end-4])))
            newmtl_id = id;
          endif
          if (isfield (MTL(id), "map_Kd") &&
              ! isempty (getfield (MTL(id), "map_Kd")))
            map_Kd_id = id;
          endif
        endfor
        if (exist ("newmtl_id"))
          MTL = MTL(newmtl_id);
          ## rename material to the OBJ's updated name (just in case!!)
          MTL.newmtl = obj_name([1:end-4]);
        else
          MTL = MTL(map_Kd_id);
          ## rename material to the OBJ's updated name (just in case!!)
          MTL.newmtl = obj_name([1:end-4]);
        endif
        if (exist ("map_Kd_id"))
          IMG = imread (MTL.map_Kd);
        else
          IMG = [];
        endif
      else
        VT = [];
        FT = [];
        MTL = [];
        IMG = [];
        mtl_name = "";
      endif
      ## Create the LongBone object
      LBM = LongBone (V, F, C, VT, FT, VN, FN, MTL, IMG, obj_name, mtl_name);
    endfunction

  endmethods

endclassdef

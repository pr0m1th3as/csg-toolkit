# csg-toolkit

## A GNU Octave package for analyzing long bone diaphyseal cross sectional geometry

The present package is based on the long-bone-diaphyseal-CSG-Toolkit and has
been created to simplify its installation and usage from within GNU Octave. It
is based on novel and robust algorithms for calculating the cross-sectional
geometric properties of the diaphyses of humerus, ulna, femur, and tibia bones
represented as a triangular mesh in a Wavefront OBJ 3D model file format.

You can find its documentation at [https://pr0m1th3as.github.io/csg-toolkit/](https://pr0m1th3as.github.io/csg-toolkit/).

## Citation

If you find this package usefull for your research, please include the following citations:

Bertsatos A, Chovalopoulou M-E. 2019. A novel method for analyzing long bone diaphyseal 
cross-sectional geometry. A GNU Octave CSG Toolkit. Forensic Science International 297: 65–71. 
https://doi.org/10.1016/j.forsciint.2019.01.041

Bertsatos A, Garoufi N, Koliaraki M, Chovalopoulou M-E. 2023. Paving new ways in forensic contexts with virtual osteology applications: csg-toolkit – a 3D osteology package for cross-sectional geometry analysis. Annals of 3D Printed Medicine, 9, 100094. https://doi.org/10.1016/j.stlm.2022.100094


## Installation

The CSG Toolkit is compatible with later versions of `Octave >= v7.2.0` and depends
on the `io >= 2.6.4` package.

To download and install the latest version issue the following command in Octave:

 `pkg install -forge csg-toolkit`

The package can be loaded on demand in Octave with the following commmand:

 `pkg load csg-toolkit`
 
If you wish to install the latest development version, type:

 `pkg install "https://github.com/pr0m1th3as/csg-toolkit/archive/refs/heads/main.zip"`

Happy long bone analysis!

## Usage

The CSG Toolkit can be called in batch processing mode using the `longbone_Analysis`
script which can handle all available 3D models found in OBJ format in the working
directory or any other folder path. Although this functionality is still available,
there is no longer need to work with a specific bone being present in a folder for
batch proccessing. The user may select a specific bone, i.e. Humerus, or may choose
to process different bones in a single batch process. Not selected bones are omitted
from batch processing.

After batch processing, the function `inspect_CSG` can be used for visual inspection
of the calculated cross-sectional contours according to the files present in the
working directory. Furthermore, it will assemble all calculated values from the
the samples (bones) that have been verified by the user and save them in a CSV file.

The need for MeshLab .pp side car files has been also relaxed. If they are present
along with their OBJ counterparts, they are utilized, if not, the initial alignment
points are automatically registered with the `longbone_Registration` function.
Moreover, the new function `longbone_CustomGeometry` allows for arbitrary user
defined cross sections along the longitudinal axis of the bones.

The CSG Toolkit also provides functionality for graphical representation of the
cross-sectional contours and their respective CSG properties, which can be ploted
with the `visualize_CrossSections` function, which utilizes the output results of
either `longbone_Analysis` stored in relevant CSV files, or the return structures
of the `longbone_Geometry` function and the generated CSV files from the
`longbone_CustomGeometry` function.

## Functions available in csg-toolkit

These are the functions available when csg-toolkit is loaded.
At the Octave command prompt, type `help function-name` for more information
about the function and it's usage.

* `inspect_CSG`
* `longbone_Analysis`
* `longbone_CustomGeometry`
* `longbone_FragmentGeometry`
* `longbone_Geometry`
* `longbone_maxDistance`
* `longbone_Registration`
* `longbone_Scaling`
* `meshArea`
* `meshBarycenter`
* `meshSection`
* `read_MeshlabPoints`
* `readObj`
* `renameObj`
* `simple_polygon3D`
* `visualize_CrossSections`
* `write_MeshlabPoints`
* `writeObj`

## Testing Datasets

The present toolkit has been extensively tested and these results are presented
in the relevant validation study available at https://doi.org/10.5281/zenodo.1466135
The initial testing dataset, which comprises three 3D mesh models of real humans
bones (a humerus, a femur and a tibia, which are part of the Athens modern reference
skeletal collection), is also freely available at https://doi.org/10.5281/zenodo.1466962
and may be used with the CSG Toolkit for demonstrating its operation. 

An additional dataset comprising eight bilateral bones (humerus, ulna, femur, and
tibia) is also available at https://doi.org/10.5281/zenodo.6614647 for testing and
demostration purposes.  Simplified versions of this dataset, reduced to 7,500 faces
without any color, texture and normals information, is also available within the
package and available once loaded in Octave.  This is only available since version
1.3.0.


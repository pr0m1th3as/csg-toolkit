# csg-toolkit

## A GNU Octave package for analyzing long bone diaphyseal cross sectional geometry

The present package is based on the [long-bone-diaphyseal-CSG-Toolkit](https://github.com/pr0m1th3as/long-bone-diaphyseal-CSG-Toolkit)
and has been created to simplify its installation and usage from within GNU Octave. Itis based on novel and robust algorithms for
calculating the cross-sectional geometric properties of the diaphyses of humerus, ulna, femur, and tibia bones represented as a
triangular mesh in a Wavefront OBJ 3D model file format.

This package provides numerous functions for analyzing CSG properties and working with 3D models. You can find its documentation at
[https://pr0m1th3as.github.io/csg-toolkit/](https://pr0m1th3as.github.io/csg-toolkit/).

## Installation

The CSG Toolkit is compatible with later versions of [`Octave >= v9.1.0`](https://octave.org/) and depends on the 
[`io >= 2.6.4`](https://gnu-octave.github.io/packages/io/), [`statistics >=1.7.4`](https://gnu-octave.github.io/packages/statistics/),
and [`datatypes >=1.0.1`](https://github.com/pr0m1th3as/datatypes/) packages.

To download and install the latest version, type:

 `pkg install -forge csg-toolkit`

The package can be loaded on demand with the following commmand:

 `pkg load csg-toolkit`

Happy long bone analysis!

## Usage

The CSG Toolkit can be called in batch processing mode using the `longbone_Analysis`
function which can handle all available 3D models found in OBJ format in any chosen
directory.  Batch processing is available for all three modes of Geometry:

1. Default:  utilizing the `longbone_Geometry` function for analyzing cross sections
at the standard 20%, 35%, 50%, 65%, and 80% intervals.
2. Custom:   utilizing the `longbone_CustomGeomtery` function for analyzing arbitrary
cross sections at any user specified interval between 10% and 90%.
3. Fragment: utilizing the `longbone_FragmentGeometry` function for analysing cross
sections from arbitrary long bone fragments.

Although this functionality is still available, there is no longer need to work with
a specific bone being present in a folder for batch proccessing. The user may select
a specific bone, i.e. Humerus, or may choose to process different bones in a single
batch process. Bone models that do not match the user selection are omitted from
batch processing.

After batch processing, the function `inspect_CSG` can be used for visual inspection
of the calculated cross-sectional contours according to the CSV files present in the
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

You may also use the CSG Toolkit for biological sex estimation.
Type `help longbone_Sex` to find out how to utilize CSG properties for biological
profiling from the long bones. More functionality comming soon!

## Testing Datasets

The present toolkit has been extensively tested and these results are presented
in the relevant validation study available at https://doi.org/10.5281/zenodo.1466135
The initial testing dataset, which comprises three 3D mesh models of real humans
bones (a humerus, a femur and a tibia, which are part of the Athens modern reference
skeletal collection), is also freely available at https://doi.org/10.5281/zenodo.1466962
and may be used with the CSG Toolkit for demonstrating its operation. 

An additional dataset comprising eight bilateral bones (humerus, ulna, femur, and
tibia) is also available at https://doi.org/10.5281/zenodo.6614647 for testing and
demostration purposes.

## Citation

If you find this package usefull for your research, please include the following citations:

Bertsatos A, Garoufi N, Koliaraki M, Chovalopoulou M-E. 2023. Paving new ways in forensic contexts with virtual osteology applications: csg-toolkit – a 3D osteology package for cross-sectional geometry analysis. Annals of 3D Printed Medicine, 9, 100094. https://doi.org/10.1016/j.stlm.2022.100094

Bertsatos A, Chovalopoulou M-E. 2019. A novel method for analyzing long bone diaphyseal 
cross-sectional geometry. A GNU Octave CSG Toolkit. Forensic Science International 297: 65–71. 
https://doi.org/10.1016/j.forsciint.2019.01.041

## Acknowledgements

Release [1.4.0](https://github.com/pr0m1th3as/csg-toolkit/releases/tag/v1.4.0) of the **csg-toolkit** package has incorporated part of the research output of my [RECONSTRUCT](https://www.physicalanthropology.gr/reconstruct.php) project, which has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie Grant agreement No 101104702.

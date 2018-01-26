# Redesign3D
- Makes local changes in 3D models to improve manufacturability
- Uses CGAL and armadillo libraries
- For 2D version, please see <https://github.com/prakhar0402/redesign2d>

### Abstract
The recent advancement of additive manufacturing technology has allowed for fabrication of highly complex parts that were difficult to produce using traditional manufacturing techniques. However, designers and novice users of additive technologies often lack the awareness of manufacturing considerations leading to lower quality parts and failures. Therefore, it is crucial to develop an automated system that can analyze and correct designs before they are send for fabrication so that wasteful iterations of *build-test-redesign* can be minimized. In this project, a variation of **shape diameter function** and **morphological operations** are used to identify critical regions in 2D slices of 3D models that can potentially cause errors in printing. Printability of slices has been quantified using an area-based printability index. In addition, a **physics-based mesh deformation scheme** is adopted to make localized corrections to slices for improving printability. The approach has also been extended to three-dimensions for correcting critically thin regions of 3D models.

Please look into `requirements.txt` for the list of required libraries and installation instructions.

Please leave your comments and feedback at <prakharj@buffalo.edu>.

If you use our code or part of our code, please cite us.  A reference to our publication to cite will be included here soon.

MAD Lab, University at Buffalo
Copyright (C) 2018  Prakhar Jaiswal <prakharj@buffalo.edu>
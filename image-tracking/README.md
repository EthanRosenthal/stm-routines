image-tracking
===========

This routine warps a series of moving images to match them to a fixed image. The routine's primary purpose was for overlaying scanning tunneling microscopy images of the same region of the sample acquired at different temperatures. The multiple hours required for each image acquisition produces significant drift and skewing of images. Combined with a significant temperature dependence to the image scaling (due to temperature dependent piezoelectric coefficients), every image has to be warped to a reference image.

"Map" in the routine refers to a dI/dV map. This is a 3-D matrix where the first two dimension correspond to the spatial dimensions, and the last dimension can be thought of as an array in which a particular imaging parameter (the bias, V) is varied. Consequently, each fixed value of the last dimension corresponds to a spatial "slice" for that particular value. 

"Topo" refers to a topographic image (2D matrix).

Directions for using the routine are included in the comments at the beginning.

# timepixcam

This C code was created for processing data from Timepix3Cam for Fluorescence Lifetime Imaging (FLIM).
It can also calculate centroid positions photon events.

Reads ASCII files with columns of photon event data, and organises the data into an xyt cube.
The output is an .ics file that can be loaded in Tri2 for fluorescence lifetime analysis.
Also outputs an ASCII file containing a histogram of the photon arrival times and an intensity image.
Besides sum (all pixels), also outputs centroided data, where only the centroid pixel of each cluster (photon event) is counted.

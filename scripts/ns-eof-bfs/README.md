For the postprocessing of the BFS results following scripts were created:

1) extractdata.m (Matlab)
It extracts velocities in x-direction from the VTK-files along a line (parallel to y, x position can be specified) for every timestep. The values are written to a csv-file.

2) (Phython)
Does the same as 1), except that it can be run on the cluster. It is not necessary to copy the files from the cluster. (see evaluation/evaluate)

2) postkbatch.m (Matlab)
Creates plots (mean, variance, skewness, kurtosis)

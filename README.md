# Voronoi-Delaunay_tesselation

This method provides a scale independent way to measure the underlying density field of a sample of data points. This can be used to calculate the 2D surface density of a given sample of data.
The Voronoi tessellation divides a 2D distribution of points into convex cells, with each cell containing only one point and a set of vertices which are closer to that point than to any other in the plane. It has the property that the local density of each cell is the inverse of the cell area. Thus, to estimate the overdensity of each cell, I first calculate the average density of the cells in the entire plane and I then compute the contrast of each cell.

In my case, I apply the method to a sample of galaxies with 3D positions i.e., x=right ascension, y=declination, z=redshift. Thus, I first slice the sample into z slices to obtain 2D planes of data.

See NEWS file for updates in latest version.

The reason for the Non-FOSS licence is because of the dependency on
RTriangle:
https://cran.r-project.org/web/packages/RTriangle/index.html a wrapper
for the Triangle library http://www.cs.cmu.edu/~quake/triangle.html
which has a non-FOSS licence. As maintainer of RTriangle, I tried, but
failed to persuade the author of Triangle to make it available under a
FOSS licence.

I've not found an R library that does constrained 2D mesh generation
as efficiently with the required flexibility to omit holes and
constrain angles for high quality meshes.

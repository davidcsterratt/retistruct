
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title>Retistruct - computational reconstruction and transformation of flattened retinae</title>
  </head>

<body>

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<img src="folding-small.png" style="float: left"><h1>Retistruct</h1>
<p class="slogan">-- computational reconstruction and transformation of flattened retinae</p>

<!-- end of project description -->

<p style="clear: both">Retistruct is
 an <a href="http://www.r-project.org">R</a> package to morph a flat
 surface with cuts (a dissected flat-mount retina) onto a
 curvilinear surface (the a standard retinal shape).  It can estimate
 the position of a point on the intact adult retina to within
 8&deg; of arc
 (3.6% of nasotemporal axis). The coordinates in reconstructed retinae
 can be transformed to visuotopic coordinates.
 </p>

<h2>How Retistruct works</h2>

<img src="retistruct.small.png" style="float: right"><p>Reconstruction
is achieved by: stitching the marked-up cuts of the flat-mount
outline; dividing the stitched outline into a mesh whose vertices then
are mapped onto a curtailed sphere; and finally moving the vertices so
as to minimise a physically-inspired deformation energy function. <ul>
      <li>Full details are in the paper: Sterratt, D. C., Lyngholm,
        D., Willshaw, D. J. and Thompson, I. D. (2013).  Standard
        Anatomical and Visual Space for the Mouse Retina:
        Computational Reconstruction and Transformation of Flattened
        Retinae with the Retistruct
        Package. <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002921"><em>PLoS
        Computational Biology</em> 9:e1002921</a>.</li>
      <li><a href="2012-09-neuroinf.pdf">This poster</a>
and <a href="http://www.youtube.com/watch?v=LpuqXo8NEOo">this YouTube
video</a>, presented at the 2012 Neuroinformatics Meeting in Munich,
has more details and examples of reconstructions and projections into
visual space. </li>
      </ul></p>

<h2>Installation and documentation</h2>

<p>Retistruct has been tested on GNU/Linux (Ubuntu 12.04), Mac OS X
  10.8 and Microsoft Windows Vista. <strong>(Note that as of 30th
  October 2013, it is not possible to install the graphical user
  interface on MacOS X 10.9 (Mavericks). This is due to the RGtk2
  package upon which the Retistruct GUI depends.)</strong></p>

<h3>Stable version</h3>

<p>To install
  the <a href="http://cran.r-project.org/web/packages/retistruct/index.html">stable
  version of Retistruct hosted on CRAN</a>, follow the instructions in
  the <a href="retistruct-user-guide.pdf">User Guide</a>. The
  installation contains a number of demonstration retinae, and
  instructions for how to handle retinal flat-mount images in
  Retistruct.</p>

<h3>Development version</h3>

<p>The development version of Retistruct contains the most recent bug
  fixes and improvements, but is not stable. Builds of the package are
  usually <a href="https://r-forge.r-project.org/R/?group_id=1436">available
  on R-forge</a> and can be installed using the R install command
  given there. If the R-forge packages are not available, you will
  have to build from the source code using <tt>R CMD build</tt>.</p>

<h3>Source code</h3>

<p>The source code can be checked out from
the <a href="https://r-forge.r-project.org/scm/?group_id=1436">R-forge
subversion repository</a>
or <a href="https://github.com/davidcsterratt/retistruct">Github</a>.</p>

<h3>Reference publication code</h3>

<p>For reference purposes, this <a href="retistruct_0.5.7.zip"> zip
file</a> contains the version of Retistruct that generated the
reconstructions in Sterratt &amp;
al. (2013; <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002921"><em>PLoS
Computational Biology</em> 9</a>). The file also contains some Matlab
code to read data directories contained by Retistruct.</p>

<h2>Sample data</h2>

<p>As well as the built-in demo data, there are some sample images to
  practise on:<ol>
  <li>Beginner: <a href="data/image.png">SMI-32 stained retina</a>. As described in
  the <a href="retistruct-user-guide.pdf">User Guide</a>, the outline
  can be marked up in <a href="http://rsb.info.nih.gov/ij/">ImageJ</a>
  and imported into Retistruct.</li>
  <li>More advanced: TIFF files
  (<a href="data/left-5x-small.tif">left</a> and
  <a href="data/right-5x-small.tif">right</a>), each containing a stack of three images
  corresponding to Figure 6 of the manuscript: retinae labelled with
  Fluoro-Emerald, Fluoro-Ruby and a brightfield image. As described in
  the <a href="retistruct-user-guide.pdf">User Guide</a>,
  use <a href="http://rsb.info.nih.gov/ij/">ImageJ</a> to mark up the
  outline on the brightfield image, and then
  use <a href="http://rsb.info.nih.gov/ij/">ImageJ's</a> particle
  analyser to find the locations of the stained cells.</li>
  </ol></p>

<h2>Problems?</h2>

If you encounter issues using Retistruct please either:
<ul>
  <li><a href="https://github.com/davidcsterratt/retistruct/issues">Submit
  an issue on Github</a></li>
  <li><a href="mailto:david.c.sterratt@ed.ac.uk">Email David Sterratt</a></li>
</ul>

<h2>Work using Retistruct</h2>

<ul>
<li>Flinn, J. M., Kakalec P., Tappero, R., Jones, B. and Lengyel,
I. (2014) "Correlations in distribution and concentration of calcium,
copper and iron with zinc in isolated extracellular deposits
associated with age-related macular
degeneration". <em>Metallomics</em>. <a href="http://dx.doi.org/10.1039/c4mt00058g">DOI:
10.1039/c4mt00058g</a>
 </li>
<li>Bleckert, A.,  Schwartz, G. W., Turner, M.H., Rieke, F., Wong,
  R. O. L. (2014) "Visual space is represented by nonmatching topographies of
  distinct mouse retinal ganglion cell types". <em>Current
  Biology</em> <strong>24</strong>:
  310-315.  <a href="http://dx.doi.org/10.1016/j.cub.2013.12.020">DOI: 10.1016/j.cub.2013.12.020</a>
</li>
</ul>

<h2>Authors and funding</h2>

<p>Retistruct was written
  by <a href="http://homepages.inf.ed.ac.uk/sterratt/">David
  Sterratt</a> at the <a href="http://www.ed.ac.uk/">University of
  Edinburgh</a>, and tested by Daniel Lyngholm and Ian Thompson at the
  <a href="http://www.kcl.ac.uk/depsta/biomedical/mrc/">MRC Centre for
  Developmental Neurobiology, KCL</a>.
 </p>

<p>This work was supported by a Programme Grant from
the <a href="http://www.wellcome.ac.uk">Wellcome Trust</a> (G083305).
</p>



</body>
</html>

<!--  LocalWords:  Retistruct YouTube Sterratt Lyngholm MRC KCL PLoS
 -->
<!--  LocalWords:  visuotopic Willshaw Ubuntu th MacOS RGtk CRAN CMD
 -->
<!--  LocalWords:  Github Matlab SMI ImageJ Fluoro brightfield DOI
 -->
<!--  LocalWords:  ImageJ's Kakalec Tappero Lengyel macular Bleckert
 -->
<!--  LocalWords:  Metallomics Rieke nonmatching
 -->

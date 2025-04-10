% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TriangulatedOutline.R
\name{TriangulatedOutline}
\alias{TriangulatedOutline}
\title{Class containing functions and data relating to Triangulation}
\description{
A TriangulatedOutline contains a function to create a
  triangulated mesh over an outline, and fields to hold the mesh
  information. Note that areas and lengths are all scaled using
  the value of the \code{scale} field.
}
\examples{
P <- rbind(c(1,1),   c(2,1),  c(2,-1),
           c(1,-1),  c(1,-2), c(-1,-2),
           c(-1,-1), c(-2,-1),c(-2,1),
           c(-1,1),  c(-1,2), c(1,2))
o <- TriangulatedOutline$new(P)
o$addTear(c(3, 4, 5))
o$addTear(c(6, 7, 8))
o$addTear(c(9, 10, 11))
o$addTear(c(12, 1, 2))
flatplot(o)

P <- list(rbind(c(1,1), c(2,1), c(2.5,2), c(3,1), c(4,1), c(1,4)),
              rbind(c(-1,1), c(-1,4), c(-2,3), c(-2,2), c(-3,2), c(-4,1)),
              rbind(c(-4,-1), c(-1,-1), c(-1,-4)),
              rbind(c(1,-1), c(2,-1), c(2.5,-2), c(3,-1), c(4,-1), c(1,-4)))
o <- TriangulatedOutline$new(P)
##' o$addTear(c(2, 3, 4))
o$addTear(c(17, 18, 19))
o$addTear(c(9, 10, 11))
o$addFullCut(c(1, 5, 16, 20))
flatplot(o)
}
\author{
David Sterratt
}
\section{Super classes}{
\code{\link[retistruct:OutlineCommon]{retistruct::OutlineCommon}} -> \code{\link[retistruct:Outline]{retistruct::Outline}} -> \code{\link[retistruct:PathOutline]{retistruct::PathOutline}} -> \code{\link[retistruct:AnnotatedOutline]{retistruct::AnnotatedOutline}} -> \code{TriangulatedOutline}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{Tr}}{3 column matrix in which each row contains IDs of
points of each triangle}

\item{\code{A}}{Area of each triangle in the mesh - has same number of
elements as there are rows of \code{T}}

\item{\code{A.tot}}{Total area of the mesh}

\item{\code{Cu}}{2 column matrix in which each row contains IDs of}

\item{\code{L}}{Length of each edge in the mesh - has same number of
elements as there are rows of \code{Cu}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TriangulatedOutline-triangulate}{\code{TriangulatedOutline$triangulate()}}
\item \href{#method-TriangulatedOutline-mapTriangulatedFragment}{\code{TriangulatedOutline$mapTriangulatedFragment()}}
\item \href{#method-TriangulatedOutline-clone}{\code{TriangulatedOutline$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="OutlineCommon" data-id="clearFeatureSets"><a href='../../retistruct/html/OutlineCommon.html#method-OutlineCommon-clearFeatureSets'><code>retistruct::OutlineCommon$clearFeatureSets()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="OutlineCommon" data-id="getFeatureSet"><a href='../../retistruct/html/OutlineCommon.html#method-OutlineCommon-getFeatureSet'><code>retistruct::OutlineCommon$getFeatureSet()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="OutlineCommon" data-id="getFeatureSetTypes"><a href='../../retistruct/html/OutlineCommon.html#method-OutlineCommon-getFeatureSetTypes'><code>retistruct::OutlineCommon$getFeatureSetTypes()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="OutlineCommon" data-id="getFeatureSets"><a href='../../retistruct/html/OutlineCommon.html#method-OutlineCommon-getFeatureSets'><code>retistruct::OutlineCommon$getFeatureSets()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="OutlineCommon" data-id="getIDs"><a href='../../retistruct/html/OutlineCommon.html#method-OutlineCommon-getIDs'><code>retistruct::OutlineCommon$getIDs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="addFeatureSet"><a href='../../retistruct/html/Outline.html#method-Outline-addFeatureSet'><code>retistruct::Outline$addFeatureSet()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getDepth"><a href='../../retistruct/html/Outline.html#method-Outline-getDepth'><code>retistruct::Outline$getDepth()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getFragment"><a href='../../retistruct/html/Outline.html#method-Outline-getFragment'><code>retistruct::Outline$getFragment()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getFragmentIDs"><a href='../../retistruct/html/Outline.html#method-Outline-getFragmentIDs'><code>retistruct::Outline$getFragmentIDs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getFragmentIDsFromPointIDs"><a href='../../retistruct/html/Outline.html#method-Outline-getFragmentIDsFromPointIDs'><code>retistruct::Outline$getFragmentIDsFromPointIDs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getFragmentPointIDs"><a href='../../retistruct/html/Outline.html#method-Outline-getFragmentPointIDs'><code>retistruct::Outline$getFragmentPointIDs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getFragmentPoints"><a href='../../retistruct/html/Outline.html#method-Outline-getFragmentPoints'><code>retistruct::Outline$getFragmentPoints()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getImage"><a href='../../retistruct/html/Outline.html#method-Outline-getImage'><code>retistruct::Outline$getImage()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getOutlineLengths"><a href='../../retistruct/html/Outline.html#method-Outline-getOutlineLengths'><code>retistruct::Outline$getOutlineLengths()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getOutlineSet"><a href='../../retistruct/html/Outline.html#method-Outline-getOutlineSet'><code>retistruct::Outline$getOutlineSet()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getPoints"><a href='../../retistruct/html/Outline.html#method-Outline-getPoints'><code>retistruct::Outline$getPoints()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getPointsScaled"><a href='../../retistruct/html/Outline.html#method-Outline-getPointsScaled'><code>retistruct::Outline$getPointsScaled()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="getPointsXY"><a href='../../retistruct/html/Outline.html#method-Outline-getPointsXY'><code>retistruct::Outline$getPointsXY()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="mapFragment"><a href='../../retistruct/html/Outline.html#method-Outline-mapFragment'><code>retistruct::Outline$mapFragment()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="mapPids"><a href='../../retistruct/html/Outline.html#method-Outline-mapPids'><code>retistruct::Outline$mapPids()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="Outline" data-id="replaceImage"><a href='../../retistruct/html/Outline.html#method-Outline-replaceImage'><code>retistruct::Outline$replaceImage()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="PathOutline" data-id="insertPoint"><a href='../../retistruct/html/PathOutline.html#method-PathOutline-insertPoint'><code>retistruct::PathOutline$insertPoint()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="PathOutline" data-id="nextPoint"><a href='../../retistruct/html/PathOutline.html#method-PathOutline-nextPoint'><code>retistruct::PathOutline$nextPoint()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="PathOutline" data-id="stitchSubpaths"><a href='../../retistruct/html/PathOutline.html#method-PathOutline-stitchSubpaths'><code>retistruct::PathOutline$stitchSubpaths()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="addFullCut"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-addFullCut'><code>retistruct::AnnotatedOutline$addFullCut()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="addPoints"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-addPoints'><code>retistruct::AnnotatedOutline$addPoints()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="addTear"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-addTear'><code>retistruct::AnnotatedOutline$addTear()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="checkTears"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-checkTears'><code>retistruct::AnnotatedOutline$checkTears()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="computeFullCutRelationships"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-computeFullCutRelationships'><code>retistruct::AnnotatedOutline$computeFullCutRelationships()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="computeTearRelationships"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-computeTearRelationships'><code>retistruct::AnnotatedOutline$computeTearRelationships()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="ensureFixedPointInRim"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-ensureFixedPointInRim'><code>retistruct::AnnotatedOutline$ensureFixedPointInRim()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getBoundarySets"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getBoundarySets'><code>retistruct::AnnotatedOutline$getBoundarySets()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getFixedPoint"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getFixedPoint'><code>retistruct::AnnotatedOutline$getFixedPoint()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getFullCut"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getFullCut'><code>retistruct::AnnotatedOutline$getFullCut()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getFullCuts"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getFullCuts'><code>retistruct::AnnotatedOutline$getFullCuts()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getRimLengths"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getRimLengths'><code>retistruct::AnnotatedOutline$getRimLengths()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getRimSet"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getRimSet'><code>retistruct::AnnotatedOutline$getRimSet()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getTear"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getTear'><code>retistruct::AnnotatedOutline$getTear()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="getTears"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-getTears'><code>retistruct::AnnotatedOutline$getTears()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="initialize"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-initialize'><code>retistruct::AnnotatedOutline$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="labelFullCutPoints"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-labelFullCutPoints'><code>retistruct::AnnotatedOutline$labelFullCutPoints()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="labelTearPoints"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-labelTearPoints'><code>retistruct::AnnotatedOutline$labelTearPoints()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="removeFullCut"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-removeFullCut'><code>retistruct::AnnotatedOutline$removeFullCut()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="removeTear"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-removeTear'><code>retistruct::AnnotatedOutline$removeTear()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="setFixedPoint"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-setFixedPoint'><code>retistruct::AnnotatedOutline$setFixedPoint()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="whichFullCut"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-whichFullCut'><code>retistruct::AnnotatedOutline$whichFullCut()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="AnnotatedOutline" data-id="whichTear"><a href='../../retistruct/html/AnnotatedOutline.html#method-AnnotatedOutline-whichTear'><code>retistruct::AnnotatedOutline$whichTear()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TriangulatedOutline-triangulate"></a>}}
\if{latex}{\out{\hypertarget{method-TriangulatedOutline-triangulate}{}}}
\subsection{Method \code{triangulate()}}{
Triangulate (mesh) outline
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TriangulatedOutline$triangulate(n = 200, suppress.external.steiner = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{Desired number of points in mesh}

\item{\code{suppress.external.steiner}}{Boolean variable describing
whether to insert external Steiner points - see
\link{TriangulatedFragment}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TriangulatedOutline-mapTriangulatedFragment"></a>}}
\if{latex}{\out{\hypertarget{method-TriangulatedOutline-mapTriangulatedFragment}{}}}
\subsection{Method \code{mapTriangulatedFragment()}}{
Map the point IDs of a \link{TriangulatedFragment} on the
  point IDs of this Outline
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TriangulatedOutline$mapTriangulatedFragment(fragment, pids)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fragment}}{\link{TriangulatedFragment} to map}

\item{\code{pids}}{Point IDs in TriangulatedOutline of points in \link{TriangulatedFragment}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TriangulatedOutline-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TriangulatedOutline-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TriangulatedOutline$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

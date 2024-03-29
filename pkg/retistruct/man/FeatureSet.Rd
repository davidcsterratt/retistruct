% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSet.R
\name{FeatureSet}
\alias{FeatureSet}
\title{Superclass containing functions and data relating to sets of
features in flat \code{\link{Outline}}s}
\description{
A FeatureSet contains information about features
  located on \code{\link{Outline}}s. Each FeatureSet contains a
  list of matrices, each of which has columns labelled \code{X}
  and \code{Y} describing the cartesian coordinates of points on
  the Outline, in the unscaled coordinate frame. Derived classes,
  e.g. a \code{\link{CountSet}}, may have extra columns. Each matrix
  in the list has an associated label and colour, which is used by
  plotting functions.
}
\author{
David Sterratt
}
\section{Super class}{
\code{\link[retistruct:FeatureSetCommon]{retistruct::FeatureSetCommon}} -> \code{FeatureSet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-FeatureSet-new}{\code{FeatureSet$new()}}
\item \href{#method-FeatureSet-clone}{\code{FeatureSet$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="getCol"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-getCol'><code>retistruct::FeatureSetCommon$getCol()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="getFeature"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-getFeature'><code>retistruct::FeatureSetCommon$getFeature()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="getFeatures"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-getFeatures'><code>retistruct::FeatureSetCommon$getFeatures()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="getIDs"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-getIDs'><code>retistruct::FeatureSetCommon$getIDs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="getIndex"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-getIndex'><code>retistruct::FeatureSetCommon$getIndex()</code></a></span></li>
<li><span class="pkg-link" data-pkg="retistruct" data-topic="FeatureSetCommon" data-id="setID"><a href='../../retistruct/html/FeatureSetCommon.html#method-FeatureSetCommon-setID'><code>retistruct::FeatureSetCommon$setID()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FeatureSet-new"></a>}}
\if{latex}{\out{\hypertarget{method-FeatureSet-new}{}}}
\subsection{Method \code{new()}}{
Constructor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FeatureSet$new(data = NULL, cols = NULL, type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{List of matrices describing data. Each matrix
should have columns named \code{X} and \code{Y}}

\item{\code{cols}}{Named vector of colours for each data set. The name is
used as the ID (label) for the data set. The colours should be names
present in the output of the \code{\link{colors}} function}

\item{\code{type}}{String}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FeatureSet-clone"></a>}}
\if{latex}{\out{\hypertarget{method-FeatureSet-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FeatureSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

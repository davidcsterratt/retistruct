% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OutlineCommon.R
\name{OutlineCommon}
\alias{OutlineCommon}
\title{Class containing functionality common to flat and reconstructed outlines}
\description{
An OutlineCommon has functionality for retrieving
  sets of features (e.g. points or landmarks associated with an
  outline)
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{version}}{Version of reconstruction file data format}

\item{\code{featureSets}}{List of feature sets associated with the outline, which may be of various types, e.g. a \link{PointSet} or \link{LandmarkSet}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-OutlineCommon-getFeatureSets}{\code{OutlineCommon$getFeatureSets()}}
\item \href{#method-OutlineCommon-getFeatureSet}{\code{OutlineCommon$getFeatureSet()}}
\item \href{#method-OutlineCommon-clearFeatureSets}{\code{OutlineCommon$clearFeatureSets()}}
\item \href{#method-OutlineCommon-getIDs}{\code{OutlineCommon$getIDs()}}
\item \href{#method-OutlineCommon-getFeatureSetTypes}{\code{OutlineCommon$getFeatureSetTypes()}}
\item \href{#method-OutlineCommon-clone}{\code{OutlineCommon$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-getFeatureSets"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-getFeatureSets}{}}}
\subsection{Method \code{getFeatureSets()}}{
Get all the feature sets
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$getFeatureSets()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
List of \link{FeatureSet}s associated with the outline
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-getFeatureSet"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-getFeatureSet}{}}}
\subsection{Method \code{getFeatureSet()}}{
Get all feature sets of a particular type, e.g. \link{PointSet} or \link{LandmarkSet}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$getFeatureSet(type)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{The type of the feature set as a string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
All \link{FeatureSet}s of that type
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-clearFeatureSets"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-clearFeatureSets}{}}}
\subsection{Method \code{clearFeatureSets()}}{
Clear all feature sets from the outline
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$clearFeatureSets()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-getIDs"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-getIDs}{}}}
\subsection{Method \code{getIDs()}}{
Get all the distinct IDs contained in the \link{FeatureSet}s
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$getIDs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Vector of IDs
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-getFeatureSetTypes"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-getFeatureSetTypes}{}}}
\subsection{Method \code{getFeatureSetTypes()}}{
Get all the distinct types of \link{FeatureSet}s
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$getFeatureSetTypes()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Vector of types as strings, e.g. \var{PointSet}, \var{LandmarkSet}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-OutlineCommon-clone"></a>}}
\if{latex}{\out{\hypertarget{method-OutlineCommon-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{OutlineCommon$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

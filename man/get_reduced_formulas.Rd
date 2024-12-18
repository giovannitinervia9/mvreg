% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvreg_utilities.R
\name{get_reduced_formulas}
\alias{get_reduced_formulas}
\title{Generate reduced formulas by removing terms iteratively}
\usage{
get_reduced_formulas(response, formula)
}
\arguments{
\item{response}{A character string representing the response variable to be used
on the left-hand side (LHS) of the generated formulas.}

\item{formula}{A formula object representing the full model from which the reduced formulas will be derived.}
}
\value{
A list of formula objects, each representing a reduced version of the input formula.
The formulas are ordered from the most complex (original formula) to the simplest (only intercept).
}
\description{
This function takes an initial formula and generates a list of reduced formulas by
successively removing one term at a time from the right-hand side (RHS) of the formula.
}
\details{
The \code{get_reduced_formulas} function systematically removes the last term
from the RHS of the input formula, creating a series of progressively simpler formulas.
It continues until only the intercept remains on the RHS. This approach allows for
easy generation of nested models, which can be used for comparing different model sizes.
}
\examples{
# Example with interaction terms
formula <- y ~ x1 * x2 * x3
reduced_formulas <- get_reduced_formulas("y", formula)
print(reduced_formulas)
}
\seealso{
\code{\link[stats]{formula}} for more details on formula objects.
}

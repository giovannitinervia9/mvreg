#-------------------------------------------------------------------------------


#' Generate reduced formulas by removing terms iteratively
#'
#' This function takes an initial formula and generates a list of reduced formulas by
#' successively removing one term at a time from the right-hand side (RHS) of the formula.
#'
#' @param response A character string representing the response variable to be used
#'        on the left-hand side (LHS) of the generated formulas.
#' @param formula A formula object representing the full model from which the reduced formulas will be derived.
#'
#' @return A list of formula objects, each representing a reduced version of the input formula.
#'         The formulas are ordered from the most complex (original formula) to the simplest (only intercept).
#'
#' @export
#'
#'
#' @details
#' The \code{get_reduced_formulas} function systematically removes the last term
#' from the RHS of the input formula, creating a series of progressively simpler formulas.
#' It continues until only the intercept remains on the RHS. This approach allows for
#' easy generation of nested models, which can be used for comparing different model sizes.
#'
#' @seealso \code{\link[stats]{formula}} for more details on formula objects.
#'
#'
#' @examples
#' # Example with interaction terms
#' formula <- y ~ x1*x2*x3
#' reduced_formulas <- get_reduced_formulas("y", formula)
#' print(reduced_formulas)
get_reduced_formulas <- function(response, formula) {

  # Extract the terms from the formula
  terms_formula <- attr(terms(formula), "term.labels")

  # Initialize a list to store the resulting formulas
  formulas <- list()

  # Loop through each term and remove the last one successively
  for (i in seq_along(terms_formula)) {
    terms_subset <- terms_formula[1:(length(terms_formula) - i + 1)]
    formula_str <- paste(paste0(response, " ~"), paste(terms_subset, collapse = " + "))
    formulas[[i]] <- as.formula(formula_str)
  }

  # Add the formula with only the intercept
  formulas[[length(formulas) + 1]] <- as.formula(paste0(response, " ~ 1"))

  return(formulas)
}

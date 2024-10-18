#-------------------------------------------------------------------------------


#' Generate formulas by dropping terms one by one
#'
#' This function generates a list of formulas by iteratively removing the last term
#' from the right-hand side of the provided formula, returning a series of reduced formulas.
#'
#' @param response A character string representing the response variable to be used on the left-hand side of the formula.
#' @param formula A formula object from which to generate the reduced formulas.
#'
#' @return A list of formula objects, each representing a reduced version of the input formula by removing terms one by one.
#' @export
#'
#' @examples
#' formula <- y ~ x1 * x2 * x3
#' response <- "y"
#' get_reduced_formulas(response, formula)
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

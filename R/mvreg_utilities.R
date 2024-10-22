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


#-------------------------------------------------------------------------------


#' Check if Two Models are Nested
#'
#' This function determines whether one model is nested within another. A model
#' is considered nested if all of its terms are a subset of the terms in the
#' other model.
#'
#' @param model1 An object of class \code{lm}, \code{glm}, or another model
#'        type with a \code{coef} method. Represents the first model to be
#'        compared.
#' @param model2 An object of class \code{lm}, \code{glm}, or another model
#'        type with a \code{coef} method. Represents the second model to be
#'        compared.
#'
#' @return A logical value: \code{TRUE} if one model is nested within the other
#'         (i.e., all the terms in one model are contained within the other),
#'         and \code{FALSE} otherwise.
#'
#' @examples
#' # Example usage with linear models
#' model1 <- lm(mpg ~ wt, data = mtcars)
#' model2 <- lm(mpg ~ wt + hp, data = mtcars)
#' is_nested(model1, model2) # TRUE, model1 is nested within model2
#'
#' model3 <- lm(mpg ~ qsec, data = mtcars)
#' is_nested(model1, model3) # FALSE, no nesting between model1 and model3
#'
#' @export
is_nested <- function(model1, model2) {
  terms1 <- names(coef(model1))
  terms2 <- names(coef(model2))
  all(terms1 %in% terms2) | all(terms2 %in% terms1)
}


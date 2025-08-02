## Experimental function --- generic but crude approach for experimentation

## Main inputs:
## - dataset to be analyzed
## - baseline model formula as character string
## - character vector of exposure variables

## add var to the model and get its coefficients from the lm
summarize_var <- function(var, base_model, data, trans = c("none", "log", "invnorm"))
{
  trans <- match.arg(trans)
  var <- switch(trans,
                none = var,
                log = sprintf("log(%s)", var), # some have 0 values
                invnorm = sprintf("invNorm(%s)", var))
  fm <- lm(paste0(base_model, " + ", var), data = data)
  smf = summary(fm)
  c(coef(summary(fm))[var, ], r.squared = smf$r.squared)
}

## a version that uses the residuals from the baseline model
## so the model is much simpler
summarize_varR = function (var, data, trans = c("none", "log", "invnorm")) 
{
    trans <- match.arg(trans)
    var <- switch(trans, none = var, log = sprintf("log(%s)", 
        var), invnorm = sprintf("invNorm(%s)", var))
    fm <- lm(paste0("Residuals ~", var), data = data)
    smf = summary(fm)
    c(coef(summary(fm))[var, ], r.squared = smf$r.squared)
}

add_rownames <- function(d, var = "rownames")
{
  d[[var]] <- rownames(d)
  d
}

#' Generic EnWAS 
#'
#' @param data data.frame that the base_model was fit to and all exposures
#' @param base_model character representation of the base_model
#' @param expvars vector of the names of the exposures to be analyzed
#' @param useResiduals logical vector indicating whether or not to fit the model to the residuals from the base_model fit
#'
#' @details
#' This function examines the effects of the exposures listed in `expvars` in a model that is described by the `base_model`.
#' It uses an internal function to perform per variable regressions and returns a data.frame with summary statistics for 
#' each of those fits. If `useResiduals=TRUE` then the base model is fit, and the residuals obtained, each regression is
#' then against those residuals. If `useResiduals=FALSE` then for each exposure the model is augmented by that exposure and
#' fit.
#'
#' @return A data.frame where each row consists of a set of summary statistics for each exposure variable.
#'  The estimates from the model, upper and lower confidence bounds (UCL and LCL), BH adjusted p-values and
#' the multiple Rsquared value from the regression that was performed for that variable (either whole model or
#' against the residuals).
#'
generic_enwas <- function(data, base_model, expvars, ..., useResiduals=FALSE)
{
  if( useResiduals) {
     Residuals = lm(base_model, data = data)$residuals
     data$Residuals = Residuals
    EV_summary_invnorm <-
     sapply(expvars, summarize_varR, 
     data = data, ...) |> t() |> as.data.frame() |> add_rownames("EV") |>
    dplyr::mutate(LCL = `Estimate` - 2 * `Std. Error`,
                  UCL = `Estimate` + 2 * `Std. Error`,
                  adj_pval = p.adjust(`Pr(>|t|)`, "BH"))
  } else {
  EV_summary_invnorm <- 
    sapply(expvars, summarize_var, base_model = base_model, data = data, ...) |>
    t() |> as.data.frame() |> add_rownames("EV") |>
    dplyr::mutate(LCL = `Estimate` - 2 * `Std. Error`,
                  UCL = `Estimate` + 2 * `Std. Error`,
                  adj_pval = p.adjust(`Pr(>|t|)`, "BH"))
    }
  return(EV_summary_invnorm)
 }

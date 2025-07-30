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
  ## print(fm)
  coef(summary(fm))[var, ]
}

## a version that uses the residuals from the baseline model
## so the model is much simpler
summarize_varR = function (var, data, trans = c("none", "log", "invnorm")) 
{
    trans <- match.arg(trans)
    var <- switch(trans, none = var, log = sprintf("log(%s)", 
        var), invnorm = sprintf("invNorm(%s)", var))
    fm <- lm(paste0("Residuals ~", var), data = data)
    coef(summary(fm))[var, ]
}

add_rownames <- function(d, var = "rownames")
{
  d[[var]] <- rownames(d)
  d
}

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



#' Environment‐Wide Association Study (EnWAS)
#'
#' @param base_model a string of base model
#' @param exposure_vars the phenotype list
#' @param data_set data set
#' @param trans "log" :log and scale transformation or inv: inverse normal transformation
#'
#' @return the model list and EnWAS result
#' @export
#'
#' @examples linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
#' linear_res <- enwas(linear_model, exposure_vars, data)
enwas <-
  function(base_model,
           exposure_vars,
           data_set,
           trans = "none") {

    num_var <- length(exposure_vars)

    model_list <- vector(mode = "list", length = num_var)
    names(model_list) <- exposure_vars
    num_cols <- sapply(data_set[, exposure_vars], is.numeric)
    num_cols <- names(num_cols[num_cols == TRUE])

    qc_cols <- c('terms',"logLik","AIC","BIC", # basic info form the model
                 "Deviance","Ratio","P_Chi", # for ANOVA LRT
                 "LR","p_LRT") # Vuong's test
    qc_mtx <- matrix(0, nrow = num_var, ncol = length(qc_cols))
    colnames(qc_mtx) <- qc_cols
    qc_mtx[,1] <- exposure_vars


    association_list <- data.frame()
    # factor_terms <- c() # to hold the factors terms
    for (i in 1:num_var) {
      exposure <- exposure_vars[i]

      # if (is.factor(data_set[, exposure])) {
      #   factor_terms <-
      #     c(factor_terms, paste0(exposure, levels(data_set[, exposure])))
      # }

      # print(sapply(data[!is.na(data[,exposure]),], function(x) sum(is.na(x))))

      # We need to remove the NAs for current phenotype before invNorm(),
      # otherwise it fill NA, which may undesired values.
      ph_dat <- data_set[!is.na(data_set[,exposure]),]
      if (trans=="inv" & exposure %in% num_cols) {
        ph_dat[,exposure] <- invNorm(ph_dat[,exposure])
      }

      if (trans=="log" & exposure %in% num_cols) {
        ph_dat[,exposure] <- scale(log(ph_dat[,exposure]+1e-5))
      }

      base_m <- lm(formula = as.formula(base_model),ph_dat)

      model <- build_formula(base_model, exposure)
      mod <- lm(model, ph_dat)
      model_list[[i]] <- mod
      mod_df <- broom::tidy(mod)
      mod_df$count <- nrow(ph_dat)

      association_list <- rbind(association_list, mod_df)

      qc_mtx[i,2:4] <- round(unlist(broom::glance(mod))[c('logLik','AIC','BIC')],3)
      #-----------------ANOVA LRT---------------------------
      ano_res <- anova(base_m,mod,test="LRT")
      qc_mtx[i,5] <- round(unlist(ano_res[2,"Sum of Sq"]),3)
      qc_mtx[i,6] <- round(ano_res$"Sum of Sq"[2]/ano_res$"RSS"[2]*100,3)
      qc_mtx[i,7] <- round(unlist(ano_res[2,"Pr(>Chi)"]),4)
      #-----------------vuongtest---------------------------
      vong <- nonnest2::vuongtest(base_m,mod)
      qc_mtx[i,8] <- round(vong$LRTstat,3)
      qc_mtx[i,9] <- round(vong$p_LRT$A,4)


    }

    # print(knitr::kable(association_list))
    # xwas_list <-
    #   association_list[association_list$term %in% c(exposure_vars, factor_terms), ]

    xwas_list <-
      association_list[association_list$term %in% exposure_vars, ]


    if(trans =="none"){
      sd_x_list <-  sapply(data_set, function(x)
        sd(as.numeric(x),na.rm = TRUE))
      for (var in exposure_vars) {
        xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] <-
          xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] * sd_x_list[var]
      }

    }



    xwas_list$fdr <-
      p.adjust(xwas_list$p.value, method = 'BY')

    xwas_list$lower <- xwas_list$estimate - 1.96 * xwas_list$std.error
    xwas_list$upper <- xwas_list$estimate + 1.96 * xwas_list$std.error

    qc_mtx <- as.data.frame(qc_mtx)
    qc_mtx <- qc_mtx[qc_mtx$terms %in%xwas_list$term,]
    qc_mtx[,2:9] <- sapply(qc_mtx[,2:9],as.numeric)

    return (list(qc_mtx = qc_mtx,enwas_res = xwas_list,model_list=model_list[xwas_list$term]))



  }



#' Environment‐Wide Association Study (EnWAS) with testing
#'
#' @param base_model a string of base model
#' @param exposure_vars the phenotype list
#' @param train_set taring set
#' @param test_set test set
#' @param lab_col label column
#' @param inv_norm
#'
#' @return the model list and EnWAS result
#' @export
#'
#' @examples linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
#' linear_res <- enwas_cv(linear_model, exposure_vars, train_set,test_set)
enwas_cv <-
  function(base_model,
           exposure_vars,
           train_set,
           test_set,
           lab_col="DIASTOLIC",
           inv_norm = FALSE) {

    num_var <- length(exposure_vars)

    model_list <- vector(mode = "list", length = num_var)
    names(model_list) <- exposure_vars
    # num_cols <- sapply(train_set[, exposure_vars], is.numeric)
    # num_cols <- names(num_cols[num_cols == TRUE])
    if (inv_norm) {
      if (length(exposure_vars)>1){
        train_set[, exposure_vars] <- lapply(train_set[, exposure_vars], invNorm)
        test_set[, exposure_vars] <- lapply(test_set[, exposure_vars], invNorm)
      }else{
        train_set[, exposure_vars] <- invNorm(train_set[, exposure_vars])
        test_set[, exposure_vars] <- invNorm(test_set[, exposure_vars])
        }

    }

    mse <- rep(NA,num_var)

    for (i in 1:num_var) {
      model_str <- build_formula(base_model, exposure_vars[i])
      model <- lm(model_str, train_set,na.action = na.omit)
      pred_vals <- predict(model, test_set)
      mse[i] <- mean((pred_vals - test_set[,lab_col]) ^ 2,na.rm = TRUE)
    }

    mse

  }




#' Build formula for EnWAS
#'
#' @param base_model the formula of the model with the covariates
#' @param exposure_var exposure variables (aka. phenotype lists or predictors)
#' @param inv Boolean flag for use inverse normal transformation
#'
#' @return formula for EnWAS
#' @export
#'
#' @examples build_formula("BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI", phenotype_lists)

build_formula <- function(base_model,exposure_var,inv=FALSE) {
  var <- "+ %s"
  if(inv) var <- " + invNorm(%s)"

  as.formula(
    sprintf(paste0(base_model, var),
            exposure_var
    )
  )
}

#' Inverse Normal Transformation
#'
#' @param x a numeric, complex, character or logical vector.
#' @param na.last for controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed; if "keep" they are kept with rank NA.
#' @param ties.method a character string specifying how ties are treated, it could be "average", "first", "last", "random", "max", "min"
#'
#' @return transformed data
#' @export
#'
#' @examples invNorm(nhanes$BMXWAIST)
invNorm <- function(x,na.last = "keep", ties.method="average") {
  qnorm((rank(x,na.last = na.last, ties.method=ties.method)
         - 3/8)/(length(x) +1 - 6/8))
  }




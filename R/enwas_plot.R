



#' Forest Plot for single model
#'
#' The top n records are selected by the absolute values of the estimators.
#'
#' @param xwas_result the EnWAS result data frame
#' @param top_n the top to n phenotype
#'
#' @return plot result
#' @export
#'
#' @examples forest_plot(linear_res$enwas_res)

forest_plot <- function(xwas_result,top_n=20) {


  xwas_top <- xwas_result |> filter(lower*upper > 0) |>
    dplyr::top_n(top_n,abs(estimate))
  if (nrow(xwas_top)==0){
    xwas_top <- xwas_result
  }

  n <- nrow(xwas_top)
  xwas_top <- xwas_top |> dplyr::arrange(dplyr::desc(estimate)) |>
    dplyr::mutate(xmin = seq(0.5, n - 0.5, by = 1),
                  xmax = seq(1.5, n + 0.5, by = 1))
  xwas_top$col <- as.numeric(rownames(xwas_top)) %% 2

  xwas_top |> ggplot(aes(x = reorder(term,estimate),
               y = estimate,
               colour = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.5,
                  cex = 1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = +Inf,
        fill = factor(col)
      ),
      color = 'black',
      alpha = 0.3
    ) +
    scale_fill_manual(values = c('white', 'grey78'), guide = 'none') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") + ylab("Estimate (95% CI)") +
    theme_bw()  # use a white background
}






#' Plot Forest for Multiple Models
#'
#' @param xwas_result_list list of the EnWAS result data frames
#' @param top_n select top n record by the differences of the estimates
#'
#' @return plot multiple model results in a single image
#' @export
#'
#' @examples forest_plot_mult(list(linear = lm_enwas$enwas_res,
#'                                ns = en_enwas$enwas_res))
forest_plot_mult <- function(xwas_result_list,top_n=20) {

  xwas_result <- do.call("rbind", xwas_result_list)
  xwas_result$EnWAS <-
    rep(names(xwas_result_list), each = nrow(xwas_result_list[[1]]))

  tem_df <- xwas_result_list[[1]]
  terms <- tem_df$term
  len <- length(terms)
  difference <- rep(0,len)
  names(difference) <- terms
  for (i in 1:len){
    difference[i] <- diff(range(xwas_result[xwas_result$term==terms[i],]$estimate))
  }
  top_diff <- names(sort(difference,decreasing=TRUE)[1:top_n])

  tem_df <- tem_df[tem_df$term %in% top_diff,]

  n <- nrow(tem_df)
  tem_df$col <- as.numeric(rownames(tem_df)) %% 2
  tem_df <- tem_df |> dplyr::mutate(xmin = seq(0.5, n - 0.5, by = 1),
                             xmax = seq(1.5, n + 0.5, by = 1))
  xwas_result[xwas_result$term %in% top_diff,] |>
    ggplot(aes(x = reorder(term,estimate),
                             y = estimate,
                             colour = EnWAS))  +
    geom_point(size = 2, position = position_dodge(width = 1)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 1),
      width = 0.5,
      cex = 1
    ) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_rect(
      data = tem_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = +Inf,
        fill = factor(col)
      ),
      color = 'black',
      alpha = 0.3
    ) +
    scale_fill_manual(values = c('white', 'grey78'), guide = 'none') +
    scale_color_manual(
      values = c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7"
      )
    )+
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") +  ylab("Standardized Effect Estimate")+
    theme_bw()  # use a white background
}



#' Make Bins for y values based on x
#'
#' @param x x values
#' @param y y values
#' @param nbin expected number of data point in each bin
#'
#' @return a data frame with binned data frame
#' @export
#'
#' @examples make_bins(x,y,200)
make_bins <- function(x, y, nbin) {
  if (nbin == 0)
    nbin <- floor(sqrt(length(x)))

  df <-
    data.frame(y = y,
               breaks = Hmisc::cut2(x, m = nbin, levels.mean = TRUE))
  df <- df |> na.omit() |>
    dplyr::group_by(breaks) |>
    dplyr::summarize(mean = mean(y),
                     std = sd(y),
                     cnt = dplyr::n()) |> dplyr::mutate(y_min = mean - 1.96 * std /
                                                   sqrt(cnt),
                                                 y_max = mean + 1.96 * std / sqrt(cnt))

  df$breaks <- as.numeric(as.character(df$breaks))

  df

}


#' Plot binned value for single data frame
#'
#' The y values are binned with respected to x values
#'
#' @param x x values
#' @param y y values
#' @param data data frame
#' @param xlab x label text
#' @param ylab y label text
#' @param title plot title
#' @param nobsperbin expected number of data points in per bins
#'
#' @return plot results
#' @export
#'
#' @examples demo = nhanes("DEMO_C")
#' plot_bins(x=RIDAGEYR,y=INDFMPIR,data = demo)
plot_bins<- function(x,y,data,xlab="Value",ylab="Binned y",title="linear",nobsperbin=600){
  arguments <- as.list(match.call())
  y = eval(arguments$y, data)
  x = eval(arguments$x, data)
  df = make_bins(x,y,nobsperbin)
  g = ggplot(df,aes(breaks,mean)) + geom_point() +
    geom_errorbar(aes(ymin=y_min,ymax=y_max))+
    geom_smooth(aes(breaks,mean),method = "lm", formula = y ~  ns(x, df=7))+
    ylab(ylab) + xlab(xlab) + labs(title = title)

  g


}





#' Plot multiple Bins
#'
#' @param df_list a list contains data frames
#' @param xlab x-axis label
#' @param ylab x-axis label
#' @param is_facet Boolean flag to decide create facets
#'
#' @return plots binned data
#' @export
#'
#' @examples plot_bins2(df_list)
plot_bins2 <-
  function(df_list,
           xlab = "x",
           ylab = "Binned y",
           is_facet = TRUE) {
    df <- do.call("rbind", df_list)
    df$Model <-
      rep(names(df_list), each = nrow(df_list[[1]]))

    g <-
      ggplot(df, aes(breaks, mean, color = Model)) + geom_point(size=1.5) +
      geom_errorbar(aes(ymin = y_min, ymax = y_max)) +
      geom_smooth(aes(breaks, mean), method = "loess",
                  formula=y~x,
                  se = FALSE) +
      scale_color_manual(
        values = c(
          "#E69F00",
          "#56B4E9",
          "#009E73",
          "#F0E442",
          "#0072B2",
          "#D55E00",
          "#CC79A7"
        )
      ) +
      ylab(ylab) + xlab(xlab) +
      theme_minimal()
    if (is_facet) {
      g <- g + facet_grid(cols = vars(Model))
    }

    g

  }


#' Print ANOVA LRT test results
#'
#' @param anov_obj result of ANOVA LRT test result
#'
#' @return print the results in a readable format
#' @export
#'
#' @examples print_anova(anova(lm_base,ns_base,test="LRT"))
print_anova <- function(anov_obj){
  for (a in attr(anov_obj,"heading")){
    cat(gsub(pattern = "\n", replacement = "  \n", x = a))
  }
  anov_df <- as.data.frame(anov_obj)
  anov_df <- round(anov_df,3)
  anov_df[!is.na(anov_df[,'Pr(>Chi)']) & anov_df[,'Pr(>Chi)'] <= 0.001, 'Pr(>Chi)'] <- "< 0.001"
  anov_df[is.na(anov_df)] <- ""

  knitr::kable(anov_df)
}


#' Plot P values
#'
#' @param xwas_result_list list of data frames of EnWAS
#' @param is_fdr a flag value decides to plot p-values or FDR
#'
#' @return plot results
#' @export
#'
#' @examples plot_p(list(linear = lm_enwas$enwas_res,
#'                                ns = ns_enwas$enwas_res))
plot_p <- function(xwas_result_list,is_fdr=FALSE){

  xwas_result <- do.call("rbind", xwas_result_list)
  xwas_result$EnWAS <-
    rep(names(xwas_result_list), each = nrow(xwas_result_list[[1]]))

  p <- if (is_fdr) "fdr" else "p.value"
  g_p.value <- xwas_result|> ggplot(aes_string("term", p ,color="EnWAS")) +
    geom_segment(aes_string(x="term", xend="term", y=0, yend="p.value")) +
    geom_point(size=4)+geom_hline(yintercept = 0.05,linetype="dashed")+
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )

  g_p.value

}



#' Plot the outcomes against terms
#'
#' @param data data
#' @param x terms or phenotype
#' @param y outcomes
#' @param gender gender of the cohort for plot point color
#' @param xlab x label text
#' @param ylab y label text
#' @param sample_ratio sample points to plot points, plot all the points is not helpful and time consuming when data set is large
#'
#' @return plot results
#' @export
#'
#' @examples g_raw(nhanes)
#' @examples g_raw(nhanes,x="RIDAGEYR",xlab="Age (Years)")
g_raw <- function(data,
                  x = "BMXBMI",
                  y = "DIASTOLIC",
                  gender = "RIAGENDR",
                  xlab = "BMI (kg/m²)",
                  ylab = "Diastolic (mm Hg)",
                  sample_ratio = 0.3) {
  g <- nhanes |> ggplot(aes_string(x = x,
                                   y = y)) +
    geom_point(
      data = ~ dplyr::group_by(.x) |> dplyr::sample_frac(sample_ratio),
      aes_string(fill = gender),
      alpha = 0.2,
      shape = 21
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      aes(colour = "Linear"),
      size = 1.5
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ splines::ns(x, df = 7),
      aes(colour = "Spline"),
      size = 1.5
    ) +
    xlab(xlab) + ylab(ylab) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9")) +
    theme_minimal() +
    labs(fill = "Gender", colour = "Model")

  g
}


#' Plot for EnWAS AC/AC
#'
#' @param data AC/AC matrix
#' @param y metrics
#' @param x term
#' @param top_n top n terms
#' @param is_desc flag to decide desc
#'
#' @return plot results
#' @export
#'
#' @examples lollipop(ns_enwas$qc_mtx,y="p_LRT")
lollipop <- function(data,
                     y = "AIC",
                     x = "terms",
                     top_n = 20,
                     is_desc = FALSE) {
  x_order <-
    if (is_desc)
      paste0("reorder(", x, ",dplyr::desc(", y, "))")
  else
    paste0("reorder(", x, ",", y, ")")
  data[order(data[, y], decreasing = is_desc), ] |> head(top_n) |>
    ggplot(aes_string(x = x_order, y = y, colour = x)) +
    geom_point(size = 4) +
    geom_segment(aes_string(
      x = x,
      xend = x,
      y = paste0("min(", y, ")"),
      yend = y
    )) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + xlab(x) +
    theme_minimal() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}


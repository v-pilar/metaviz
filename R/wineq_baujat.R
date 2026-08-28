#'@title Baujat plot for the change of influence on the overall effect of each
#'study between a fixed-effect and random-effects model

#'@description Creates a baujat plot which shows the direction and magnitude of
#'change in influence of each study on the overall effect between a fixed-effect
#'and random-effects model based on the same data.

#'@param x metafor rma.uni object conducted with method “FE” or “REML” (the chosen
#'method of the input model makes no difference in the resulting plot)
#'@param labs a vector of study labels
#'@param nr_labs specifies the number of studies to be labelled starting from the
#'righthand side of the x-axis.
#'@param col boolean argument that determines whether the study markers are colored or not.
#'@param method determines whether x-axis values are based on the fixed-effect
#'(“FE”) or the random-effects model (“REML”). If no value is entered, the method
#'is extracted from model x.
#'@param SPR boolean argument that determines whether the x-axis shows the
#'squared pearson residual for the random-effects model (when method == “REML”)
#'instead of the contribution to the Cochran Q-test for heterogeneity to account
#'for differences in the study variances.

#'@details The function wineq_baujat creates a baujat plot (Baujat et al., 2002)
#'with the added feature that influence on the overall result is shown for both
#'the fixed-effect and the random-effects model connected by an arrow. This allows
#'the user to identify both direction and magnitude of the change in a study’s
#'influence on the overall result between the two models. Subsequently, one can
#'easily identify studies that most strongly drive a shift in the overall result
#'between the models. These are the studies that both gain a lot in influence in
#'the random effects model as well as contribute a lot to heterogeneity
#'@return A baujat plot which contains y-axis values for both a fixed-effect and
#'a random-effects model connected by an arrow is created using ggplot2.
#'@author Verena Pilar <verena.pilar@outlook.com>
#'@references Baujat, B., Mahé, C., Pignon, J. - P. & Hill, C. (2002). A graphical
#'method for exploring heterogeneity in meta-analyses: Application to a meta-analysis
#'of 65 trials. \emph{Statistics in Medicine}, 21(\emph{18}): 2641–2652.

#'@examples
#' library(metafor)
#'
#' # Calculating a random-effects model based on the mozart data
#' mozart_r <- rma(yi = d,
#'                 sei = se,
#'                 data = mozart,
#'                 method = "REML")
#'
#' # Plotting the wineq baujat plot based on the mozart data
#' # using a matrix as input
#' wineq_baujat(x = mozart[, c("d", "se")])
#' # using a rma.uni model as input
#' wineq_baujat(x = mozart_r)
#'
#' # Adding study name labels to the right-most studies on the x-axis and
#' # determining how many studies shall be labeled
#' wineq_baujat(x = mozart_r, labs = mozart$study_name, nr_labs = 5)

#'@export

wineq_baujat <- function (x, labs = NULL, nr_labs = 10,
                          col = TRUE,
                          method = NULL,
                          SPR = TRUE) {

  #'@import ggplot2
  #'@import dplyr
  #'@import metaviz
  #'@import metafor
  #'@importFrom magrittr %>%
  NULL

  # preprocessing

  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }

  if("rma" %in% class(x)) {
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    n <- length(es)


    if(is.null(method) || !(method %in% c("REML", "FE", "REM", "FEM"))) {
      method <- x$method
    }

    if(method == "REM") {
      method <- "REML"
    } else if (method == "FEM") {
      method <- "FE"
    }

    if(any(!x$not.na)) {
      warning("The dataset contains at least one missing value, only complete cases are used.")
      if(!is.null(labs)) {
        if(NROW(labs) == nrow(x$data)) {    # using NROW to account for labs possibly being a vector or a data frame
          labs <- labs[x$not.na]
        } else {
          warning("Length of the labs argument doesn't match the length of the dataset. Argument is ignored")
          labs <- NULL
        }

      } else {
        labs <- NULL
      }
    }

  } else {

    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        if(!is.null(labs)) {
          labs <- labs[stats::complete.cases(x[, c(1, 2)])]
        }

        x <- x[stats::complete.cases(x), ]
      }
      # check if input is numeric
      if(!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
        stop("Input argument has to be numeric.")
      }
      # check if there are any negative standard errors
      if(!all(x[, 2] >= 0)) {
        stop("Negative standard errors supplied")
      }
      # extract effects and standard errors
      es <- x[, 1]
      se <- x[, 2]
      n <- length(es)
    } else {
      stop("Unknown input argument. See help ('metaviz').")
    }

    if(is.null(method) || !(method %in% c("REML", "FE"))) {
      method <- "REML"
    }
  }

  # creating df to use as data for the models
  df <- data.frame (
    es = es,
    se = se
  )


  # calculating models

  model_fem <- rma(yi = es,
                   sei = se,
                   data = df,
                   method = "FE")

  model_rem <- rma(yi = es,
                   sei = se,
                   data = df,
                   method = "REML")


  # calculating x-values

  beta_fem <- summary(model_fem)$beta[1]
  beta_rem <- summary(model_rem)$beta[1]



  # x-values
  if(SPR == FALSE) {
    x_fem <- (((model_fem$yi) - c(beta_fem))^2) / model_fem$vi
    x_rem <- (((model_rem$yi) - c(beta_rem))^2) / model_rem$vi
    x_lab <- "Contribution to Overall Heterogeneity"
  } else {
    x_fem <- stats::resid(model_fem)^2 / (model_fem$tau2 + model_fem$vi)
    x_rem <- stats::resid(model_rem)^2 / (model_rem$tau2 + model_rem$vi)
    x_lab <- "Squared Pearson Residual"
  }



  # y-values for REM
  y_rem <- numeric(n)
  y_fem <- numeric(n)

  for (row in 1:nrow(df)) {
    data_int <- df[-row,]      # intermediate data
    model_rem_int <- rma(yi = es,
                         sei = se,
                         data = data_int,
                         method = "REML")
    beta_rem_int <- summary(model_rem_int)$beta[1]
    var_rem_int <- (summary(model_rem_int)$se)^2

    model_fem_int <- rma(yi = es,
                         sei = se,
                         data = data_int,
                         method = "FE")
    beta_fem_int <- summary(model_fem_int)$beta[1]
    var_fem_int <- (summary(model_fem_int)$se)^2


    y_fem_int <- ((beta_fem_int - beta_fem)^2)/var_fem_int
    y_rem_int <- ((beta_rem_int - beta_rem)^2)/var_rem_int
    y_fem[row] <- y_fem_int
    y_rem[row] <- y_rem_int
  }



  # creating plotdata dataframe
  plotdata <- data.frame(
    x_fem = as.numeric(x_fem),
    x_rem = as.numeric(x_rem),
    y_fem = as.numeric(y_fem),
    y_rem = as.numeric(y_rem)
  )

  if (!is.null(labs) && length(labs) == n) {
    plotdata <- plotdata %>%
      mutate(labs = labs)
  } else if (!is.null(labs)) {
    stop("Error: 'labs' must correspond to the number of effect sizes.")
  }


  # creating column for later labs plotting placement (above or below arrow)
  # basis is y_vals_fem (= start of arrow), plus/minus a certain value
  plotdata$labs_pos <- ifelse(plotdata$y_rem <= plotdata$y_fem, -0.8, 1.5)


  x_model <- NULL # to avoid no visible binding for global variable note

  if (method == "REML") {
    plotdata$x_model <- x_rem
  } else if (method == "FE") {
    plotdata$x_model <- x_fem
  } else {
    plotdata$x_model <- x_rem
  }

  max_labs <- plotdata %>%
    arrange(desc(x_model)) %>%
    utils::head(nr_labs)

  # creating plot
  p <- ggplot(data=plotdata, aes(x=x_model))+

    geom_segment(aes(x = x_model, xend = x_model,
                     y = y_fem, yend = y_rem,
                     color = ifelse(y_rem > y_fem, "greater", "smaller")),
                 linewidth=.7,
                 arrow=arrow(angle=30, length=unit(.2, "cm")))

  if (!is.null(labs)) {
    p <- p + geom_text(data = max_labs, aes(y=y_fem, label=labs), vjust=max_labs$labs_pos, size=3)
  }

  if (col == FALSE) {
    p <- p + scale_color_manual(values=c("greater"="grey30", "smaller"="grey70"))
  } else {
    p <- p + scale_color_manual(values=c("greater"="red3", "smaller"="steelblue3"))
  }

  p <- p +
    labs(y = expression("Influence on Overall Effect (FEM " %->% " REM)"), x = x_lab) +
    theme_minimal() +
    theme(legend.position="none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    coord_cartesian(clip = "off")


  p



}



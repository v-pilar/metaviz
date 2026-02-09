
#'@title Forest plot for the comparison of fixed-effect and random-effects models
#'
#'@description Creates a rain, thick or classic forest plot variant that integrates
#'study weights and overall effects of both a fixed-effect and a random-effects model based on the same data.
#'
#'
#'@param x metafor rma.uni object conducted with method “FE” or “REML” (the chosen method of the input model makes
#'no difference in the resulting plot)
#'@param group factor indicating the group membership of each study
#'@param variant “classic” (default), “thick” or “rain” to create a classic, thick or rainforest plot variant
#'@param method determines whether x-axis values are based on the fixed-effect (“FE”) or the random-effects model
#'(“REML”; default). If no value is entered, the method is extracted from model x.
#'@param study_labels y-axis labels for the study effects
#'@param summary_label_FE y-axis label for the fixed-effect model summary effect
#'@param summary_label_REML y-axis label for the random-effects model summary effect
#'@param confidence_level confidence level for the study effects and summary effects
#'@param summary_line adds dashed vertical lines intersecting each summary effect
#'@param summary_col_FE determines color of the fixed-effect model summary effect
#'@param summary_col_REML determines color of the random-effects model summary effect
#'@param col “weights”: colors point estimates (classic variant), errorbars (thick variant) or raindrops (rain variant)
#'according to weight change on a gradient from blue to red, “BW”: greyscale version
#'@param errorbar_col boolean argument that determines whether the errorbars of the classic variant are colored
#'according to weight change (“TRUE”) or in black (“FALSE”)
#'@param text_size determines text size within the plot
#'@param xlab x-axis label
#'@param x_limit determines the limits of the x-axis. Input is a numeric vector of length 2 (min, max).
#'@param x_trans_function function which transforms x-axis labels back to their original scale when data consists,
#'for example, of log-odds-ratios or Fisher’s z values
#'@param x_breaks option to costumize the number of breaks on the x-axis. Input is a numeric vector specifying the breaks
#'@param annotate_CI adds a right-hand side table to the plot containing the confidence intervals of each effect
#'@param study_table custom table on the left-hand side of the plot that contains study information. Takes a dataframe
#'as input which has to be of a length equal to the number of studies.
#'@param summary_table custom table on the left-hand side of the plot that contains information about the summary
#'effects. Takes a dataframe as input which contains one row for each summary effect.
#'@param table_headers headers for each column of the left-hand side table. Takes a character vector as input
#'@param table_layout numeric layout matrix to customize the arrangement of the plot and tables
#'@param show_legend shows a color legend below the plot which corresponds to the change in weight between the fixed and random-effects model
#'@param ... further arguments passed to the internal helper functions for the classic, thick and rainforest variants of the wineq forest plot
#'
#'@details The function wineq_forest creates a forest plot by use of ggplot2 that integrates a fixed-effect
#'and a random-effects model based on the same data. It thus provides insight into a sample-specific as well
#'as generalizable result (Borenstein et al., 2010). Color-coded point estimates denote a gain or loss in
#'study weight between the two models and facilitate the identification of small-study effects. Overall effects
#'are shown for both models. This forest plot comes in three variants: classic forest plot, thick forest plot
#'and rainforest plot (Schild & Voracek, 2015). Optional tables provide additional statistical information
#'about the sample.
#'Note: This function was developed on the basis of the viz_forest code which was created by Michael Kossmeier.
#'@references
#'Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2010).
#'A basic introduction to fixed‐effect and random‐effects models for meta‐analysis.
#'\emph{Research synthesis methods}, 1(\emph{2}), 97-111. https://doi.org/10.1002/jrsm.12
#'
#'Schild, A. H., & Voracek, M. (2015). Finding your way out of the
#'forest without a trail of bread crumbs: Development and evaluation of two
#'novel displays of forest plots. \emph{Research Synthesis Methods}, \emph{6},
#'74-86.
#'@return A wineq forest plot is created by use of ggplot2.
#'@author Verena Pilar <verena.pilar@univie.ac.at>
#'@examples
#' library(metafor)
#' # Arranging the data according to effect size to faciliate the identification of small-study effects
#' mozart <- mozart[order(mozart$d),]
#'
#' # Calculating a random-effects model based on the mozart data
#' mozart_r <- rma(yi = d,
#'                 sei = se,
#'                 data = mozart,
#'                 method = "REML")
#'
#' # Plotting a wineq forest plot based on the mozart data
#' # using a matrix as input
#' wineq_forest(x = mozart[, c("d", "se")], study_labels = mozart$study_name)
#' # using a rma.uni model as input
#' wineq_forest(x = mozart_r, study_labels = mozart$study_name)
#'
#' # Thick and rainforest plot variants of the wineq forest plot
#' wineq_forest(mozart_r, variant = "thick")
#' wineq_forest(mozart_r, variant = "rain")

#'@export
wineq_forest <- function(x, group = NULL, variant = "classic", method = "REML",
                         study_labels = NULL, summary_label_FE = NULL, summary_label_REML = NULL,
                         confidence_level = 0.95, summary_line = TRUE,
                         summary_col_FE = "grey70", summary_col_REML = "grey50",
                         col = "weights", errorbar_col = TRUE,
                         text_size = 3, xlab = "Effect", x_limit = NULL,
                         x_trans_function = NULL, x_breaks = NULL,
                         annotate_CI = FALSE, study_table = NULL, summary_table = NULL,
                         table_headers = NULL, table_layout = NULL, show_legend = FALSE, ...) {

  #'@import ggplot2
  #'@import dplyr
  #'@import metafor

  # Handle input object -----------------------------------------------------
  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }


  if("rma" %in% class(x)) {
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    n <- length(es)


    # check if group argument has the right length
    if(!is.null(group) & (length(group) != length(es))) {
      warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored.")
      group <- NULL
    }


    # If No group is supplied try to extract group from input object of class rma.uni (metafor)
    if(is.null(group) && ncol(x$X) > 1) {
      #check <- if only categorical moderators were used
      if(!all(x$X == 1 || x$X == 0) || any(apply(as.matrix(x$X[, -1]), 1, sum) > 1))  {
        stop("Can not deal with metafor output object with continuous and/or more than one categorical moderator variable(s).")
      }
      # extract group vector from the design matrix of the metafor object
      no.levels <- ncol(x$X) - 1
      group <- factor(apply(as.matrix(x$X[, -1])*rep(1:no.levels, each = n), 1, sum))
    }
  } else {
    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        study_labels <- study_labels[stats::complete.cases(x[, c(1, 2)])]
        if(!is.null(group)) {
          group <- group[stats::complete.cases(x)]
        }
        x <- x[stats::complete.cases(x), ]
      }
      # check if input is numeric
      if(!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
        stop("Input argument has to be numeric; see help(viz_forest) for details.")
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
  }

  # Preprocess data ---------------------------------------------------------
  # check if group is a factor
  if(!is.null(group) && !is.factor(group)) {
    group <- as.factor(group)
  }
  # check if group vector has the right length
  if(!is.null(group) && (length(group) != length(es))) {
    warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored")
    group <- NULL
  }

  # if no group argument is supplied, use all cases
  if(is.null(group)) {
    group <- factor(rep(1, times = n))
  }

  # drop unused levels of group factor
  group <- droplevels(group)
  k <- length(levels(group))

  # main data
  x <- data.frame(es, se, group)


  # check col is of length 1, or nrow(x) in case of variant classic or thick
  if(variant == "rain") {
    stopifnot(length(col) == 1)
  } else {
    if(variant == "thick" || variant == "classic") {
      stopifnot(length(col) == 1 || length(col) == nrow(x))
    }
  }

  # check summary_col_FE is of length 1, or length(levels(group)) in case of variant classic or thick
  if(variant == "rain") {
    stopifnot(length(summary_col_FE) == 1)
  } else {
    if(variant == "thick" || variant == "classic") {
      stopifnot(length(summary_col_FE) == 1)
    }
  }

  # check summary_col_REML is of length 1, or length(levels(group)) in case of variant classic or thick
  if(variant == "rain") {
    stopifnot(length(summary_col_REML) == 1)
  } else {
    if(variant == "thick" || variant == "classic") {
      stopifnot(length(summary_col_REML) == 1)
    }
  }



  # Compute meta-analytic summary effect estimates
  M <- NULL # To avoid "no visible binding for global variable" warning for non-standard evaluation
  # compute meta-analytic summary effect for each group
  get_bse <- function(es, se, type = "b") {
    res_f <- metafor::rma.uni(yi = es, sei = se, method = "FE")
    res_r <- metafor::rma.uni(yi = es, sei = se, method = "REML")
    if(type == "b") {
      res_b <- c(res_f$b[[1]],res_r$b[[1]])
    } else {
      if(type == "se") {
        res_se <- c(res_f$se[[1]], res_r$se[[1]])
      } else {
        stop()
      }
    }
  }
  M <- x %>%
    group_by(group) %>%
    summarise(M_FE = get_bse(es, se, type = "b")[1],
              M_REML = get_bse(es, se, type = "b")[2])

  M_FE <- M %>%
    select(M_FE)
  M_REML <- M %>%
    select(M_REML)
  summary_es_FE <- unlist(M_FE)
  summary_es_REML <- unlist(M_REML)


  # compute standard error of the meta-analytic summary effect for each group
  M <- x %>%
    group_by(group) %>%
    summarise(M_FE = get_bse(es, se, type = "se")[1],
              M_REML = get_bse(es, se, type = "se")[2])
  M_FE <- M %>%
    select(M_FE)
  M_REML <- M %>%
    select(M_REML)
  summary_se_FE <- unlist(M_FE)
  summary_se_REML <- unlist(M_REML)



  # Compute tau^2 estimate
  # compute tau squared for each group
  get_tau2 <- function(es, se) {
    metafor::rma.uni(yi = es, sei = se, method = "REML")$tau2[[1]]
  }


  M <- x %>%
    group_by(group) %>%
    summarise(M = get_tau2(es, se)) %>%
    select(M)
  summary_tau2_REML <- unlist(M)
  summary_tau2_FE <- rep(0, times = k)


  # if not exactly one name for every study is supplied the default is used (numbers 1 to the number of studies)
  if(is.null(study_labels) || length(study_labels) != n) {
    if(!is.null(study_labels) && length(study_labels) != n) {
      warning("Argument study_labels has wrong length and is ignored.")
    }
    study_labels <- 1:n
  }

  # if not exactly one name for every subgroup is suppied the default is used
  if(is.null(summary_label_FE) || length(summary_label_FE) != k || is.null(summary_label_REML) || length(summary_label_REML) != k) {
    if(!is.null(summary_label_REML) && length(summary_label_REML) != k) {
      warning("Argument summary_label_REML has wrong length and is ignored.")
    }
    if(!is.null(summary_label_FE) && length(summary_label_FE) != k) {
      warning("Argument summary_label_FE has wrong length and is ignored.")
    }
    if(k != 1) {
      summary_label_FE <- paste("Subgroup: ", levels(group), " FEM", sep = "")
      summary_label_REML <- paste("Subgroup: ", levels(group), " REM", sep = "")
    } else {
      summary_label_FE <- "Summary FEM"
      summary_label_REML <- "Summary REM"

    }
  }

  if(confidence_level <= 0 || confidence_level >= 1) {
    stop("Argument confidence_level must be larger than 0 and smaller than 1.")
  }

  if(!is.null(x_trans_function) && !is.function(x_trans_function)) {
    warning("Argument x_trans_function must be a function; input ignored.")
    x_trans_function <- NULL
  }

  if(!is.null(table_layout) && !is.matrix(table_layout)) {
    warning("Agument of table_layout is not a matrix and is ignored.")
    table_layout <- NULL
  }

  # Determine IDs for studies and summary effects which correspond to plotting y coordinates
  ids <- function(group, n) {
    k <- length(levels(group))
    ki_start <- cumsum(c(1, as.numeric(table(group))[-k] + 3))
    ki_end <- ki_start + as.numeric(table(group)) -1
    study_IDs <- numeric(n)
    for(i in 1:k) {
      study_IDs[group==levels(group)[i]] <- ki_start[i]:ki_end[i]
    }

    summary_IDs_FE <- ki_end + 2
    summary_IDs_REML <- ki_end + 3


    data.frame("ID" = (n + 3*k-2) - c(study_IDs, summary_IDs_FE, summary_IDs_REML),
               "type" = factor(c(rep("study", times = length(study_IDs)),
                                 rep ("summary_FE", times = length(summary_IDs_FE)),
                                 rep ("summary_REML", times = length(summary_IDs_REML)))))


  }

  ID <- ids(group, n = n)




  plotdata <- data.frame("x" = es, "se" = se,
                         "ID" = ID$ID[ID$type == "study"],
                         "labels" = study_labels,
                         "group"= group,
                         "x_min" = es - stats::qnorm(1 - (1 - confidence_level)/2)*se,
                         "x_max" = es + stats::qnorm(1 - (1 - confidence_level)/2)*se)


  madata <- data.frame("summary_es_FE" = summary_es_FE,
                       "summary_es_REML" = summary_es_REML,
                       "summary_se_FE" = summary_se_FE,
                       "summary_se_REML" = summary_se_REML,
                       "summary_tau2_FE" = summary_tau2_FE,
                       "summary_tau2_REML" = summary_tau2_REML,
                       "ID_FE" = ID$ID[ID$type == "summary_FE"],
                       "ID_REML" = ID$ID[ID$type == "summary_REML"])






  # Create forest plot variant ------------------------------------------------------
  args <- c(list(plotdata = plotdata, madata = madata, col = col,
                 method = method, summary_line = summary_line,
                 study_labels = study_labels, summary_label_FE = summary_label_FE,
                 summary_label_REML = summary_label_REML,
                 study_table = study_table, summary_table = summary_table,
                 annotate_CI = annotate_CI, confidence_level = confidence_level,
                 summary_col_FE = summary_col_FE, summary_col_REML = summary_col_REML,
                 errorbar_col = errorbar_col,
                 text_size = text_size, xlab = xlab, x_limit = x_limit,
                 x_trans_function = x_trans_function, x_breaks = x_breaks), list(...))

  if(variant == "rain") {
    p <- do.call(internal_wineq_forest_rain, args)
  } else {
    if(variant == "thick") {
      p <- do.call(internal_wineq_forest_thick, args)
    } else {
      if(variant == "classic") {
        p <- do.call(internal_wineq_forest_classic, args)
      } else {
        stop("The argument of variant must be one of rain, thick or classic.")
      }
    }
  }




  #####################################################
  # TABLEPLOTS
  #####################################################

  #### create RAIN FP legend ##########################################################################
  if (variant == "rain") {
    legbreaks <- c(0, 0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9, 2.2, 2.5, 3)
    legbreaks_labs <- c("0   ", 0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9, 2.2, 2.5, "3<")

    leg_colors <- c("#00c0e0", "#2be0ff", "#c7faff", "#dbd9d9", "#dbd9d9", "#c7adb7", "#ad5c6a", "#a82342", "#70020f", "#000000")

    # Create data frame for segments
    leg_df <- data.frame(
      xmin = legbreaks[-length(legbreaks)],
      xmax = legbreaks[-1],
      ymin = 0,
      ymax = 1,
      fill = leg_colors
    )

    # Data for labels at break points
    leg_labels_df <- data.frame(
      x = legbreaks,
      y = -0.1,
      label = as.character(legbreaks_labs)
    )

    leg <- ggplot() +
      geom_rect(data = leg_df,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), color = "#ffffff"
      ) +
      scale_fill_identity() +
      geom_text(data = leg_labels_df, aes(x = x, y = y, label = label), vjust = 1, size = 2.5) +
      scale_y_continuous(expand = c(0,0), limits = c(-0.3, 1.2)) +
      coord_fixed(ratio = 0.4, clip = "off") +
      theme_void() +
      theme(legend.position = "none")



    # make titel its own grob
    legend_title <- grid::textGrob(
      "Relative weight change\nloss < 1.0 > gain",
      x = 0, y = 0.5,
      just = "right",
      hjust = 0,
      gp = grid::gpar(fontsize = 10)
    )

    # gtable to combine the two
    legend_table <- gtable::gtable(
      widths = grid::unit.c(
        grid::unit(4, "cm"),
        grid::unit(0.5, "npc")
      ),
      heights = grid::unit(1.5, "npc")
    )

    # add the grobs to the table
    legend_table <- gtable::gtable_add_grob(
      legend_table,
      legend_title,
      t = 1, l = 1
    )

    legend_table <- gtable::gtable_add_grob(
      legend_table,
      ggplotGrob(leg),
      t = 1, l = 2
    )


    leg_grob <- legend_table
  }


  ####################################################

  # Construct tableplots with study and summary information --------
  if(annotate_CI == TRUE || !is.null(study_table) || !is.null(summary_table)) {

    # set limits for the y axis of the table plots
    y_limit <- c(min(plotdata$ID) - 4, max(plotdata$ID) + 1.5)


    # Function to create table plots
    table_plot <- function(tbl, ID, r = 5.5, l = 5.5, tbl_titles = NULL) {
      # all columns and column names are stacked to a vector
      df_to_vector <- function(df) {
        v <- vector("character", 0)
        for(i in 1:ncol(df)) v <- c(v, as.vector(df[, i]))
        v
      }
      if(!is.data.frame(tbl)) tbl <- data.frame(tbl)
      tbl <- data.frame(lapply(tbl, as.character), stringsAsFactors = FALSE)
      if(is.null(tbl_titles)) {
        tbl_titles <- names(tbl)
      }
      v <- df_to_vector(tbl)

      # For study labels with newlines in it, the width of the column is now set according to longest line and not the whole label
      nchar2<-function(x){unlist(sapply(strsplit(x,"\n"), function(x) max(nchar(x, keepNA = FALSE))))}
      area_per_column <- cumsum(c(1, apply(rbind(tbl_titles, tbl), 2, function(x) max(round(max(nchar2(x))/100, 2),  0.03))))

      x_values <- area_per_column[1:ncol(tbl)]
      x_limit <- range(area_per_column)


      lab <- data.frame(y = rep(ID, ncol(tbl)),
                        x = rep(x_values,
                                each = length(ID)),
                        value = v, stringsAsFactors = FALSE)

      lab_title <- data.frame(y = rep(max(plotdata$ID) + 1, times = length(tbl_titles)),
                              x = x_values,
                              value = tbl_titles)

      # To avoid "no visible binding for global variable" warning for non-standard evaluation
      y <- NULL
      value <- NULL
      ggplot(lab, aes(x = x, y = y)) +
        geom_text(aes(label = value), size = text_size, hjust = 0, vjust = 0.5) +
        geom_text(data = lab_title, aes(x = x, y = y, label = value), size = text_size, hjust = 0, vjust = 0.5) +
        coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F) +
        geom_hline(yintercept = max(plotdata$ID) + 0.5) +
        theme_bw() +
        theme(text = element_text(size = 1/0.352777778*text_size),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              axis.text.x = element_text(colour="white"),
              axis.text.y = element_blank(),
              axis.ticks.x = element_line(colour="white"),
              axis.ticks.y = element_blank(),
              axis.line.x = element_line(colour="white"),
              axis.line.y = element_blank(),
              plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt")) +
        labs(x = "", y = "")
    }

    # Study and/or summary table left
    if(!is.null(study_table) || !is.null(summary_table)) {
      # Case study table and summary table are both supplied
      if(!is.null(study_table) && !is.null(summary_table)) {
        if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
        if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
        study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
        summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
        if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
        if(nrow(summary_table) != k*2) stop('summary_table must be a data.frame with one row for each summary effect.')
        if(ncol(summary_table) < ncol(study_table)) {
          n_fillcol <- ncol(study_table) - ncol(summary_table)
          summary_table <- data.frame(summary_table, matrix(rep("", times = nrow(summary_table) * n_fillcol), ncol = n_fillcol))
          summary_table<- stats::setNames(summary_table, names(study_table))
        } else {
          if(ncol(summary_table) > ncol(study_table)) {
            n_fillcol <- ncol(summary_table) - ncol(study_table)
            study_table <- data.frame(study_table, matrix(rep("", times = nrow(study_table) * n_fillcol), ncol = n_fillcol))
            study_table <- stats::setNames(study_table, names(summary_table))
          }
        }
        if(any(names(study_table) != names(summary_table))) summary_table <- stats::setNames(summary_table, names(study_table))
      } else {
        # Case only study table is supplied
        if(is.null(summary_table)) {
          if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
          study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
          if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
          summary_table <- as.data.frame(matrix(rep("", times = ncol(study_table) * k), ncol = ncol(study_table)), stringsAsFactors = FALSE)
          summary_table <- stats::setNames(summary_table, names(study_table))

        }
        # Case only summary table is supplied
        if(is.null(study_table)) {
          if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
          summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
          if(nrow(summary_table) != k*2) stop('summary_table must be a data.frame with one row for each summary effect.')
          study_table <- as.data.frame(matrix(rep("", times = ncol(summary_table) * n), ncol = ncol(summary_table)), stringsAsFactors = FALSE)
          study_table <- stats::setNames(study_table, names(summary_table))
        }
      }

      table_left <- data.frame(rbind(study_table, summary_table))


      # set table headers
      if(!is.null(table_headers)) {
        if(length(table_headers) >= ncol(table_left)) {
          table_headers_left <- table_headers[1:ncol(table_left)]
        } else {
          warning("Argument table_headers has not the right length and is ignored.")
          table_headers_left <- NULL
        }
      } else {
        table_headers_left <- NULL
      }

      table_left_plot <- table_plot(table_left, ID = ID$ID, r = 0, tbl_titles = table_headers_left)
    } else {
      table_left <- NULL
    }

    # Textual CI and effect size values right
    if(annotate_CI == TRUE) {

      # set table headers
      if(!is.null(table_headers)) {
        if(is.null(table_left)) {
          if(length(table_headers) == 1) {
            table_headers_right <- table_headers
          } else {
            warning("Argument table_headers has not the right length and is ignored.")
            table_headers_right <- NULL
          }
        } else {
          if(length(table_headers) == ncol(table_left) + 1) {
            table_headers_right <- table_headers[ncol(table_left) + 1]
          } else {
            table_headers_right <- NULL
          }
        }
      } else {
        table_headers_right <- NULL
      }

      if(is.null(table_headers_right)){
        table_headers_right <- paste(xlab, " [", confidence_level*100, "% CI]", sep = "")
      }

      x_hat <- c(plotdata$x, madata$summary_es_FE, madata$summary_es_REML)
      lb <- c(c(plotdata$x, madata$summary_es_FE, madata$summary_es_REML) - stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, madata$summary_se_FE, madata$summary_se_REML))
      ub <-  c(c(plotdata$x, madata$summary_es_FE, madata$summary_es_REML) + stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, madata$summary_se_FE, madata$summary_se_REML))

      if(!is.null(x_trans_function)) {
        x_hat <- x_trans_function(x_hat)
        lb <- x_trans_function(lb)
        ub <- x_trans_function(ub)
      }

      lb <- format(round(lb, 2), nsmall = 2)
      ub <- format(round(ub, 2), nsmall = 2)
      x_hat <- format(round(x_hat, 2), nsmall = 2)

      CI <- paste(x_hat, " [", lb, ", ", ub, "]", sep = "")
      CI_label <- data.frame(CI = CI, stringsAsFactors = FALSE)

      table_CI <- table_plot(CI_label, ID = c(plotdata$ID, madata$ID_FE, madata$ID_REML), l = 0, r = 11,  tbl_titles = table_headers_right)

    } else {
      table_CI <- NULL
    }



    ######################################################################################
    # Align forest plot and table(s)
    if(!is.null(table_CI) && !is.null(table_left)) {
      if(is.null(table_layout)) {
        layout_matrix <- matrix(c(rep(1, times = ncol(table_left)), rep(2, times = 3), 3), nrow = 1)

      } else {
        layout_matrix <- table_layout
      }

      if (variant == "rain") {
        legend_grob <- leg_grob


      } else{
        get_legend <- function(p) {
          tmp <- ggplotGrob(p)
          leg <- gtable::gtable_filter(tmp, "guide-box")
          return(leg)
        }

        legend_grob <- get_legend(p)
        p <- p + theme(legend.position = "none")
      }

      main_grob <- gridExtra::arrangeGrob(
        table_left_plot,
        p,
        table_CI,
        layout_matrix = layout_matrix
      )

      main_grob_gt <- gtable::gtable_add_rows(
        gtable::gtable_add_grob(
          gtable::gtable_add_rows(main_grob, grid::unit(1, "lines")),
          grobs = grid::grobTree(legend_grob),
          t = nrow(main_grob) + 1,
          l = ncol(table_left) + 1,
          r = ncol(main_grob) - 1
        ),
        grid::unit(0.2, "null")
      )

      if (show_legend==FALSE) {
        p <- main_grob
      } else {
        p <- main_grob_gt
      }

      ggpubr::as_ggplot(p)


    } else {
      if(!is.null(table_CI) && is.null(table_left)) {
        if(is.null(table_layout)) {
          layout_matrix <- matrix(c(1, 1, 1, 1, 2), nrow = 1)

        } else {
          layout_matrix <- table_layout
        }

        if (variant == "rain") {
          legend_grob <- leg_grob

        } else{
          get_legend <- function(my_plot) {
            tmp <- ggplotGrob(my_plot)
            leg <- gtable::gtable_filter(tmp, "guide-box")
            return(leg)
          }

          legend_grob <- get_legend(p)
          p <- p + theme(legend.position = "none")

        }

        main_grob <- gridExtra::arrangeGrob(
          p,
          table_CI,
          layout_matrix = layout_matrix
        )

        main_grob_gt <- gtable::gtable_add_rows(
          gtable::gtable_add_grob(
            gtable::gtable_add_rows(main_grob, grid::unit(2, "lines")),
            grobs = grid::grobTree(legend_grob),
            t = nrow(main_grob) + 1,
            l = 1,
            r = 4
          ),
          grid::unit(0.2, "null")
        )


        if (show_legend==FALSE) {
          p <- main_grob
        } else {
          p <- main_grob_gt
        }

        ggpubr::as_ggplot(p)


      } else {
        if(is.null(table_CI) && !is.null(table_left)) {
          if(is.null(table_layout)) {
            layout_matrix <- matrix(c(rep(1, times = 1 + ncol(table_left)), 2, 2, 2, 2, 2), nrow = 1)
          } else {
            layout_matrix <- table_layout
          }


          if (variant == "rain") {
            legend_grob <- leg_grob

          } else{
            get_legend <- function(my_plot) {
              tmp <- ggplotGrob(my_plot)
              leg <- gtable::gtable_filter(tmp, "guide-box")
              return(leg)
            }

            legend_grob <- get_legend(p)
            p <- p + theme(legend.position = "none")
          }

          main_grob <- gridExtra::arrangeGrob(
            table_left_plot,
            p,
            layout_matrix = layout_matrix
          )

          main_grob_gt <- gtable::gtable_add_rows(
            gtable::gtable_add_grob(
              gtable::gtable_add_rows(main_grob, grid::unit(1, "lines")),
              grobs = grid::grobTree(legend_grob),
              t = nrow(main_grob) + 1,
              l = ncol(table_left) + 1,
              r = ncol(main_grob)
            ),
            grid::unit(0.2, "null")
          )


          if (show_legend==FALSE) {
            p <- main_grob
          } else {
            p <- main_grob_gt
          }

          ggpubr::as_ggplot(p)
        }
      }
    }
  } else {
    if (variant == "rain") {
      if (show_legend == FALSE) {
        p
      } else {

        legend_grob <- leg_grob


        layout_matrix <- matrix(c(rep(1,5)), nrow = 1)

        main_grob <- gridExtra::arrangeGrob(
          p,
          layout_matrix = layout_matrix
        )

        main_grob_gt <- gtable::gtable_add_rows(
          gtable::gtable_add_grob(
            gtable::gtable_add_rows(main_grob, grid::unit(3, "lines")),
            grobs = grid::grobTree(legend_grob),
            t = nrow(main_grob) + 1,
            l = 2,
            r = ncol(main_grob)
          ),
          grid::unit(0.2, "null")
        )


        p <- main_grob_gt
        ggpubr::as_ggplot(p)
      }

    } else {
      if (show_legend == FALSE) {
        p + theme(legend.position = "none")
      } else {
        p
      }
    }

  }
}



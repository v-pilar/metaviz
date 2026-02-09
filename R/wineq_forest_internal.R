#'Internal helper function of wineq_forest to create a classic wineq forest plot.
#'
#'Creates a classic wineq forest plot. Called by wineq_forest for variant = "classic"
#'Note: This function was developed on the basis of the viz_forest_internal code which was created by Michael Kossmeier.
#'@keywords internal
internal_wineq_forest_classic <- function(plotdata, madata, summary_line = NULL,
                                          method = "REML",
                                          study_labels = NULL, summary_label_FE = NULL,
                                          summary_label_REML = NULL,
                                          study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                          confidence_level = 0.95, #col = "Blues",
                                          col = NULL, summary_col_FE = "steelblue1",
                                          summary_col_REML = "steelblue4",tick_col = "firebrick",
                                          errorbar_col = TRUE, #PI = TRUE,
                                          text_size = 3, xlab = "Effect", x_limit = NULL,
                                          x_trans_function = NULL, x_breaks = NULL) {

  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))



  weight_REML <- 1/(plotdata$se^2 + madata$summary_tau2_REML[as.numeric(plotdata$group)])
  weight_FE <- 1/plotdata$se^2

  plotdata$rel_weight_REML <- weight_REML/sum(weight_REML)
  plotdata$rel_weight_FE <- weight_FE/sum(weight_FE)
  plotdata$rel_wc <- plotdata$rel_weight_REML/plotdata$rel_weight_FE


  if(method == "FE") {
    plotdata$plot_meth <- plotdata$rel_weight_FE
  } else if (method == "REML") {
    plotdata$plot_meth <- plotdata$rel_weight_REML
  } else {
    stop("Plotting_method has to be one of 'FE' and 'REML'")
  }


  # set limits and breaks for the y axis and construct summary diamond

  y_limit <- c(min(plotdata$ID) - 4, max(plotdata$ID) + 1.5)
  y_tick_names <- c(as.vector(study_labels), as.vector(summary_label_FE), as.vector(summary_label_REML))[order(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)]
  y_breaks <- sort(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)
  summarydata1 <- data.frame("x.diamond" = c(madata$summary_es_FE - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE,
                                             madata$summary_es_FE + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE),
                             "y.diamond" = c(madata$ID_FE,
                                             madata$ID_FE  + 0.3,
                                             madata$ID_FE,
                                             madata$ID_FE - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )
  summarydata2 <- data.frame("x.diamond" = c(madata$summary_es_REML - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML,
                                             madata$summary_es_REML + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML),
                             "y.diamond" = c(madata$ID_REML,
                                             madata$ID_REML + 0.3,
                                             madata$ID_REML,
                                             madata$ID_REML - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )


  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }

  # Set plot margins. If table is aligned on the left, no y axis breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  ID <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL

  # create classic forest plot

  plotdata$errorbar_color <- ifelse(plotdata$rel_wc > 1, "#75001b", "#8ab2b8")
  plotdata$errorbar_BW <- ifelse(plotdata$rel_wc > 1, "grey30", "grey70")

  p <-
    ggplot(data = plotdata, aes(y = ID, x = x, fill = rel_wc)) +
    geom_vline(xintercept = 0, linetype = 2)
  if(errorbar_col == TRUE) {
    if(col=="weights") {
      p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                   col = errorbar_color))
    } else {
      p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                   col = errorbar_BW))
    }
  } else {
    p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                 col = "black"))
  }

  p <- p + scale_color_identity()  # added so ggplot2 interprets hex color codes as the colors and not factor levels



  p <- p + geom_point(data = plotdata, aes(size = plot_meth, fill = rel_wc), shape = 22, col = "black")



  if(summary_line == TRUE) {
    p <- p + geom_vline(xintercept = madata$summary_es_FE, linetype = 2, linewidth = 0.5, color="grey70")+
      geom_vline(xintercept = madata$summary_es_REML, linetype = 2, linewidth = 0.5, color="grey50")
  }
  p <- p + geom_polygon(data = summarydata1, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_FE, linewidth = 0.1)
  p <- p + geom_polygon(data = summarydata2, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_REML, linewidth = 0.1)

  p <- p +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
  if(!is.null(x_trans_function)) {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)})
    } else {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)},
                           breaks = x_breaks)
    }
  } else {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab)
    } else {
      p <- p +
        scale_x_continuous(breaks = x_breaks,
                           name = xlab)
    }
  }

  if (is.null(col)) {
    col <- "weights"
  }

  if(col=="weights"){
    p <- p +
      scale_fill_gradientn(
        colors = c("#2be0ff","#a6f2ff", "#c7eaf0","#ededed", "#e6e6e6","#de8e8e","#b8001f","#450008", "#1a0003"),
        values = scales::rescale(c(0, 0.25, 0.50, 0.75, 1, 1.45, 1.75, 2.05, 2.35)),
        limits = c(0,3),
        n.breaks = 6,
        guide = "colorbar",
        name = "Relative weight change\nloss < 1.0 < gain",
        oob = scales::squish  # Ensure out-of-bounds values take the extreme colors
      )

  } else if (col=="BW") {
    p <- p +
      scale_fill_gradientn(n.breaks = 6, colors = c("#ffffff","#e6e6e6", "#9c9c9c","#575757","#000000"))
  }else {
    p <- p +
      scale_fill_gradientn(n.breaks = 6, colors = c(rep(col,5)))
  }

  p <- p +
    guides(size = "none") +
    scale_size_area(max_size = 3) +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))
  p

}














#'Internal helper function of wineq_forest to create a thick wineq forest plot.
#'
#'Creates a thick wineq forest plot. Called by wineq_forest for variant = "thick"
#'Note: This function was developed on the basis of the viz_forest_internal code which was created by Michael Kossmeier.
#'@keywords internal
internal_wineq_forest_thick <- function(plotdata, madata, summary_line = NULL,
                                        method = "REML",
                                        study_labels = NULL, summary_label_FE = NULL,
                                        summary_label_REML = NULL,
                                        study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                        confidence_level = 0.95, #col = "Blues",
                                        col = NULL, summary_col_FE = "steelblue1",
                                        summary_col_REML = "steelblue4",tick_col = "firebrick",
                                        errorbar_col = TRUE,
                                        text_size = 3, xlab = "Effect", x_limit = NULL,
                                        x_trans_function = NULL, x_breaks = NULL) {

  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))

  weight_REML <- 1/(plotdata$se^2 + madata$summary_tau2_REML[as.numeric(plotdata$group)])
  weight_FE <- 1/plotdata$se^2
  plotdata$rel_weight_REML <- weight_REML/sum(weight_REML)
  plotdata$rel_weight_FE <- weight_FE/sum(weight_FE)
  plotdata$rel_wc <- plotdata$rel_weight_REML/plotdata$rel_weight_FE


  if (method == "REML") {
    rel_weight <- weight_REML/sum(weight_REML)
  } else {
    rel_weight <- weight_FE/sum(weight_FE)
  }
  plotdata$rel_weight <- rel_weight
  plotdata <- plotdata %>%
    mutate(y_max = ID + rel_weight/(4*max(rel_weight)),
           y_min = ID - rel_weight/(4*max(rel_weight))
    )


  # set limits and breaks for the y axis and construct summary diamond

  y_limit <- c(min(plotdata$ID) - 4, max(plotdata$ID) + 1.5)
  y_tick_names <- c(as.vector(study_labels), as.vector(summary_label_FE), as.vector(summary_label_REML))[order(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)]
  y_breaks <- sort(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)
  summarydata1 <- data.frame("x.diamond" = c(madata$summary_es_FE - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE,
                                             madata$summary_es_FE + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE),
                             "y.diamond" = c(madata$ID_FE,
                                             madata$ID_FE  + 0.3,
                                             madata$ID_FE,
                                             madata$ID_FE - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )
  summarydata2 <- data.frame("x.diamond" = c(madata$summary_es_REML - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML,
                                             madata$summary_es_REML + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML),
                             "y.diamond" = c(madata$ID_REML,
                                             madata$ID_REML + 0.3,
                                             madata$ID_REML,
                                             madata$ID_REML - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )



  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }


  # Set plot margins. If table is aligned on the left, no y axus breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL
  y_min <- NULL
  y_max <- NULL
  ID <- NULL


  ### error bar color
  plotdata$errorbar_color <- ifelse(plotdata$rel_wc > 1, "#75001b", "#8ab2b8")
  plotdata$errorbar_BW <- ifelse(plotdata$rel_wc > 1, "grey30", "grey70")
  plotdata$point_color <- ifelse(plotdata$rel_wc > 1, "firebrick", "#76c7de")
  plotdata$point_BW <- ifelse(plotdata$rel_wc > 1, "firebrick", "#76c7de")

  # tickdata
  tick_size <- max(plotdata$rel_weight/(6*max(plotdata$rel_weight)))
  tickdata <- data.frame(x = c(plotdata$x, plotdata$x), ID = c(plotdata$ID, plotdata$ID),
                         y = c(plotdata$ID + tick_size,
                               plotdata$ID - tick_size),
                         point_color = c(plotdata$point_color, plotdata$point_color),
                         point_BW = c(plotdata$point_BW, plotdata$point_BW))

  # Create thick forest plot
  p <-
    ggplot(data = plotdata, aes(y = ID, x = x, fill = rel_wc))

  if(errorbar_col == TRUE) {
    if(col=="weights") {
      p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                   col = errorbar_color))
    } else {
      p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                   col = errorbar_BW))
    }
  } else {
    p <- p + geom_errorbarh(data = plotdata, aes(xmin = x_min, xmax = x_max, y = ID, width = 0,
                                                 col = "black"))
  }

  p <- p + scale_color_identity()

  p <- p + geom_rect(data = plotdata, aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max,
                                          group = ID), linewidth = 0.1)



  if(col=="weights") {
    p <- p + geom_line(data = tickdata, aes(x = x, y = y, group = ID, color = point_color),
                       linewidth = 1.5, inherit.aes = FALSE)
  } else {
    p <- p + geom_line(data = tickdata, aes(x = x, y = y, group = ID, color = point_BW),
                       linewidth = 1.5, inherit.aes = FALSE)
  }




  if(summary_line == TRUE) {
    p <- p + geom_vline(xintercept = madata$summary_es_FE, linetype = 2, linewidth = 0.5, color="grey70")+
      geom_vline(xintercept = madata$summary_es_REML, linetype = 2, linewidth = 0.5, color="grey50")
  }
  p <- p + geom_polygon(data = summarydata1, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_FE, linewidth = 0.1)
  p <- p + geom_polygon(data = summarydata2, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_REML, linewidth = 0.1)


  p <- p +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
  if(!is.null(x_trans_function)) {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)})
    } else {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)},
                           breaks = x_breaks)
    }
  } else {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab)
    } else {
      p <- p +
        scale_x_continuous(breaks = x_breaks,
                           name = xlab)
    }
  }


  if (is.null(col)) {
    col <- "weights"
  }

  if(col=="weights"){
    p <- p +
      scale_fill_gradientn(
        colors = c("#2be0ff","#a6f2ff", "#c7eaf0","#ededed", "#e6e6e6","#de8e8e","#b8001f","#450008", "#1a0003"),
        values = scales::rescale(c(0, 0.25, 0.50, 0.75, 1, 1.45, 1.75, 2.05, 2.35)),
        limits = c(0,3),
        n.breaks = 6,
        guide = "colorbar",
        name = "Relative weight change\nloss < 1.0 < gain",
        oob = scales::squish  # Ensure out-of-bounds values adopt the extreme colors
      )

  } else if (col=="BW") {
    p <- p +
      scale_fill_gradientn(n.breaks = 6, colors = c("#e6e6e6","#c4c4c4", "#9c9c9c","#575757","#000000"))
  }else {
    p <- p +
      scale_fill_gradientn(n.breaks = 6, colors = c(rep(col,5)))
  }


  p <- p +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          legend.position = "bottom",
          #legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))
  p
}


















#'Internal helper function of wineq_forest to create a wineq rainforest plot.
#'
#'Creates a wineq rainforest plot. Called by wineq_forest for variant = "rain"
#'Note: This function was developed on the basis of the viz_forest_internal code which was created by Michael Kossmeier.
#'@keywords internal
internal_wineq_forest_rain <- function(plotdata, madata, summary_line = NULL,
                                       method = "REML",
                                       study_labels = NULL, summary_label_FE = NULL,
                                       summary_label_REML = NULL,
                                       study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                       confidence_level = 0.95, col = NULL,
                                       summary_col_FE = NULL, summary_col_REML = NULL,
                                       detail_level = 1,
                                       text_size = 3, xlab = "Effect", x_limit = NULL,
                                       x_trans_function = NULL, x_breaks = NULL,
                                       errorbar_col = TRUE) {


  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))

  # weight of each study used to scale the height of each raindrop

  weight_REML <- 1/(plotdata$se^2 + madata$summary_tau2_REML[as.numeric(plotdata$group)])
  weight_FE <- 1/plotdata$se^2

  plotdata$rel_weight_REML <- weight_REML/sum(weight_REML)
  plotdata$rel_weight_FE <- weight_FE/sum(weight_FE)
  plotdata$rel_wc <- plotdata$rel_weight_REML/plotdata$rel_weight_FE



  if (method == "FE"){
    plotdata$rel_weight <- plotdata$rel_weight_FE
  } else {
    plotdata$rel_weight <- plotdata$rel_weight_REML
  }


  plotdata <- plotdata %>%
    mutate(tick_col = ifelse(rel_wc < 1, "firebrick", "#76c7de"))


  tick_size <- max(plotdata$rel_weight/(6 * max(plotdata$rel_weight)))
  tickdata <- data.frame(x = c(plotdata$x, plotdata$x), ID = c(plotdata$ID, plotdata$ID),
                         tick_col = c(plotdata$tick_col, plotdata$tick_col),
                         y = c(plotdata$ID + tick_size,
                               plotdata$ID - tick_size))

  # function ll constructs a likelihood raindrop of a study. Each raindrop is built out of
  # several distinct segments (to color shade the raindrop)
  ll <- function(x, max.range, max.weight) {
    # width of the region over which the raindop is built
    se.factor <- ceiling(stats::qnorm(1 - (1 - confidence_level)/2))
    width <- abs((x[1] - se.factor * x[2]) - (x[1] + se.factor * x[2]))

    # max.range is internally determined as the width of the broadest raindop.
    # The number of points to construct the raindrop is chosen proportional to
    # the ratio of the width of the raindrop and the max.range,
    # because slim/dense raindrops do not need as many support points as very broad ones.
    # Minimum is 200 points (for the case that width/max.range gets very small)
    length.out <- max(c(floor(1000 * width / max.range), 200))

    # Create sequence of points to construct the raindrop. The number of points is chosen by length.out
    # and can be changed by the user with the parameter detail_level.
    # At least 50 support points (for detail level << 1) per raindrop are chosen
    support <- seq(x[1] - se.factor * x[2], x[1] + se.factor * x[2], length.out = max(c(length.out * detail_level, 50)))

    # The values for the likelihood drop are determined: The likelihood for different hypothetical true values
    # minus likelihood for the observed value (i.e. the maximum likelihood) plus the confidence.level quantile of the chi square
    # distribution with one degree of freedom divided by two.
    # Explanation: -2*(log(L(observed)/L(hypothetical))) is an LRT test and approx. chi^2 with 1 df and significance threshold
    # qchisq(confidence.level, df = 1).
    # That means by adding the confidence.level quantile of the chi square
    # distribution with one degree of freedom (the significance threshold) divided by two,
    # values with l_mu < 0 differ significantly from the observed value.
    threshold <- stats::qchisq(confidence_level, df = 1)/2
    l_mu <- log(stats::dnorm(x[1], mean = support, sd = x[2])) - log(stats::dnorm(x[1], mean = x[1], sd = x[2])) + threshold

    #scale raindrop such that it is proportional to the meta-analytic weight and has height smaller than 0.5
    l_mu <- l_mu/max(l_mu) * x[3]/max.weight * 0.45

    # Force raindrops of studies to have minimum height of 0.05 (i.e. approx. one tenth of the raindrop with maximum height)
    if(max(l_mu) < 0.05) {
      l_mu <- l_mu/max(l_mu) * 0.05
    }

    # mirror values for raindrop
    l_mu_mirror <- -l_mu

    # select only likelihood values that are equal or larger than zero,
    # i.e. values that also lie in the confidence interval (using normality assumption)
    sel <- which(l_mu >= 0)

    # Construct data.frame
    d <- data.frame("support" = c(support[sel], rev(support[sel])), "log_density" = c(l_mu[sel], rev(l_mu_mirror[sel])))

    # The number of segments for shading is chosen as follows: 40 segements times the detail_level per drop
    # as default. The minimum count of segments is 20 (If detail_level is << 1), with the exception that
    # if there are too few points for 20 segments then nrow(d)/4 is used (i.e. at least 4 points per segment)
    data.frame(d, "segment" = cut(d$support, max(c(40*detail_level), min(c(20, nrow(d)/4)))))
  }

  # compute the max range of all likelihood drops for function ll.
  max.range <- max(abs((plotdata$x + stats::qnorm(1 - (1 - confidence_level)/2) * plotdata$se) -
                         (plotdata$x - (stats::qnorm(1 - (1 - confidence_level)/2) * plotdata$se))))

  # computes all likelihood values and segments. The output is a list, where every element
  # constitutes one study raindop
  res <- apply(cbind(plotdata$x, plotdata$se, plotdata$rel_weight), 1,  FUN = function(x) {ll(x, max.range = max.range, max.weight = max(plotdata$rel_weight))})
  ### try adding color information here!


  # name every list entry, i.e. raindrop, and add id column
  names(res) <- plotdata$ID
  for(i in 1:length(res)) {
    res[[i]] <- data.frame(res[[i]], .id = plotdata$ID[i])
  }

  # The prep.data function prepares the list of raindrops in three ways for plotting (shading of segments):
  # 1) the values are sorted by segments, such that the same segments of each raindrop are joined together
  # 2) segments are renamed with integer values from 1 to the number of segments per raindrop
  # 3) to draw smooth raindrops the values at the right hand boundary of each segment have to be the first
  # values at the left hand boundary of the next segment on the right.
  prep.data <- function(res) {
    res <- lapply(res, FUN = function(x) {x <- x[order(x$segment), ]})
    res <- lapply(res, FUN = function(x) {x$segment <- factor(x$segment, labels = 1:length(unique(x$segment))); x})
    res <- lapply(res, FUN = function(x) {
      seg_n <- length(unique(x$segment))
      first <- sapply(2:seg_n, FUN = function(n) {min(which(as.numeric(x$segment)==n))})
      last <-  sapply(2:seg_n, FUN = function(n) {max(which(as.numeric(x$segment)==n))})
      neighbor.top <-   x[c(stats::aggregate(support~segment, FUN = which.max, data=x)$support[1],
                            cumsum(stats::aggregate(support~segment, FUN = length, data=x)$support)[-c(seg_n-1, seg_n)] +
                              stats::aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]), c("support", "log_density")]
      neighbor.bottom <-   x[c(stats::aggregate(support~segment, FUN = which.max, data=x)$support[1],
                               cumsum(stats::aggregate(support~segment, FUN = length, data=x)$support[-c(seg_n-1, seg_n)])+
                                 stats::aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]) + 1, c("support", "log_density")]
      x[first, c("support", "log_density")] <- neighbor.top
      x[last, c("support", "log_density")] <- neighbor.bottom
      x
    }
    )
    res
  }
  res <- prep.data(res)

  # merge the list of raindops in one dataframe for plotting
  res <- do.call(rbind, res)

  # set limits and breaks for the y axis and construct summary diamond

  y_limit <- c(min(plotdata$ID) - 4, max(plotdata$ID) + 1.5) ### -3 --> -4
  y_tick_names <-  c(as.vector(study_labels), as.vector(summary_label_FE), as.vector(summary_label_REML))[order(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)] ###
  y_breaks <- sort(c(plotdata$ID, madata$ID_FE, madata$ID_REML), decreasing = T)
  summarydata1 <- data.frame("x.diamond" = c(madata$summary_es_FE - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE,
                                             madata$summary_es_FE + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_FE,
                                             madata$summary_es_FE),
                             "y.diamond" = c(madata$ID_FE,
                                             madata$ID_FE  + 0.3,
                                             madata$ID_FE,
                                             madata$ID_FE - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )
  summarydata2 <- data.frame("x.diamond" = c(madata$summary_es_REML - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML,
                                             madata$summary_es_REML + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se_REML,
                                             madata$summary_es_REML),
                             "y.diamond" = c(madata$ID_REML,
                                             madata$ID_REML + 0.3,
                                             madata$ID_REML,
                                             madata$ID_REML - 0.3),
                             "diamond_group" = rep(1:k, times = 4)
  )


  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }

  # To shade all segments of each raindop symmetrically the min abs(log_density) per raindrop is used
  # as aesthetic to fill the segments. This is necessary because otherwise the first log_density value per
  # segment would be used leading to asymmetrical shading
  min.ld <- stats::aggregate(log_density ~ segment + .id, FUN  = function(x) {min(abs(x))}, data = res)
  names(min.ld) <- c("segment", ".id", "min_log_density")
  res <- merge(res, min.ld, sort = F)


  ### specify thresholds for differently colored raindrops
  # above 1 = red/dark, below 1 = blue/light
  # several intermittent steps

  plotdata$col_group <- NA

  plotdata <- plotdata %>%
    mutate(col_group = case_when(
      rel_wc > 2.5                 ~ "gain5",
      rel_wc > 2.2 & rel_wc <= 2.5 ~ "gain4",
      rel_wc > 1.9 & rel_wc <= 2.2 ~ "gain3",
      rel_wc > 1.6 & rel_wc <= 1.9 ~ "gain2",
      rel_wc > 1.3 & rel_wc <= 1.6 ~ "gain1",
      rel_wc < 0.1                 ~ "loss3",
      rel_wc < 0.4 & rel_wc >= 0.1 ~ "loss2",
      rel_wc < 0.7 & rel_wc >= 0.4 ~ "loss1",
      TRUE                         ~ "mid"
    ))



  # prepare color info for matching with res
  col_info <- plotdata %>%
    select(ID, col_group)

  # create corresponding ID columns to avoid errors due to different col names
  res$ID <- as.numeric(res$.id)
  col_info$ID <- as.numeric(col_info$ID)

  # match res and col_info
  res_col <- left_join(res, col_info, by = c("ID"))

  # for plotting
  res <- res_col



  # Set plot margins. If table is aligned on the left, no y axis breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  support <- NULL
  segment <- NULL
  min_log_density <- NULL
  log_density <- NULL
  .id <-
    x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL
  ID <- NULL



  plotdata$point_color <- ifelse(plotdata$rel_wc > 1, "firebrick", "#76c7de")


  # Create Rainforest plot
  p <-
    ggplot(data = res, aes(y = .id, x = support)) +
    geom_errorbarh(data = plotdata, col = "steelblue2", aes(xmin = x_min, xmax = x_max, y = ID, width = 0), inherit.aes = FALSE) +

    geom_polygon(data = res %>% filter(col_group == "gain5"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#000000", low = "#ba001f", guide = "none") +
    scale_color_gradient(high = "#000000", low = "#ba001f", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "gain4"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#70020f", low = "#e31751", guide = "none") +
    scale_color_gradient(high = "#70020f", low = "#e31751", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "gain3"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#a82342", low = "#d9529a", guide = "none") +
    scale_color_gradient(high = "#a82342", low = "#d9529a", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "gain2"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#ad5c6a", low = "#ccc0c4", guide = "none") +
    scale_color_gradient(high = "#ad5c6a", low = "#ccc0c4", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "gain1"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#c7adb7", low = "#d6ced0", guide = "none") +
    scale_color_gradient(high = "#c7adb7", low = "#d6ced0", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "mid"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#dbd9d9", low = "#fafafa", guide = "none") +
    scale_color_gradient(high = "#dbd9d9", low = "#fafafa", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "loss1"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#c7faff", low = "#daecf5", guide = "none") +
    scale_color_gradient(high = "#c7faff", low = "#daecf5", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "loss2"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#2be0ff", low = "#f2f7ff", guide = "none") +
    scale_color_gradient(high = "#2be0ff", low = "#f2f7ff", guide = "none") +

    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +

    geom_polygon(data = res %>% filter(col_group == "loss3"),
                 aes(x = support, y = as.numeric(.id) + log_density,
                     color = min_log_density,
                     fill = min_log_density,
                     group = paste(.id, segment)), linewidth = 0.1) +
    scale_fill_gradient(high = "#00c0e0", low = "#c7faff", guide = "none") +
    scale_color_gradient(high = "#00c0e0", low = "#c7faff", guide = "none") +

    ggnewscale::new_scale_color() # for point_color in next geom_point




  p <- p + geom_line(data = tickdata, aes(x = x, y = y, group = ID, color = tick_col),
                     linewidth = 1) +
    guides(color = "none")



  if(summary_line == TRUE) {
    p <- p + geom_vline(xintercept = madata$summary_es_FE, linetype = 2, linewidth = 0.5, color="grey70")+
      geom_vline(xintercept = madata$summary_es_REML, linetype = 2, linewidth = 0.5, color="grey50")
  }

  p <- p + geom_polygon(data = summarydata1, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_FE, linewidth = 0.1)
  p <- p + geom_polygon(data = summarydata2, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col_REML, linewidth = 0.1)

  p <- p +

    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
  if(!is.null(x_trans_function)) {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)})
    } else {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)},
                           breaks = x_breaks)
    }
  } else {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab)
    } else {
      p <- p +
        scale_x_continuous(breaks = x_breaks,
                           name = xlab)
    }
  }



  p <- p +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))

  p
}



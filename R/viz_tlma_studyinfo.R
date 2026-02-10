
#'@title A dataframe containing information on the study effects of a three-level meta-analysis
#'@description Provides a dataframe with information on the study effects of a three-level meta-analysis.
#'
#'@param x metafor rma.mv object
#'@param study_ID grouping variable for the study-level of the three-level meta-analysis that was used in the rma.mv model
#'@param effect_ID grouping variable for the effect-level of the three-level meta-analysis that was used in the rma.mv model
#'@param confidence_level_ci numeric confidence level for the confidence intervals of the study effects
#'
#'@details The function viz_tlma_studyinfo creates a dataframe with information pertaining to the study effects of
#'a three-level meta-analysis. The study effects are the result of random-effects models conducted on the effects
#'contained within each respective study. The dataframe contains the number of effects originating from the study,
#'the study effect with its standard error and confidence intervals, and the weight of the study within the three-level
#'meta-analysis.
#'@return A dataframe containing statistical information about the study effects of a three-level meta-analysis is created.
#'@author Verena Pilar <verena.pilar@univie.ac.at>
#'@examples
#' library(metafor)
#' library(psymetadata)
#' # Extract wibbelink2017 data
#' testdata <- wibbelink2017
#'
#' # Determine level 1 and 2 grouping variables for the three-level meta-analysis
#' # study IDs = level 2
#' ID <- testdata$study_id
#' # effect IDs = level 1
#' ID2 <- testdata$es_id
#'
#' # Calculate the three-level meta-analytic model
#' testmodel <- rma.mv(yi,
#'                     vi,
#'                     random = ~ 1 | ID/ID2,
#'                     tdist = TRUE,
#'                     data = testdata,
#'                     method = "REML")
#'
#' # Create a table with study-level information
#' viz_tlma_studyinfo(testmodel, ID, ID2)

#' @export

viz_tlma_studyinfo <- function (x, study_ID, effect_ID, confidence_level_ci = 0.95) {

  group = NULL

  ID <- study_ID
  ID2 <- effect_ID

  n_ID <- max(ID)

  if("rma.mv" %in% class(x)) {
    yi <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    vi <- as.numeric(x$vi)
    n <- length(yi)
    n_ID <- length(unique(ID))

  } else {
    stop("The input has to be a rma.mv model.")
  }

  if(is.null(ID) || is.null(ID2)) {
    stop("Please provide the arguments study_ID and effect_ID")
  } else if (length(ID) != length(yi) || length(ID2) != length(yi)) {
    stop("study_ID and effect_ID have to be of the same length as you model input dataset.")
  }


  if(is.null(group)) {
    group <- as.factor(rep(1, times = length(yi)))
  } else {
    group <- as.factor(group)
  }

  # drop unused levels of group factor
  group <- droplevels(group)
  k <- length(levels(group))

  data <- data.frame(yi, se, vi, ID, ID2, group)


  model <- x

  estimate <- round(model$b[1], 2)
  var_bs <- model$sigma2[1]    # between-studies variance
  var_ws <- model$sigma2[2]  # within-study variance


  # CIs
  data <- data %>%
    mutate(ci_lb = yi - se * qnorm(0.975)) %>%
    mutate(ci_ub = yi + se * qnorm(0.975))

  ##############################################################################
  ##############################################################################
  ####        #########  #####################################  ################
  ####        #########  #####################################  ################
  ####        #########                                         ######      ####
  ####        #########  #####################################  ######      ####
  #####################  #####################################  ######      ####
  ##############################################################################
  ##############################################################################

  # creating a separate dataset for study-level information
  studydata = data.frame(
    yi_ID = numeric(n_ID),
    se_ID = numeric(n_ID),
    k = numeric(n_ID),
    ci_lb_ID = numeric(n_ID),
    ci_ub_ID = numeric(n_ID),
    ci_lb_ES = numeric(n_ID),
    ci_ub_ES = numeric(n_ID),
    weight_ID = numeric(n_ID),
    type = numeric(n_ID)
  )

  # add ID column for merging later on
  studydata$ID <- unique(data$ID)


  row <- 1


  ###############################


  for (i in 1:max(data$ID)){
    subdata<-subset(data, ID==i)
    uni=nrow(subdata)

    if (uni==1) {
      studydata$yi_ID[row] <- subdata$yi
      studydata$se_ID[row] <- subdata$se
      studydata$ci_lb_ID[row] <-  subdata$yi - (subdata$se * 1.96)
      studydata$ci_ub_ID[row] <-  subdata$yi + (subdata$se * 1.96)
      studydata$ci_lb_ES[row] <- subdata$yi - (subdata$se * 1.96)
      studydata$ci_ub_ES[row] <- subdata$yi + (subdata$se * 1.96)
      studydata$weight_ID[row] <- 1 / subdata$se^2
      studydata$type[row] <- "singleES"
    }
    else {
      model_ID <- metafor::rma.uni(yi = subdata$yi, sei = subdata$se, method = "REML", data = subdata)

      diagonal <- 1/(subdata$vi + var_ws)
      D <- diag(diagonal)
      obs <- nrow(subdata)
      I <- matrix(c(rep(1, (obs^2))), nrow = obs)
      M <- D%*%I%*%D
      inv_sumVar <- sum(1/(subdata$vi + var_ws))
      O <- 1/((1/var_bs) + inv_sumVar)
      V <- D - (O*M)
      T <- as.matrix(subdata$yi)
      X <- matrix(c(rep(1, obs)), ncol=1)
      var_effect <- solve(t(X)%*%V%*%X)


      studydata$yi_ID[row] <- model_ID$b
      studydata$se_ID[row] <- model_ID$se
      studydata$ci_lb_ID[row] <-  model_ID$ci.lb
      studydata$ci_ub_ID[row] <-  model_ID$ci.ub
      studydata$ci_lb_ES[row]<-model_ID$b - 1.96*median(subdata$se)
      studydata$ci_ub_ES[row]<-model_ID$b + 1.96*median(subdata$se)
      studydata$weight_ID[row]<-1/ var_effect
      studydata$type[row] <- "multiES"
    }

    studydata$k[row]<-nrow(subdata)
    studydata$J[row] <- c(paste("J =",studydata$k[i]))

    row <-  row + 1
  }


  # arrange studydata for table plotting later
  studydata <- studydata %>%
    arrange(yi_ID)

  studylevel_info <- studydata %>%
    select(
      k = k,
      yi_study = yi_ID,
      se_study = se_ID,
      ci_lb_study = ci_lb_ID,
      ci_ub_study = ci_ub_ID,
      weight_study = weight_ID
    )

  return(studylevel_info)


}

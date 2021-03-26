##################
# Load libraries #
##################

library(SuperLearner)
library(tidyverse)
library(caret)

################
# INPUT FORMAT #
################

# `feature_table` 
# - should be a data frame with features in rows and samples in columns.

# `sample_metadata` 
# - should be a data frame containing sample-specific metadata. Must have a column 
# named 'subjectID' describing per-subject unique identifiers. For longitudinal designs, 
# this variable is expected to have non-unique values. Additionally, a column named 'Y' must be present 
# which is the outcome of interest (can be binary or continuous). 
# Row names of sample_metadata must match the column names of feature_table.

# `feature_metadata`:
# - should be a data frame containing feature-specific metadata. Must have a column 
# named 'featureID' describing per-feature unique identifiers. Additionally, if multiple omics layers 
# are present, a column named 'featureType' should describe the corresponding source layers  
# (e.g. metagenomics, metabolomics, etc.). Row names must match that of feature_table.

##############################
# Integrated Learner Wrapper #
##############################

run_integrated_learner_nested_CV<-function(feature_table,
                                           sample_metadata, 
                                           feature_metadata,
                                           outerCVfolds = 10, # How many folds in the external V-fold CV? Defaul is 10.
                                           innerCVfolds = 0, # How many folds per iteration of the V-fold CV? Default is 0 (no inner CV).
                                           seed = 1234, # Specify the arbitrary seed value for reproducibility. Default is 1234.
                                           base_learner = 'SL.BART', # Base learner for single layers and concatenation.
                                           meta_learner = 'SL.BART', # Meta learner for stacked generalization.
                                           run_concat = TRUE, # Should vanilla concatenated base learner be run? Default is TRUE.
                                           verbose = TRUE # Should detailed message be printed? Default is TRUE.
){ 
                                    
  
  #######################
  # Basic sanity checks #
  #######################
  
  ############################
  # Check dimension mismatch #
  ############################
  
  if(all(rownames(feature_table)==rownames(feature_metadata))==FALSE)
    stop("Both feature_table and feature_metadata should have the same rownames.")
  
  if(all(colnames(feature_table)==rownames(sample_metadata))==FALSE)
    stop("Row names of sample_metadata must match the column names of feature_table.")
  
  #########################
  # Check missing columns #
  #########################
  
  if (!'subjectID' %in% colnames(sample_metadata)){
    stop("sample_metadata must have a column named 'subjectID' describing per-subject unique identifiers.")
  }
  
  if (!'Y' %in% colnames(sample_metadata)){
    stop("sample_metadata must have a column named 'Y' describing the outcome of interest.")
  }
  
  if (!'featureID' %in% colnames(feature_metadata)){
    stop("feature_metadata must have a column named 'featureID' describing per-feature unique identifiers.")
  }
  
  if (!'featureType' %in% colnames(feature_metadata)){
    stop("feature_metadata must have a column named 'featureType' describing the corresponding source layers.")
  }

  ###############################################################
  # Set parameters and extract subject IDs for sample splitting #
  ###############################################################
  
  set.seed(seed)
  subjectID <- unique(sample_metadata$subjectID)
  
  ##################################
  # Trigger V-fold CV (Outer Loop) #
  ##################################
  
  subjectCvFoldsIN <- caret::createFolds(1:length(subjectID), k = outerCVfolds, returnTrain=TRUE)
  
  ########################################
  # Curate subect-level samples per fold #
  ########################################
  
  obsIndexIn <- vector("list", outerCVfolds) 
  for(k in 1:length(obsIndexIn)){
    x <- which(!sample_metadata$subjectID %in%  subjectID[subjectCvFoldsIN[[k]]])
    obsIndexIn[[k]] <- x
  }
  names(obsIndexIn) <- sapply(1:outerCVfolds, function(x) paste(c("fold", x), collapse=''))
  
  ###############################
  # Set up data for SL training #
  ###############################
  
  outerCVcontrol = list(V = outerCVfolds, shuffle = FALSE, validRows = obsIndexIn)

  #################################################
  # Stacked generalization input data preparation #
  #################################################
  
  feature_metadata$featureType<-as.factor(feature_metadata$featureType)
  name_layers<-with(droplevels(feature_metadata), list(levels = levels(featureType)), nlevels = nlevels(featureType))$levels
  SL_fit_predictions<-vector("list", length(name_layers)) 
  SL_fit_layers<-vector("list", length(name_layers)) 
  names(SL_fit_layers)<-name_layers
  names(SL_fit_predictions)<-name_layers
  
  ##################################################################
  # Carefully subset data per omics and run each individual layers #
  ##################################################################
  
  for (i in seq_along(name_layers)){
    
    if (verbose) cat('Running base model for layer ', i, "...", "\n")
    
    ##################################
    # Prepate single-omic input data #
    ##################################
    
    include_list<-feature_metadata %>% filter(featureType == name_layers[i]) 
    t_dat_slice<-feature_table[rownames(feature_table) %in% include_list$featureID, ]
    dat_slice<-as.data.frame(t(t_dat_slice))
    Y = sample_metadata$Y
    X = dat_slice
    
    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_layers[[i]] <- SuperLearner::CV.SuperLearner(Y = Y, 
                                                        X = X, 
                                                        cvControl = outerCVcontrol,    
                                                        innerCvControl = list(list(V = innerCVfolds)),
                                                        control = list(saveFitLibrary = TRUE),
                                                        verbose = verbose, 
                                                        SL.library = base_learner) 
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_layers[[i]]$Y<-sample_metadata['Y']
    SL_fit_layers[[i]]$X<-X
    
    ###############################################################
    # Remove redundant data frames and save pre-stack predictions #
    ###############################################################
    
    rm(include_list); rm(t_dat_slice); rm(dat_slice); rm(X) 
    SL_fit_predictions[[i]]<-SL_fit_layers[[i]]$SL.predict
  }
  
  ####################
  # Stack all models #
  ####################
  
  if (verbose) cat('Running stacked model...\n')
  
  ##############################
  # Prepate stacked input data #
  ##############################
  
  combo <- as.data.frame(do.call(cbind, SL_fit_predictions))
  names(combo)<-name_layers

  ###################################
  # Run user-specified meta learner #
  ###################################
  
  SL_fit_stacked<- SuperLearner::CV.SuperLearner(Y = Y, 
                                                 X = combo, 
                                                 cvControl = outerCVcontrol,    
                                                 innerCvControl = list(list(V = innerCVfolds)),
                                                 control = list(saveFitLibrary = TRUE),
                                                 verbose = verbose, 
                                                 SL.library = meta_learner) 
  
  ###################################################
  # Append the corresponding y and X to the results #
  ###################################################
  
  SL_fit_stacked$Y<-sample_metadata['Y']
  SL_fit_stacked$X<-combo

  #######################################
  # Run concatenated model if specified #
  #######################################
  
  if(run_concat){
    
    if (verbose) cat('Running concatenated model...\n')
    
    ###################################
    # Prepate concatenated input data #
    ###################################
    
    fulldat<-as.data.frame(t(feature_table))

    ###################################
    # Run user-specified base learner #
    ###################################
    
    SL_fit_concat<-SuperLearner::CV.SuperLearner(Y = Y, 
                                                 X = fulldat, 
                                                 cvControl = outerCVcontrol,    
                                                 innerCvControl = list(list(V = innerCVfolds)),
                                                 control = list(saveFitLibrary = TRUE),
                                                 verbose = verbose, 
                                                 SL.library = base_learner) 
    
    ###################################################
    # Append the corresponding y and X to the results #
    ###################################################
    
    SL_fit_concat$Y<-sample_metadata['Y']
    SL_fit_concat$X<-fulldat

    ######################
    # Save model results #
    ######################
    
    SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                     SL_fit_stacked = SL_fit_stacked, 
                     SL_fit_concat = SL_fit_concat)
    } else{ 
      SL_fits<-list(SL_fit_layers = SL_fit_layers, 
                    SL_fit_stacked = SL_fit_stacked)
    }
  
  ##########
  # Return #
  ##########
  
  return(SL_fits)
}  


#####################################################################
# Rename after adding serialize = TRUE in bartMachine2 (from ck37r) #
#####################################################################

# Temporary wrapper that needs to be fixed in SuperLearner
#' Wrapper for bartMachine learner
#'
#' Support bayesian additive regression trees via the bartMachine package.
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Optional dataframe to predict the outcome
#' @param obsWeights Optional observation-level weights (supported but not tested)
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param family "gaussian" for regression, "binomial" for binary
#'   classification
#' @param num_trees The number of trees to be grown in the sum-of-trees model.
#' @param num_burn_in Number of MCMC samples to be discarded as "burn-in".
#' @param num_iterations_after_burn_in Number of MCMC samples to draw from the
#'   posterior distribution of f(x).
#' @param alpha Base hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param beta Power hyperparameter in tree prior for whether a node is
#'   nonterminal or not.
#' @param k For regression, k determines the prior probability that E(Y|X) is
#'   contained in the interval (y_{min}, y_{max}), based on a normal
#'   distribution. For example, when k=2, the prior probability is 95\%. For
#'   classification, k determines the prior probability that E(Y|X) is between
#'   (-3,3). Note that a larger value of k results in more shrinkage and a more
#'   conservative fit.
#' @param q Quantile of the prior on the error variance at which the data-based
#'   estimate is placed. Note that the larger the value of q, the more
#'   aggressive the fit as you are placing more prior weight on values lower
#'   than the data-based estimate. Not used for classification.
#' @param nu Degrees of freedom for the inverse chi^2 prior. Not used for
#'   classification.
#' @param verbose Prints information about progress of the algorithm to the
#'   screen.
#' @param serialize If TRUE, bartMachine results can be saved to a file, but
#'   will require additional RAM.
#' @param ... Additional arguments (not used)
#'
#' @encoding utf-8
#' @export
SL.BART <- function(Y, X, newX, family, obsWeights, id,
                    num_trees = 50, num_burn_in = 250, verbose = F,
                    alpha = 0.95, beta = 2, k = 2, q = 0.9, nu = 3,
                    num_iterations_after_burn_in = 1000,
                    serialize = TRUE,
                    ...) {
  #.SL.require("bartMachine")
  
  ################
  ### CK changes:
  if (family$family == "binomial") {
    # Need to convert Y to a factor, otherwise bartMachine does regression.
    # And importantly, bartMachine expects the first level to be the positive
    # class, so we have to specify levels.
    Y = factor(Y, levels = c("1", "0"))
  }
  model = bartMachine::bartMachine(X, Y, num_trees = num_trees,
                                   num_burn_in = num_burn_in, verbose = verbose,
                                   alpha = alpha, beta = beta, k = k, q = q, nu = nu,
                                   num_iterations_after_burn_in = num_iterations_after_burn_in,
                                   serialize = serialize)
  # pred returns predicted responses (on the scale of the outcome)
  #pred <- bartMachine:::predict.bartMachine(model, newX)
  pred <- predict(model, newX)
  
  fit <- list(object = model)
  class(fit) <- c("SL.bartMachine")
  
  out <- list(pred = pred, fit = fit)
  return(out)
}


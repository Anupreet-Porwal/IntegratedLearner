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

run_integrated_learner_repeatedCV<-
  function(feature_table,
           sample_metadata,
           feature_metadata,
           feature_table_valid = NULL, # Feature table from validation set. Must have the exact same structure as feature_table. If missing, uses feature_table for feature_table_valid.
           sample_metadata_valid = NULL, # Sample-specific metadata table from independent validation set. Must have the exact same structure as sample_metadata.
           folds = 5, # How many folds in the K-fold CV? Defaul is 5.
           reps = 20, # How many repeats of the K-fold CV? Default is 20.
           seed = 1234, # Specify the arbitrary seed value for reproducibility. Default is 1234.
           base_learner = 'SL.BART', # Base learner for single layers and concatenation.
           meta_learner = 'SL.BART', # Meta learner for stacked generalization.
           run_concat = TRUE, # Should vanilla concatenated base learner be run? Default is TRUE.
           run_stacked = TRUE, # Should stacked model be run? Default is TRUE.
           verbose = TRUE # Should detailed message be printed? Default is TRUE.
           ){ 
    
    ####################################
    # Created Repeated CV Placeholders #
    ####################################
    
    SL_fits<-vector("list", length = reps)
    
    ####################################
    # Repeat Regular IntegratedLearner #
    ####################################
     
    for (r in 1:length(SL_fits)){
      if (verbose) cat('Running K-fold CV iteration', r, "...", "\n")
      SL_fits[[r]]<-run_integrated_learner_CV(feature_table = feature_table,
                                              sample_metadata = sample_metadata,
                                              feature_metadata = feature_metadata,
                                              feature_table_valid = feature_table_valid, 
                                              sample_metadata_valid = sample_metadata_valid, 
                                              folds = folds, 
                                              seed = seed + r, 
                                              base_learner = base_learner, 
                                              meta_learner = meta_learner, 
                                              run_concat = run_concat,  
                                              run_stacked = run_stacked,
                                              verbose = verbose)
    }
    
    
    ##########
    # Return #
    ##########
    
    names(SL_fits)<-paste('Replication', 1:reps, sep ='')
    return(SL_fits)
  
}  


### Function to predict drug lethality/viability on cancer cell-lines choosing an PharmacoGx dataset
### Inputs: ###
# 1. dataset : c("NCI60_2021","GDSC_2020(v2-8.2)","FIMM_2016","UHNBreast_2019","PRISM_2020","gCSI_2019","CTRPv2_2015","GRAY_2017","CCLE_2015")
# 2. model : c('lm','rf','xgTree','neuralnet','gaussprLinear','knn','svmLinear','lasso','ridge','elasticnet')
# 3. no_models : number of models used to infer drug-target interactions (default = 50 models, when we used 50 models)
# 4. data_used : pre=processed data that will be used fro training
# 5. lethality_data_path : the path were every available lethality datase is saved in rds format
# 6. available_targets " the available targets that are going to be used.
viability_model <- function(data_used,
                            available_targets,
                            dataset=c("NCI60_2021","GDSC_2020(v2-8.2)","FIMM_2016",
                                      "UHNBreast_2019","PRISM_2020","gCSI_2019",
                                      "CTRPv2_2015","GRAY_2017","CCLE_2015"),
                            lethality_data_path ='',
                            model = c('lm','rf','xgbTree','neuralnet','gaussprLinear','knn','svmLinear','lasso','ridge','elasticnet'),
                            no_models=50){
  library(tidyverse)
  library(readr)
  library(caret)
  
  # Initialize some further parameters
  if (model=='lasso'){
    lambda_grid <- c(0,10^seq(-5, -1, length = 1000),1)
    alpha_grid <- 1
  }
  if (model=='ridge'){
    lambda_grid <- 10^seq(-3, 1, length = 1000)
    alpha_grid <- 0
  }
  if (model=='elasticnet'){
    lambda_grid <- 10^seq(-3, 1, length = 100)
    alpha_grid <- seq(0.01,0.99,0.05)
  }
  
  ### Load EC50 data
  data <- readRDS(paste0(lethality_data_path,dataset,'.rds'))
  # Make test set based on A549
  # test_smiles_data <- data %>% filter(A549==1)
  test_smiles_data <- data
  available_targets <- available_targets[which(!(available_targets %in% c('drugid','smiles','EC50')))]
  ## make sure that for each drug its targets appear at least once in some other
  targets_analysis <- test_smiles_data[,c('smiles',available_targets)]
  targets_analysis <- targets_analysis %>% gather('target','value',-smiles)
  targets_analysis <- targets_analysis %>% filter(value==1) %>% select(target,smiles) %>% unique()
  targets_analysis <- targets_analysis %>% filter(target %in% available_targets)
  targets_analysis <- targets_analysis %>% group_by(target) %>% mutate(drugs_per_target = n_distinct(smiles)) %>% ungroup()
  targets_analysis <- targets_analysis %>% group_by(smiles) %>% mutate(no_targets_only_there=sum(drugs_per_target==1)) %>% ungroup()
  ### exclude from test sets the drugs with even one of their targets only found with them
  targets_analysis <- targets_analysis %>% filter(no_targets_only_there>0)
  exclude_smiles <- unique(targets_analysis$smiles)
  # repeat to exclude a drug if it appears only in one cell-line
  targets_analysis <- data %>% select(-all_of(c('drugid','EC50',available_targets)))
  targets_analysis <- targets_analysis %>% gather('cell','value',-smiles)
  targets_analysis <- targets_analysis %>% filter(value==1) %>% select(cell,smiles) %>% unique()
  targets_analysis <- targets_analysis %>% group_by(cell) %>% mutate(drugs_per_cell = n_distinct(smiles)) %>% ungroup()
  targets_analysis <- targets_analysis %>% group_by(smiles) %>% mutate(no_cells_only_there=sum(drugs_per_cell==1)) %>% ungroup()
  ### final_exclusion
  targets_analysis <- targets_analysis %>% filter(no_cells_only_there>0)
  exclude_smiles_2 <- unique(targets_analysis$smiles)
  exclude_smiles <- unique(exclude_smiles,exclude_smiles_2)
  test_smiles_data <- test_smiles_data %>% filter(!(smiles %in% exclude_smiles))
  test_smiles_data <- test_smiles_data %>% select(drugid,smiles,all_of(available_targets))
  test_smiles_data <- distinct(test_smiles_data)
  # cv_indices <- createFolds(y = seq_len(nrow(test_smiles_data)), k = 3, list = TRUE, returnTrain = FALSE)
  
  ### LOOCV Training validation procedure
  train_corr_loocv <- NULL
  test_corr_all_loocv <- NULL
  test_preds_loocv <- NULL
  test_preds_shuffled <- NULL
  test_true_loocv <-  NULL
  test_corr_all_shuffled <- NULL
  ypreds_all <- NULL
  ypreds_all_shuffle <- NULL
  ytrue_all <- NULL
  mdls <- NULL
  # k <- 1
  for (loocv_index in 1:nrow(test_smiles_data)){
    loocv_smile <- test_smiles_data$smiles[loocv_index]
    # find test_index
    test_index <- which(data_used$smiles==loocv_smile)
    test_index_a549 <- which(data_used$A549==1 & data_used$smiles==loocv_smile)
    # Subset your data into training and testing sets based on the indices
    train_data_new <- data_used[-test_index, ]
    train_drugs <- train_data_new$smiles
    train_data_new <- train_data_new %>% select(-drugid,-smiles)
    rownames(train_data_new) <- NULL
    test_data_new <- data_used[test_index, ]
    test_drugs <- test_data_new$smiles
    test_data_new <- test_data_new %>% select(-drugid,-smiles)
    rownames(test_data_new) <- NULL
    ### keep separately only LOOCV A549
    test_data_a549 <- data_used[test_index_a549, ]
    test_drugs <- test_data_a549$smiles
    test_data_a549 <- test_data_a549 %>% select(-drugid,-smiles)
    rownames(test_data_a549) <- NULL
    gc()
    
    # build new model
    ctrl <- trainControl(method = "cv", number = 10)
    if (model %in% c('lasso','ridge','elasticnet')){
      mdl <- train(EC50 ~ ., data = train_data_new, method = 'glmnet', trControl = ctrl,tuneGrid = expand.grid(alpha = alpha_grid, lambda = lambda_grid))
    }else{
      invisible(capture.output(mdl <- train(EC50 ~ ., data = train_data_new, method = model, trControl = ctrl)))
    }
    mdls[[loocv_smile]] <- mdl
    y_train_new <- predict(mdl,newdata = train_data_new)
    train_corr_loocv[loocv_smile] <- cor(y_train_new,train_data_new$EC50)
    # first correlation for all test samples
    y_new <- predict(mdl,newdata = test_data_new)
    ypreds_all <- c(ypreds_all,y_new)
    ytrue_all <- c(ytrue_all,test_data_new$EC50)
    test_corr_all_loocv[loocv_smile] <- cor(y_new,test_data_new$EC50)
    # only A549
    y_new <- predict(mdl,newdata = test_data_a549)
    test_preds_loocv <- c(test_preds_loocv,y_new)
    test_true_loocv <- c(test_true_loocv,test_data_a549$EC50)
    
    ### repeat for shuffled labels
    train_data_shuffled <- train_data_new[sample(1:nrow(train_data_new)),]
    # train_data_shuffled[,1:(ncol(train_data_shuffled)-9)] <- tmp_shuffled
    train_data_shuffled$EC50 <- train_data_new$EC50
    if (model %in% c('lasso','ridge','elasticnet')){
      mdl <- train(EC50 ~ ., data = train_data_shuffled, method = 'glmnet', trControl = ctrl,tuneGrid = expand.grid(alpha = alpha_grid, lambda = lambda_grid))
    }else{
      invisible(capture.output(mdl <- train(EC50 ~ ., data = train_data_shuffled, method = model, trControl = ctrl)))
    }
    # first correlation for all test samples
    y_new <- predict(mdl,newdata = test_data_new)
    ypreds_all_shuffle <- c(ypreds_all_shuffle,y_new)
    test_corr_all_shuffled[loocv_smile] <- cor(y_new,test_data_new$EC50)
    # only A549
    y_new <- predict(mdl,newdata = test_data_a549)
    test_preds_shuffled <- c(test_preds_shuffled,y_new)
    
    # print progress
    message(paste0('Finished LOOCV :',loocv_index,' out of ',nrow(test_smiles_data)))
    # message(paste0('Finished fold :',k,' out of ',length(cv_indices)))
    # k <- k+1
  }
  final_results <- data.frame(train=train_corr_loocv,
                                  test = test_corr_all_loocv,
                                  A549_LOOCV = rep(cor(test_preds_loocv,test_true_loocv),length(train_corr_loocv)),
                                  shuffled=test_corr_all_shuffled,
                                  A549_LOOCV_shuffled=rep(cor(test_preds_shuffled,test_true_loocv),length(train_corr_loocv)))
  
  ### LOOCV predictions only for A549
  res_test_new <- data.frame(predicted=test_preds_loocv,EC50 = test_true_loocv)
  res_test_new_shuffled <- data.frame(predicted=test_preds_shuffled,EC50 = test_true_loocv)
  res_test <- rbind(res_test_new %>% mutate(data = 'test'),
                              res_test_new_shuffled %>% mutate(data = 'shuffled'))
  ### LOOCV predictions only for all
  res_test_new_all <- data.frame(predicted=ypreds_all,EC50 = ytrue_all)
  res_test_new_shuffled_all <- data.frame(predicted=ypreds_all_shuffle,EC50 = ytrue_all)
  res_test_all <- rbind(res_test_new_all %>% mutate(data = 'test'),
                              res_test_new_shuffled_all %>% mutate(data = 'shuffled'))
  
  return(list(mdls=mdls,final_results=final_results,res_test=res_test,res_test_all=res_test_all))
}
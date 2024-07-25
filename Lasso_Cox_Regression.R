#Machine learning models on single sample network data
library(tidyr)
library(igraph)
library(glmnet)
library(caret)
library(data.table)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ggsurvfit)
library(pamr)
library(RegParallel)
library(survAUC)
library(survivalROC)
library(risksetROC)
library(c060)
library(rstatix)


#import functions for data wrangling:
setwd("/datos/ot/folder1")
model_format <- function(x){
  rownames(x) <- x[,1]
  x[,1] <- NULL
  x <- t(x) %>% scale()
  return(x)
  
}
cooler_row_merge <- function(x,y){
  merged_df <- merge(x,y,by="row.names")
  rownames(merged_df) <- merged_df[,1]
  merged_df[,1] <- NULL
  return(merged_df)
}

#import weighted degree matrix and clinical data
gts <- fread("lioness_survival_luadtp_gene_targeting_score.tsv", sep = "\t") %>%  as.data.frame()
clinical_data <- readRDS("rnaseq-luad-clinicstg.rds") %>% colData() %>% as.data.frame() %>% filter(., rownames(.) %in% colnames(gts))

#extract covariate info and format
covariates <- clinical_data[c("gender", "age_at_diagnosis")]
covariates$gender <- ifelse(covariates$gender == "female",1,0)

#extract surv info and format
survival <- clinical_data[c("days_to_last_follow_up", "vital_status")]
survival <- survival[survival$days_to_last_follow_up != 0,]
survival$vital_status <- ifelse(survival$vital_status == "Alive", 0, 1) 



#format degree and expr data
gts <- model_format(gts)

#add clinical covariates 
gts <- cooler_row_merge(covariates, gts)
gts <- gts[rownames(survival),]

#make feature matrix (one hot encode categorical values and impute NAs)
gts <- gts[rownames(gts) %in% rownames(survival),] %>% makeX(.,na.impute = T)

#Penalization factors will make sure that clinical covariates are not penalized 
pen_factors_gts <- c(rep(0,2),rep(1,ncol(gts)-length(covariates))) 

#Merge survival object with gts matrix and creaate training test split
data <- cooler_row_merge(gts, survival)

## 70% of the sample size
sample_size <- floor(0.70 * nrow(data))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = sample_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

X_train <- train[,-which(names(train) %in% c("days_to_last_follow_up","vital_status"))] %>% makeX() #%>% scale()
Y_train <- train[,which(names(train) %in% c("days_to_last_follow_up","vital_status"))]
Y_train <- Surv(time=Y_train$days_to_last_follow_up, event = Y_train$vital_status, type = "right")

X_test <- test[,-which(names(test) %in% c("days_to_last_follow_up","vital_status"))] %>% makeX() #%>% scale()
Y_test <- test[,which(names(test) %in% c("days_to_last_follow_up","vital_status"))]
Y_test <- Surv(time=Y_test$days_to_last_follow_up, event = Y_test$vital_status, type = "right")

y_train <- train[,which(names(train) %in% c("days_to_last_follow_up","vital_status"))]
y_test <- test[,which(names(test) %in% c("days_to_last_follow_up","vital_status"))]


#We will now try  testing different LASSO modalities


#Adaptive LASSO Cox Regression Model--------------------------------------------
set.seed(123)
fit <- cv.glmnet(X_train, Y_train, family = "cox", alpha = 0, type.measure = "C",nfolds = 5)
weights <- abs(1 / as.vector(coef(fit, s = "lambda.min")))
weights[c(1,2)] = 0 # don't penalize age and gender

# adaptive Lasso Cox by using cv.glmnet to obtain the 5-fold CV optimal lambda.min or lambda.1se
cvfit <- cv.glmnet(X_train, Y_train, family = "cox", nfolds = 5,type.measure = "C", penalty.factor = weights)
mod <- cvfit$glmnet.fit
lambda_optimal <- cvfit$lambda.min # optimal lambda

predict(mod, s=lambda_optimal, newx = X_test) %>% Cindex(., Y_test)
#-------------------------------------------------------------------------------


#Regular LASSO Model------------------------------------------------------------
fit_lasso <- cv.glmnet(X_train, Y_train, family = "cox", type.measure = "C",alpha = 1, nfolds = 5, penalty.factor = pen_factors_gts)
mod_rl <- fit_lasso$glmnet.fit 
lambda_optimal <- fit_lasso$lambda.1se # optimal lambda
betas <- as.vector(coef(mod_rl, s = lambda_optimal))
beta.positive <- colnames(X_train)[betas > 0]
beta.negative <- colnames(X_train)[betas < 0]
allnames <- names(coef(mod_rl)[, ncol(coef(mod_rl))]
                  [order(coef(mod_rl)[, ncol(coef(mod_rl))], decreasing = TRUE)])
cols <- rep("grey80", length(allnames))
cols[allnames %in% beta.positive] <- "tomato"
cols[allnames %in% beta.negative] <- "darkblue"

plot <- plotmo::plot_glmnet(mod_rl,
                    label = TRUE, s = lambda_optimal, col = cols,
                    xlab = expression(log ~ ~lambda), ylab = expression(beta),
                    main="Standard LASSO Regularization Path"
                    
)

title("Standard LASSO Regularization Path     \n\n")
predict(mod_rl, s=lambda_optimal, newx = X_test) %>% Cindex(., Y_test)



#Relaxed LASSO------------------------------------------------------------------
relaxed_rl <- relax.glmnet(x=X_train, y=Y_train, family = "cox", fit = mod_rl, path = T)
predict(relaxed_rl, s=lambda_optimal,gamma = 0, newx = X_test) %>% Cindex(., Y_test)



#Univariate filtering + penalized cox models------------------------------------
#We perfrom univariate cox regression on the training data only

cox_data <- cooler_row_merge(survival,X_train)
#Univariate filtering:
res <- RegParallel(
  data = cox_data,
  formula = 'Surv(days_to_last_follow_up, vital_status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(cox_data)[5:ncol(cox_data)],
  blocksize = 2000,
  cores = 10,
  nestedParallel = FALSE,
  conflevel = 95)
res <- res[!is.na(res$P),]
res <- res[order(res$LogRank, decreasing = FALSE),]

#Use significance cutoff of 0.2 to avoid false negatives
final <- subset(res, LogRank < 0.05)

#Filter training and test data and keep only selected genes 
train_reselect_univariate <- X_train[,(colnames(X_train) %in% final$Variable) | colnames(X_train) == c("gender","age_at_diagnosis")]
test_reselect_univariate <- X_test[,(colnames(X_test) %in% final$Variable) | colnames(X_test) == c("gender","age_at_diagnosis")]

pen_factors_unifilt_reselect <- c(rep(0,2),rep(1,ncol(train_reselect_univariate)-length(covariates))) 



# univariate filter + Regular lasso---------------------------------------------
cv_unifilt <- cv.glmnet(x=train_reselect_univariate, y=Y_train, family="cox", type.measure = "C", penalty.factor= pen_factors_unifilt_reselect, alpha=1)

mod_url <- cv_unifilt$glmnet.fit 
lambda_optimal <- cv_unifilt$lambda.1se # optimal lambda
betas <- as.vector(coef(mod_url, s = lambda_optimal))
beta.positive <- colnames(train_reselect_univariate)[betas > 0]
beta.negative <- colnames(train_reselect_univariate)[betas < 0]
allnames <- names(coef(mod_url)[, ncol(coef(mod_url))]
                  [order(coef(mod_url)[, ncol(coef(mod_url))], decreasing = TRUE)])
cols <- rep("gray80", length(allnames))
cols[allnames %in% beta.positive] <- "seagreen3"
cols[allnames %in% beta.negative] <- "hotpink"




# png("Unifilt_LASSO_regpath.png", res=216, width = 1440, height = 1440)
# plotmo::plot_glmnet(mod_url,
#                     label = TRUE, s = lambda_optimal, col = cols,
#                     xlab = expression(log ~ ~lambda), ylab = expression(beta),
#                     main="Univariate Filter + LASSO Regularization Path"
# )
# 
# dev.off()



predict(mod_url, s=lambda_optimal, newx = test_reselect_univariate) %>% Cindex(., Y_test)

# png("Unifilt_LASSO_CV_plot.png", res=216, width = 1440, height = 1440)
# plot(cv_unifilt, main="10-fold Cross Validation: Univariate Filter + LASSO  \n\n")
# dev.off()
# 

#-------------------------------------------------------------------------------

#univariate filter + Relaxed lasso----------------------------------------------
mod_url_relax <- relax.glmnet(x=train_reselect_univariate, y=Y_train, family="cox", fit =mod_url ,path = T)
predict(mod_url_relax, s=lambda_optimal, gamma = 0, newx = test_reselect_univariate) %>% Cindex(., Y_test)
plot(mod_url_relax)

cplot_test <- coefplot::coefplot(mod_url_relax, sort="magnitude", lambda = lambda_optimal)

# Get gene symbols
gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                      filters = 'ensembl_gene_id',
                      values = df$ensembl_id,
                      mart = ensembl)

# Merge with original data frame
df <- merge(df, gene_symbols, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE)

print(df)

#Obtain scoefficients for vizualization 
coefs_u.l.rlx <- get_coefs(model = mod_url_relax,  lambda = lambda_optimal) %>% as.matrix() %>% as.data.frame()
colnames(coefs_u.l.rlx) <- "Coefficient"
coefs_u.l.rlx <- coefs_u.l.rlx %>% .[order(-.$Coefficient),,drop=F]

fwrite(coefs_u.l.rlx, "coefficients_u_l_rlx.tsv", sep = "\t", row.names = T)


#Relaxing the fit shows better performance
--------------------------------------------------------------------------------

  
  
  
#Model validation***************************************************************
  
#Measure goodness of fit by dichotomizing predicted risk scores:----------------

# univariate filter + Relaxed Lasso:
pred_lp <- predict(mod_url_relax, s=lambda_optimal, gamma = 0, newx = test_reselect_univariate, type = "link")

# dichotomize by prognostic scores (linear predictor)  by median to divide the validation patients into two groups
group_dichotomize <- as.numeric(pred_lp > median(pred_lp))

# draw two survival curves based on KM estimation and compare them by a log-rank test
dat_tmp <- data.frame(time = y_test[, 1], status = y_test[, 2], group = group_dichotomize)
sfit <- survfit(Surv(time, status) ~ group, data = dat_tmp)

ggsurv <- ggsurvplot(sfit,
                     conf.int = TRUE, risk.table = TRUE,
                     xlab = "Days to last follow up", legend = c(.2, .3),
                     legend.labs = c("Low risk", "High risk"), legend.title = "Dichotomized groups",
                     risk.table.y.text.col = TRUE, risk.table.y.text = FALSE
)
ggsurv$plot <- ggsurv$plot +
  annotate("text", x = 500.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))
ggsurv$table <- ggsurv$table + labs(y = "Dichotomized\n groups")
ggsurv


# png("Relaxed_LASSO_riskgroups_KM_TEST_plot.png", res = 216, width =1200 ,height = 1440)
# ggsurv
# dev.off()


# #KM on stratified linear predictors of TRAINING data:
# # univariate filter + Relaxed Lasso:
# pred_lp <- predict(mod_url_relax, s=lambda_optimal, gamma = 0, newx = train_reselect_univariate, type = "link")
# 
# # dichotomize by prognostic scores (linear predictor)  by median to divide the validation patients into two groups
# group_dichotomize <- as.numeric(pred_lp > median(pred_lp))
# 
# # draw two survival curves based on KM estimation and compare them by a log-rank test
# dat_tmp <- data.frame(time = y_train[, 1], status = y_train[, 2], group = group_dichotomize)
# sfit <- survfit(Surv(time, status) ~ group, data = dat_tmp)
# 
# ggsurv <- ggsurvplot(sfit,
#                      conf.int = TRUE, risk.table = TRUE,
#                      xlab = "Days to last follow up", legend = c(.2, .3),
#                      legend.labs = c("Low risk", "High risk"), legend.title = "Dichotomized groups",
#                      risk.table.y.text.col = TRUE, risk.table.y.text = FALSE
# )
# ggsurv$plot <- ggsurv$plot +
#   annotate("text", x = 500.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))
# ggsurv$table <- ggsurv$table + labs(y = "Dichotomized\n groups")
# ggsurv
# 
# 
# png("Relaxed_LASSO_riskgroups_KM_TRAINING_plot.png", res = 216, width =1200 ,height = 1440)
# ggsurv
# dev.off()


#univariate filter + regular Lasso:
pred_lp_rl <- predict(mod_url, s=lambda_optimal, newx = test_reselect_univariate, type = "link")

# dichotomize by prognostic scores (linear predictor)  by median to divide the validation patients into two groups
group_dichotomize <- as.numeric(pred_lp_rl > median(pred_lp_rl))

# draw two survival curves based on KM estimation and compare them by a log-rank test
dat_tmp <- data.frame(time = y_test[, 1], status = y_test[, 2], group = group_dichotomize)
sfit <- survfit(Surv(time, status) ~ group, data = dat_tmp)

ggsurv <- ggsurvplot(sfit,
                     conf.int = TRUE, risk.table = TRUE,
                     xlab = "Time since diagnosis (year)", legend = c(.2, .3),
                     legend.labs = c("Low risk", "High risk"), legend.title = "Dichotomized groups",
                     risk.table.y.text.col = TRUE, risk.table.y.text = FALSE
)
ggsurv$plot <- ggsurv$plot +
  annotate("text", x = 2.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))
ggsurv$table <- ggsurv$table + labs(y = "Dichotomized\n groups")
ggsurv

                            
#Relaxed fit is able to split the two risk groups significantly. Regular lasso 
#is not able to achieve significance (p < 0.05) but displays the same overall trend

#-------------------------------------------------------------------------------

#Measure discriminatory power locally (discrete times) or globally (across all time)

#survival ROC curves

png("Relaxed_LASSO_AUC.png", width = 480*4,height = 600, res = 72*4)
## 2 x 5 layout
layout(matrix(1:4, byrow = T, ncol = 4))

#univariate + relaxed lasso
ROC <- risksetROC(
  Stime = y_test[, 1], status = y_test[, 2],
  marker = pred_lp, predict.time = 365.25, method = "Cox",
  main = "1 year survival", col = "seagreen3", type = "s",
  lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity"
)
text(0.7, 0.2, paste("AUC =", round(ROC$AUC, 3)))

ROC_2 <- risksetROC(
  Stime = y_test[, 1], status = y_test[, 2],
  marker = pred_lp, predict.time = 730.5, method = "Cox",
  main = "2 year survival", col = "seagreen3", type = "s",
  lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity"
)
text(0.7, 0.2, paste("AUC =", round(ROC_2$AUC, 3)))

ROC_3 <- risksetROC(
  Stime = y_test[, 1], status = y_test[, 2],
  marker = pred_lp, predict.time = 1095.75, method = "Cox",
  main = "3 year survival", col = "seagreen3", type = "s",
  lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity"
)
text(0.7, 0.2, paste("AUC =", round(ROC_3$AUC, 3)))

ROC_5 <- risksetROC(
  Stime = y_test[, 1], status = y_test[, 2],
  marker = pred_lp, predict.time = 1826.25, method = "Cox",
  main = "5 year survival", col = "seagreen3", type = "s",
  lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity"
)
text(0.7, 0.2, paste("AUC =", round(ROC_5$AUC, 3)))

dev.off()

  
pred_lp_train <- predict(mod_url_relax, s=lambda_optimal, gamma = 0, newx = train_reselect_univariate, type = "link")

# Time dependent Brier Score:

data_validate_regular <- data.frame(time = y_test[, "days_to_last_follow_up"], status = y_test[, "vital_status"], lp = as.vector(pred_lp_rl))
data_validate_rlx <- data.frame(time = y_test[, "days_to_last_follow_up"], status = y_test[, "vital_status"], lp = as.vector(pred_lp))



lasso_validate_regular <- coxph(Surv(time, status) ~ lp, data = data_validate_regular, y = TRUE, x = TRUE)
lasso_validate_rlx <- coxph(Surv(time, status) ~ lp, data = data_validate_rlx, y = TRUE, x = TRUE)

Brier_validate_regular <- riskRegression::Score(list("Regular_Lasso" = lasso_validate_regular), formula = Surv(time, status) ~ 1, data = data_validate_regular, conf.int = FALSE, metrics = "brier", summary = "ibs", times = sort(unique(data_validate_regular$time)))$Brier$score
Brier_validate_rlx <- riskRegression::Score(list("UF_Relaxed_Lasso" = lasso_validate_rlx), formula = Surv(time, status) ~ 1, data = data_validate_rlx, conf.int = FALSE, metrics = "brier", summary = "ibs", times = sort(unique(data_validate_rlx$time)))$Brier$score



Brier_score <- rbind(Brier_validate_regular,  Brier_validate_rlx)
Brier_score <- Brier_score[Brier_score$model != "Null model", ]


ggplot(Brier_score, aes(times, Brier, group = model, color = model)) +
  xlab("Evaluation time points (year)") +
  ylab("Brier score") +
  geom_step(direction = "vh") +
  theme(legend.position = c(0.15, 0.88), legend.title = element_blank())


#Integrated Brier Score:
Brier_validate_ibs <- Brier_validate[Brier_validate$model == "Brier_validate", ]
Brier_validate_ibs$IBS[which.max(Brier_validate_ibs$times)]

Brier_validate_regular_ibs <- Brier_validate_regular[Brier_validate_regular$model == "Brier_validate_regular_lasso", ]
Brier_validate_regular_ibs$IBS[which.max(Brier_validate_regular_ibs$times)]

#-------------------------------------------------------------------------------





#Measure uncertainty in models' performance by resampling-based methods---------

#define function to extract active coeficients (i.e. coefs != 0)
get_coefs <- function(model, lambda){
coefs_active <- coef.glmnet(model, s=lambda)
coefs_sign_active <- coefs_active[coefs_active[,1] != 0,,drop=F] #get active features != 0
return(coefs_sign_active)
}

#define function for univariate filtering
univariate_filter <- function(survival,X_train,p.val){
  cox_data <- cooler_row_merge(survival,X_train)
  #Univariate filtering:
  res <- RegParallel(
    data = cox_data,
    formula = 'Surv(days_to_last_follow_up, vital_status) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(cox_data)[5:ncol(cox_data)],
    blocksize = 2000,
    cores = 10,
    nestedParallel = FALSE,
    conflevel = 95)
  res <- res[!is.na(res$P),]
  res <- res[order(res$LogRank, decreasing = FALSE),]
  
  #Use significance cutoff of 0.2 to avoid false negatives
  final <- subset(res, LogRank < p.val)
  return(final)
  
  
}

#define function to calculate integrated Brier's score  

ibs <- function(y_test,pred_lp){
  #Integrated Brier's score:
  
  # prepare data format suited for function Score() from the riskRegression package
  data_validate <- data.frame(time = y_test[, "days_to_last_follow_up"], status = y_test[, "vital_status"], lp = as.vector(pred_lp))
  lasso_validate <- coxph(Surv(time, status) ~ lp, data = data_validate, y = TRUE, x = TRUE)
  
  # calculate Brier scores based on both training and validation data
  Brier_validate <- riskRegression::Score(list("Brier_validate" = lasso_validate), formula = Surv(time, status) ~ 1, data = data_validate, conf.int = FALSE, metrics = "brier", summary = "ibs", times = sort(unique(data_validate$time)))$Brier$score
  
  #Calculate Integrated Brier Score:
  Brier_validate_ibs <- Brier_validate[Brier_validate$model == "Brier_validate", ]
  i_score <- Brier_validate_ibs$IBS[which.max(Brier_validate_ibs$times)]
return(i_score)
}


#define function for manual model relaxing
# manual_relax <- function(train_reselect_univariate,Y_train){
#   
#   evaluate_cv_perf <- cv.glmnet(x=train_reselect_univariate, y=Y_train, family="cox", type.measure = "C", penalty.factor= pen_factors_unifilt_reselect, alpha=1)
#   investigate_evaluate <- coef(evaluate_cv_perf, s=evaluate_cv_perf$lambda.1se)
#   number_of_var_evaluate <- investigate_evaluate[investigate_evaluate[,1] != 0,] %>% length()
#   
#   unpenal_factor_evaluate <- rep(0, number_of_var_evaluate)
#   vars_evaluate <-  investigate_evaluate[investigate_evaluate[,1] != 0,] %>% as.data.frame()
#   
#   x_train_filt_lasso_univar_evaluate <- X_train[,(colnames(X_train) %in% rownames(vars_evaluate))]
#   x_test_filt_lasso_univar_evaluate <- X_test[,(colnames(X_test) %in% rownames(vars_evaluate))]
#   
#   fit_unpenal_evaluate <- glmnet(x=x_train_filt_lasso_univar_evaluate, y = Y_train, family="cox", alpha=1, nlambda = 1000)
#   
# 
#   
#   selected_genes_model.u.l <-  c(selected_genes_model.u.l, list(rownames(coefs.u.l)))
#   
#   
#   
#   pr <- predict(fit_unpenal_evaluate, s=0, newx =  x_test_filt_lasso_univar_evaluate) %>% Cindex(., Y_test)
#   
#   return(pr)
#   
#   
# }





#merge data frame with target variables and weighted degrees
#Merge survival object with gts matrix and creaate training test split
data <- cooler_row_merge(gts, survival)


#We randomly split the data 100 times and perform univariate filtering, lasso fitting and discrimination measures in each iteration.





k = 100
summarized_performance <- cbind.data.frame(C.index.l=NA,C.index.u.l=NA, C.index.u.l.rlx=NA,
                                           Uno.l=NA,Uno.u.l=NA,Uno.u.l.rlx=NA, 
                                           Brier.l=NA,Brier.u.l=NA,Brier.u.l.rlx=NA,
                                           iteration=NA)
selected_genes_model.l <- c()
selected_genes_model.u.l <- c()

set.seed(123)

for (i in 1:k) {
  
#random split of training and test data: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## 70% of the sample size
  sample_size <- floor(0.70 * nrow(data))

  ## set the seed to make your partition reproducible
  train_ind <- sample(seq_len(nrow(data)), size = sample_size)

  train <- data[train_ind, ]
  test <- data[-train_ind, ]
  
  X_train <- train[,-which(names(train) %in% c("days_to_last_follow_up","vital_status"))] %>% makeX() 
  Y_train <- train[,which(names(train) %in% c("days_to_last_follow_up","vital_status"))]
  Y_train <- Surv(time=Y_train$days_to_last_follow_up, event = Y_train$vital_status, type = "right")
  
  X_test <- test[,-which(names(test) %in% c("days_to_last_follow_up","vital_status"))] %>% makeX()
  Y_test <- test[,which(names(test) %in% c("days_to_last_follow_up","vital_status"))]
  Y_test <- Surv(time=Y_test$days_to_last_follow_up, event = Y_test$vital_status, type = "right")
  
  y_train <- train[,which(names(train) %in% c("days_to_last_follow_up","vital_status"))]
  y_test <- test[,which(names(test) %in% c("days_to_last_follow_up","vital_status"))]
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  
#Lasso on all genes ------------------------------------------------------------
cv_lasso <- cv.glmnet(X_train, Y_train, family = "cox", type.measure = "C",alpha = 1, nfolds = 5, penalty.factor = pen_factors_gts)
lambda_optimal.l <- cv_lasso$lambda.1se # optimal lambda
mod.l <- cv_lasso$glmnet.fit 
coefs.l <- get_coefs(model=mod.l, lambda=lambda_optimal.l)

selected_genes_model.l <-  c(selected_genes_model.l, list(rownames(coefs.l)))

predict.l <- predict(mod.l, s=lambda_optimal.l, newx = X_test, type = "link")

#-------------------------------------------------------------------------------



#Univariate filtering-----------------------------------------------------------

#apply uf to training set only
uf <- univariate_filter(survival=y_train,X_train = X_train, p.val = 0.05)

#Filter training and test data and keep only selected genes 
train_reselect_univariate <- X_train[,(colnames(X_train) %in% uf$Variable) | colnames(X_train) == c("gender","age_at_diagnosis")]
test_reselect_univariate <- X_test[,(colnames(X_test) %in% uf$Variable) | colnames(X_test) == c("gender","age_at_diagnosis")]

#get penalization factor for filtered dat
pen_factors_unifilt_reselect <- c(rep(0,2),rep(1,ncol(train_reselect_univariate)-length(covariates))) 

#-------------------------------------------------------------------------------



#lasso on uf data---------------------------------------------------------------
cv_unifilt <- cv.glmnet(x=train_reselect_univariate, y=Y_train, family="cox", type.measure = "C", penalty.factor= pen_factors_unifilt_reselect, alpha=1)
mod.u.l <- cv_unifilt$glmnet.fit 
lambda_optimal.u.l <- cv_unifilt$lambda.1se # optimal lambda

coefs.u.l <- get_coefs(model=mod.u.l, lambda = lambda_optimal.u.l)

selected_genes_model.u.l <-  c(selected_genes_model.u.l, list(rownames(coefs.u.l)))

predict.u.l <- predict(mod.u.l, s=lambda_optimal.u.l, newx = test_reselect_univariate, type = "link")


#-------------------------------------------------------------------------------

#relaxed lasso on uf data-------------------------------------------------------

mod.u.l.rlx <- relax.glmnet(x=train_reselect_univariate, y=Y_train, family="cox", fit =mod.u.l ,path = T)

predict.u.l.rlx <- predict(mod.u.l.rlx, s=lambda_optimal.u.l, gamma=0, newx = test_reselect_univariate, type = "link")


#-------------------------------------------------------------------------------


#Harrel's C-Index
Cindex.l <- Cindex(predict.l, y=Y_test)
Cindex.u.l <- Cindex(predict.u.l, y=Y_test)
Cindex.u.l.rlx <- Cindex(predict.u.l.rlx, y=Y_test)

#Uno's C-Index
Uno.l <- UnoC(Surv.rsp = Y_train, Surv.rsp.new = Y_test, lpnew = predict.l)
Uno.u.l <- UnoC(Surv.rsp = Y_train, Surv.rsp.new = Y_test, lpnew = predict.u.l)
Uno.u.l.rlx <- UnoC(Surv.rsp = Y_train, Surv.rsp.new = Y_test, lpnew = predict.u.l.rlx)

#Integrated Brier's Score:

Ibs.l <- ibs(y_test = y_test,pred_lp = predict.l)
Ibs.u.l <- ibs(y_test = y_test, pred_lp = predict.u.l)
Ibs.u.l.rlx <- ibs(y_test = y_test, pred_lp = predict.u.l.rlx)

  
  #fill performance metrics and capture selected genes in models
  
summarized_performance[i,1] <-Cindex.l
summarized_performance[i,2]<- Cindex.u.l
summarized_performance[i,3]<- Cindex.u.l.rlx
summarized_performance[i,4]<- Uno.l
summarized_performance[i,5]<- Uno.u.l
summarized_performance[i,6]<- Uno.u.l.rlx
summarized_performance[i,7]<- Ibs.l
summarized_performance[i,8]<- Ibs.u.l
summarized_performance[i,9]<- Ibs.u.l.rlx
summarized_performance[i,10]<- i


  
  print(paste0("Iteration: ",i)) 
  
}




summarized_performance



selected_genes_l <- unlist(selected_genes_model.l) %>% table() %>% as.data.frame()
selected_genes_l <- selected_genes_l[order(selected_genes_l$Freq, decreasing = T),]


selected_genes_u.l <- unlist(selected_genes_model.u.l) %>% table() %>% as.data.frame()
selected_genes_u.l <- selected_genes_u.l[order(selected_genes_u.l$Freq, decreasing = T),]


gathered.Harrel <- summarized_performance %>% gather(key = "variable", value = "value", c("C.index.l", "C.index.u.l","C.index.u.l.rlx"))
gathered.Harrel <- gathered.Harrel %>% select(., c(variable, value))
colnames(gathered.Harrel) <- c("Model","C-Index")
gathered.Harrel <- gathered.Harrel %>%
                   mutate(Model = recode(Model, 
                        "C.index.u.l.rlx" = "U+L+R", 
                        "C.index.u.l" = "U+L", 
                        "C.index.l" = "L"))



resample_results_plot <- ggplot(gathered.Harrel, aes(x = variable, y = value, col = variable)) +
  geom_boxplot()+ geom_jitter(size=1,alpha=0.5)+
  geom_hline(yintercept = 0.5, col = "red")+
  ylim(.3,1)+
  labs(x = "Model", y = "C-Index", title = "") 

stat.test <- gathered.Harrel %>%
  wilcox_test(`C-Index` ~ Model, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Model")
bxp_cindex <- ggboxplot(gathered.Harrel, x = "Model", y = "`C-Index`", fill = "Model", 
                 palette = c("#00AFBB", "#E7B800", "#FC4E07")) + geom_jitter(col="black", alpha=0.2)+ geom_hline(yintercept = .50, col="red", linetype="dashed")+
  stat_pvalue_manual(stat.test, label = paste0("{p.adj.signif} ","{p.adj}"), tip.length = 0.01,step.increase = c(.1,.1,.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
bxp_cindex



ggsave(
  filename = "Resample_bxplt.png",
  plot = bxp_cindex,
  width = 8,
  height = 8,dpi = 150
)



gathered.ibs <- summarized_performance %>% gather(key = "variable", value = "value", c("Brier.l", "Brier.u.l","Brier.u.l.rlx"))

gathered.ibs <- gathered.ibs %>% select(., c(variable, value))
colnames(gathered.ibs) <- c("Model","Integrated Brier Score")
gathered.ibs <- gathered.ibs %>%
  mutate(Model = recode(Model, 
                        "Brier.u.l.rlx" = "U+L+R", 
                        "Brier.u.l" = "U+L", 
                        "Brier.l" = "L"))


stat.test <- gathered.ibs %>%
  wilcox_test(`Integrated Brier Score` ~ Model, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Model")
bxp_ibs <- ggboxplot(gathered.ibs, x = "Model", y = "Integrated Brier Score", fill = "Model", 
                        palette = c("#00AFBB", "#E7B800", "#FC4E07")) + geom_jitter(col="black", alpha=0.2)+ geom_hline(yintercept = .25, col="red", linetype="dashed")+
  stat_pvalue_manual(stat.test, label = paste0("{p.adj.signif} "), tip.length = 0.01,step.increase = c(.1,.1,.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
bxp_ibs


bxps_grid_h <- ggarrange(bxp_cindex, bxp_ibs, ncol = 2, nrow = 1, common.legend = T)


ggsave(
  filename = "performance_resample_cindex_ibs_vertical.png",
  plot = bxps_grid_h,
  width = 4,
  height = 8,dpi = 150
)





ggplot(gathered.ibs, aes(x = variable, y = value, col = variable)) +
  geom_boxplot() + geom_jitter(size=1,alpha=0.5)+
  geom_hline(yintercept = 0.25, col="red")+
  ylim(0,.30)+
  labs(x = "Model", y = "integrated Brier's Score", title = "") + theme_bw()

gathered.Uno <- summarized_performance %>% gather(key = "variable", value = "value", c("Uno.l", "Uno.u.l","Uno.u.l.rlx"))

ggplot(gathered.Uno, aes(x = variable, y = value, col = variable)) +
  geom_boxplot() + geom_jitter(size=0.5,alpha=0.5)+
  ylim(0,1)+
  labs(x = "Model", y = "Uno C Statistic", title = "")

  
selected_genes_l
selected_genes_u.l

fwrite(selected_genes_u.l, "genes_frequencies_resampling.tsv", sep="\t")
  


#Vizualize 2 genes

#CHRDL2---------------------------------------
viz_data <- cooler_row_merge(survival, gts)


high <- quantile(viz_data$ENSG00000054938, .75)
low <- quantile(viz_data$ENSG00000054938, .25)


high <- median(viz_data$ENSG00000054938)


viz_data$class <- ifelse(viz_data$ENSG00000054938 >= high, "High CHRDL2 Degree",
                         ifelse(viz_data$ENSG00000054938 <= low, 'Low CHRDL2 Degree', 'Mid CHRDL2 Degree')) %>% as.factor()

viz_data$class <- ifelse(viz_data$ENSG00000054938 >= high, "High CHRDL2 Weighted Degree",'Low CHRDL2 Weighted Degree') %>% as.factor()

surv.plot <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ class, data = viz_data) %>%
  ggsurvplot(conf.int = TRUE, risk.table = TRUE,
             xlab = "Days to last follow up", legend = c(.2, .3),
             legend.labs = c("Low CHRDL2 Weighted Degree"," High CHRDL2 Weighted Degree"), legend.title = "Dichotomized groups",
             risk.table.y.text.col = TRUE, risk.table.y.text = FALSE)
surv.plot

fit_chrdl <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ class, data = viz_data)  
 ggsurv_chr <- ggsurvplot(fit_chrdl,conf.int = T,
                          xlab = "Days to last follow up", legend = c(.2, .3),
                          legend.labs = c("Low CHRDL2 Weighted Degree"," High CHRDL2 Weighted Degree"), legend.title = "Dichotomized groups", palette = c("tomato","Blue"))

ggsurv_chr$plot <- ggsurv_chr$plot +
  annotate("text", x = 1000.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(fit_chrdl)$pval.txt))

CHRDL2 <- ggsurv_chr
CHRDL2
#----------------------------------------------

#SPP2------------------------------------------
viz_data <- cooler_row_merge(survival, gts)


high <- quantile(viz_data$ENSG0000007208, .75)
low <- quantile(viz_data$ENSG0000007208, .25)


high <- median(viz_data$ENSG0000007208)


viz_data$class <- ifelse(viz_data$ENSG00000054938 >= high, "High CHRDL2 Degree",
                         ifelse(viz_data$ENSG00000054938 <= low, 'Low CHRDL2 Degree', 'Mid CHRDL2 Degree')) %>% as.factor()

viz_data$class <- ifelse(viz_data$ENSG0000007208 >= high, "High SPP2 Weighted Degree",'Low SPP2 Weighted Degree') %>% as.factor()


fit_spp2 <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ class, data = viz_data)  
ggsurv_spp2 <- ggsurvplot(fit_spp2,conf.int = T,
                         xlab = "Days to last follow up", legend = c(.2, .3),
                         legend.labs = c("Low SPP2 Weighted Degree"," High SPP2 Weighted Degree"), legend.title = "Dichotomized groups", palette = c("tomato","Blue"))

ggsurv_spp2$plot <- ggsurv_spp2$plot +
  annotate("text", x = 1000.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(fit_spp2)$pval.txt))

SPP2 <- ggsurv_spp2
SPP2

png("CHRLD2_KM.png", res = 216, width =1200 ,height = 1200)
CHRDL2
dev.off()

png("SPP2_KM.png", res = 216, width =1200 ,height = 1200)
SPP2
dev.off()

#generalize----------------------------------------------

top_stable <- head(selected_genes_u.l, 14)[-(1:2),][,1] %>% as.character()


plots_list <- list()

for (i in 1:length(top_stable)) {

viz_data <- cooler_row_merge(survival, gts)

high <- as.matrix(viz_data[top_stable[i]]) %>% as.numeric %>% median()


viz_data$class <- ifelse((viz_data[top_stable[i]] %>% as.matrix() %>% as.numeric) >= high, "High Weighted Degree",'Low Weighted Degree') %>% as.factor()


fit_general <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ class, data = viz_data)  

ggsurv_general <- ggsurvplot(fit_general,conf.int = T,
                          xlab = "Days to last follow up", legend = c(.2, .3),
                          legend.labs = c("Low Weighted Degree"," High Weighted Degree"), legend.title = "Dichotomized groups", palette = c("tomato","Blue"))

ggsurv_general$plot <- ggsurv_general$plot +
  annotate("text", x = 2000.6, y = .03, label = paste0("Log-rank test:\n", surv_pvalue(fit_general)$pval.txt))

general <- ggsurv_general$plot

plots_list[[i]] <- ggsurv_general


}

grid_km <- arrange_ggsurvplots(plots_list, ncol = 4, nrow = 3, title = "Weighted Gene Degrees", print = F)

ggsave(plot = grid_km,width = 15, height = 15, dpi = 150,filename = "grid_km.png")


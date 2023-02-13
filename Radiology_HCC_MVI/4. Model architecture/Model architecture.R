#### ICC #####

library(psych)

T1 <- read.csv("radiomics features/zhejiang all radomics features.csv",header = T)
T2 <- read.csv("radiomics features/zhejiang all radomics features repeat.csv",header = T)
name <- c(T2[,1]);T1_2 <- subset(T1,ID %in% name[1:50]) 
T1_3 <- T1_2;T2_2 <- T2
T1_3$ID <- NULL;T2_2$ID <- NULL
x <- dim(T1_3)[1]
y <- dim(T2_2)[2]
T12 <- cbind(T1_3,T2_2);dim(T12)

icc <- c(1:y)
for(i in 1:y) {icc[i] <- ICC(T12[,c(i,i+y)])$results[2,2]} #inter-class [2,2]; intra-class [1,2]
mean(icc);median(icc)
h <- length(which(icc > 0.75));h
l <- length(which(icc <= 0.75));l

A <- which(icc > 0.75)
A <- A+1
b <- c(1)      
col <- c(b,A)
data1 <- T1[,col]
write.csv(data1, row.names = FALSE,file="radiomics features/zheyi all radomics features ICC.csv")

#### z-score ####

library(caret)

center1 <- read.csv("radiomics features/zheyi all radomics features ICC",header = T)
center23 <- read.csv("radiomics features/Nantongxiangya all radomics features ICC.csv",header = T)
center4 <- read.csv("radiomics features/yunnan all radomics features ICC.csv",header = T)

normal_para1 <- preProcess(x = train,method = c("center","scale"))
normal_para2 <- preProcess(x = test1,method = c("center","scale"))
normal_para3 <- preProcess(x = test2,method = c("center","scale"))

df_center1 <- predict(object = normal_para1,newdata = center1)
df_center23 <- predict(object = normal_para2,newdata = center23)
df_center4 <- predict(object = normal_para3,newdata = center4)

write.csv(df_center1, row.names = FALSE,file="radiomics features/zheyi all radomics features z.csv")
write.csv(df_center23, row.names = FALSE,file="radiomics features/nantongxiangya all radomics features z.csv")
write.csv(df_center4, row.names = FALSE,file="radiomics features/yunnan all radomics features z.csv")

#### model building ####

library(caret)
library(glmnet)
library(rms)
library(dplyr)
library(mRMRe) 

rm(list=ls())

#It is worth noting that Valida* refers to an internal test set and test* refers to an external test set. 
#This naming is done for easy coding and less error.
#Importance!!! The radiomics features selected process in five Radiomics models were same as following with different input radiomics features (types). 

#### Zhejiang center divided into the training set and internal test set ####

df <- read.csv("radiomics features/zheyi all radomics features z.csv")
df2 <- read.csv("newdata2/radiomics features/nantongxiangya all radomics features z.csv")

set.seed(100)
index <- createDataPartition(df$Label,p = 0.7)
training_set <- df[index$Resample1,]
validation_set <- df[-index$Resample1,]
# write.csv(training_set,file = "training.csv",row.names = FALSE)
# write.csv(validation_set,file = "test1.csv",row.names = FALSE)
prop.table(table(training_set$Label))
prop.table(table(validation_set$Label))

train_label <- training_set[, 2]
valid_label <- validation_set[, 2]
test_label <- df2[, 2]
train_set_nolabel <- training_set[,-c(1:2)]
dim(training_set)
dim(validation_set)

####  mRMR ####

mrmr_feature<-train_set_nolabel 
mrmr_feature$y<-train_label
target_indices = which(names(mrmr_feature)=='y')

for (m in which(sapply(mrmr_feature, class)!="numeric")){
  mrmr_feature[,m]=as.numeric(unlist(mrmr_feature[,m]))
}

data4 <- mRMR.data(data = data.frame(mrmr_feature))
mrmr=mRMR.ensemble(data = data4, target_indices = target_indices,
                   feature_count = 20, solution_count = 1)
index=mrmr@filters[[as.character(mrmr@target_indices)]]
index=as.numeric(index)
data_reduce = train_set_nolabel[,index]

#### LASSO ####

cv_x <- as.matrix(data_reduce)
cv_y <- train_label

set.seed(1000)
lasso_selection <- cv.glmnet(x=cv_x,
                             y=cv_y,
                             family = "binomial",
                             type.measure = "deviance",
                             alpha = 1,
                             nfolds = 5)

par(font.lab = 2, mfrow = c(2,1), mar = c(4.5,5,3,2))
plot(x = lasso_selection, las = 1, xlab = "Log(lambda)")
nocv_lasso <- glmnet(x = cv_x, y = cv_y, family = "binomial",alpha = 1)
plot(nocv_lasso,xvar = "lambda",las=1,lwd=2,xlab="Log(lambda)")
abline(v = log(lasso_selection$lambda.min),lwd=1,lty=3,col="black")

coefPara <- coef(object = lasso_selection,s="lambda.min")
lasso_values <- as.data.frame(which(coefPara != 0, arr.ind = T))

lasso_names <- rownames(lasso_values)[-1]
lasso_coef <- data.frame(Feature = rownames(lasso_values),
                         Coef = coefPara[which(coefPara !=0,arr.ind = T)])
lasso_coef
lasso_coef_len <- length(lasso_coef[,1])
lasso_coef_save <- lasso_coef[,1][2:lasso_coef_len]
lasso_coef_save

train_set_lasso <- data.frame(cv_x)[lasso_names]
valid_set_lasso <- validation_set[names(train_set_lasso)]
test_set_lasso <- df2[names((train_set_lasso))]
test_all = as.matrix(test_set_lasso)
Data_all = as.matrix(rbind(train_set_lasso,valid_set_lasso))
xn = nrow(Data_all)
yn = ncol(Data_all)
xn2 = nrow(test_set_lasso)
yn2 = ncol(test_set_lasso)

beta = as.matrix(coefPara[which(coefPara !=0),])
betai_Matrix = as.matrix(beta[-1])
beta0_Matrix = matrix(beta[1],xn,1)
Radcore_Matrix = Data_all %*% betai_Matrix +beta0_Matrix
radscore_all = as.numeric(Radcore_Matrix)

beta0_Matrix2 = matrix(beta[1],xn2,1)
Radcore_Matrix2 = test_all %*% betai_Matrix +beta0_Matrix2
Radscore_test = as.numeric(Radcore_Matrix2)

Radscore_train = radscore_all[1:nrow(train_set_lasso)]
Radscore_valid = radscore_all[(nrow(train_set_lasso)+1):xn]

Radscore_train_matrix <- matrix(Radscore_train,ncol = 1)
predata_train_matrix <- data.frame(training_set[,c(1:2)])
lasso_coef_train_matrix <- select(training_set,lasso_coef_save)
radscore_train_data1 <- cbind(predata_train_matrix,Radscore_train_matrix)
radscore_train_data2 <- cbind(predata_train_matrix,lasso_coef_train_matrix)
write.csv(radscore_train_data2,file = "select feature/train radiomics features.csv",row.names = FALSE)

Radscore_valid_matrix <- matrix(Radscore_valid,ncol = 1)
predata_valid_matrix <- data.frame(validation_set[,c(1:2)])
lasso_coef_valid_matrix <- select(validation_set,lasso_coef_save)
radscore_valid_data1 <- cbind(predata_valid_matrix,Radscore_valid_matrix)
radscore_valid_data2 <- cbind(predata_valid_matrix,lasso_coef_valid_matrix)
write.csv(radscore_valid_data2,file = "select feature3/test1 radiomics features.csv",row.names = FALSE)

Radscore_test_matrix <- matrix(Radscore_test,ncol = 1)
predata_test_matrix <- data.frame(df2[,c(1:2)])
lasso_coef_test_matrix <- select(df2,lasso_coef_save)
radscore_test_data1 <- cbind(predata_test_matrix,Radscore_test_matrix)
radscore_test_data2 <- cbind(predata_test_matrix,lasso_coef_test_matrix)
write.csv(radscore_test_data2,file = "select feature3/test2 radiomics features.csv",row.names = FALSE)

####  Metrics ####

library(pROC)

train <- read.csv("select feature/train radiomics features.csv",row.names = "CaseName")
test1 <- read.csv("select feature3/test1 radiomics features.csv",row.names = "CaseName")
test2 <- read.csv("select feature3/test2 radiomics features.csv",row.names = "CaseName")

model <- glm(Label ~ ., data = train,family='binomial')
summary(model)

train_predicted= predict(model,train)
test1_predicted= predict(model,test1)
test2_predicted= predict(model,test2)

trainroc = roc(train$Label,train_predicted)
test1roc = roc(test1$Label,test1_predicted)
test2roc = roc(test2$Label,test2_predicted)

auc1 = trainroc$auc
auc2 = test1roc$auc
auc3 = test2roc$auc

auc1;ci.auc(auc1)
auc2;ci.auc(auc2)
auc3;ci.auc(auc3)

roc_result_train <- coords(trainroc, "best")

#Changing the train/test1/test2, we can calculate the metrics of each set.
TP <- dim(test1[as.numeric(test1$Label)==1 & test1_predicted > roc_result_train$threshold, ])[1]
FP <- dim(test1[as.numeric(test1$Label)==0 & test1_predicted > roc_result_train$threshold, ])[1]
TN <- dim(test1[as.numeric(test1$Label)==0 & test1_predicted <= roc_result_train$threshold, ])[1]
FN <- dim(test1[as.numeric(test1$Label)==1 & test1_predicted <= roc_result_train$threshold, ])[1]

ACC <- (TP + TN) / (TP + TN + FP + FN)
Sen <- TP/(TP+FN)
Spe <- TN/(TN+FP)

#### DeLong test ####

#DeLong test for test set according to train models' predict; for example: test1_predicted= predict(logistic1,test1)

train1 <- read.csv("train Rfeature total.csv",row.names = "CaseName")
train2 <- read.csv("train Rfeature lesion.csv",row.names = "CaseName")
train3 <- read.csv("train Rfeature delta1.csv",row.names = "CaseName")
train4 <- read.csv("train Rfeature delta2.csv",row.names = "CaseName")
train5 <- read.csv("train Rfeature lesion_av.csv",row.names = "CaseName")

logistic1 <- glm(Label ~ ., data = train1,family='binomial')
train_predicted1= predict(logistic1,train1)
logistic2 <- glm(Label ~ ., data = train2,family='binomial')
train_predicted2= predict(logistic2,train2)
logistic3 <- glm(Label ~ ., data = train3,family='binomial')
train_predicted3= predict(logistic3,train3)
logistic4 <- glm(Label ~ ., data = train4,family='binomial')
train_predicted4= predict(logistic4,train4)
logistic5 <- glm(Label ~ ., data = train5,family='binomial')
train_predicted5= predict(logistic5,train5)

trainroc1 = roc(train1$Label,train_predicted1)
trainroc2 = roc(train2$Label,train_predicted2)
trainroc3 = roc(train3$Label,train_predicted3)
trainroc4 = roc(train4$Label,train_predicted4)
trainroc5 = roc(train5$Label,train_predicted5)

roc.test(trainroc1,trainroc2,method = "delong")
roc.test(trainroc1,trainroc3,method = "delong")
roc.test(trainroc1,trainroc4,method = "delong")
roc.test(trainroc1,trainroc5,method = "delong")

#### IDI test #####
#IDI test for test set according to train models' predict.

library(PredictABEL)

data_IDI <- train_predicted1%>%
  cbind(train_predicted2)%>%
  cbind(train_predicted3)%>%
  cbind(train_predicted4)%>%
  cbind(train_predicted5)%>%
  cbind(train$Label)

data_IDI <- as.matrix(data_IDI)
#replace the parameter "predrisk1" with train_predict3(4\5), we can obtain IDI results of different compare between the totol radiomic model and others. 
reclassification(data=data_IDI,cOutcome =6,predrisk1=train_predicted2,predrisk2=train_predicted1,cutoff=c(0,0.6,1)) 

#### Hybrid model building ####

library(pROC)

df <- read.csv("Clinical/zhejiang clinical.csv")
train <- subset(df, set == "1")
valid <- subset(df, set == "2")
prognosis <- read.csv("clinical/Kunming clinical.csv")
test <- read.csv("Clinical/Nantongxiangya clinical.csv")

#Rad.score is the predicted of the total radiomic model.
model <- glm(MVI ~  Pseudo.capsule + TTPVI + Peritumoral.enhance + Rad.score, data = train,family='binomial') 
model2 <- glm(MVI ~  Rad.score , data = train,family='binomial')
summary(model)
summary(model2)

train_predicted1= predict(model,train)
valid_predicted1= predict(model,valid)
test_predicted1= predict(model,test)
prognosis_predicted1 = predict(model,prognosis) # It is used for log-rank test for DFS and OS based on median in outcome cohort.

train_predicted2= predict(model2,train)
valid_predicted2= predict(model2,valid)
test_predicted2= predict(model2,test)
prognosis_predicted2 = predict(model,prognosis) # It is used for log-rank test for DFS and OS based on median in outcome cohort.

#### Metrics of Hybird model and Delong test are as same as the radiomics models section as mentioned. ####

####  DCA  ####

library(rmda)

train$train_predicted1 <- train_predicted1
train$train_predicted2 <- train_predicted2

Hybird<- decision_curve(MVI~train_predicted1,data = train, family = binomial(link ='logit'),
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95,study.design = 'case-control',
                        population.prevalence = 0.36) 

Radiomics<- decision_curve(MVI~train_predicted2,data = train, family = binomial(link ='logit'), 
                           thresholds = seq(0,1, by = 0.01),
                           confidence.intervals= 0.95,study.design = 'case-control',
                           population.prevalence= 0.36)

List<- list(Hybird,Radiomics)

plot_decision_curve(List,curve.names= c('Hybird model','Radiomics model'),
                    cost.benefit.axis =FALSE,col = c("#EE2617FF","#F2A241FF"),
                    confidence.intervals =FALSE,standardize = FALSE,ylim = c(0,0.4))

####  Calibration curve  ####

library(prodlim)
library(riskRegression)
library(ResourceSelection)

model1 <- glm(MVI ~  Pseudo.capsule + TTPVI + Peritumoral.enhance + Rad.score , data = train,family='binomial')
model2 <- glm(MVI ~  Rad.score , data = train,family='binomial')

#changing the train/valid/test2, we can plot calibration curve of each set.
xb=Score(list(model1=model1,model2=model2),MVI~1,data=train,plots="cal") 
plotCalibration(xb,brier.in.legend=T,method = "quantile",q = 6,col = c("#EE2617FF","#F2A241FF"),auc.in.legend = F,legend = T)


## 머신러닝 기반 GWAS 분석
# 참고논문 :
# 1. https://doi.org/10.3390/ijms23105538
# 2. https://doi.org/10.3390/plants12142659
# 데이터 출처 : https://github.com/miguelperezenciso/DLpipeline/tree/master/DATA

# 필요한 라이브러리 로드
library(tidyverse)
library(caret)
library(scales)
library(doMC)
registerDoMC(cores = 4)

# 데이터 불러오기
setwd("~/Google Drive/내 드라이브/Park-Jeong_Woon/개인_공부/GWAS")
X <- read.csv("wheat.X", sep = " ", header = FALSE)
y <- read.csv("wheat.Y", sep = " ", header = FALSE)

# 데이터 전처리
colnames(X) <- paste0("SNP", seq(1, ncol(X)))
colnames(y) <- paste0("Env", seq(1, ncol(y)))

X <- cbind(X, y$Env1)
colnames(X)[ncol(X)] <- "y"

# Repeated k-fold Cross Validation 설정
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
'''
    trainControl`s argument explain.
    Reference: "https://www.rdocumentation.org/packages/caret/versions/6.0-92/topics/trainControl"
    number is how many fold you want to use. You set 5, then, it means 5-fold cross cv in repeatedcv.
    repeats is only used when you set repeated CV. (repeats: For repeated k-fold cross-validation only: the number of complete sets of folds to compute)
    U set repeats as 10, that means -> 10 times repeat of 5-fold cross-validation.
    So, If you want to 1000 time iterations: set repeats as 1000. it will automately run 1000 times iteration.
    You can see iteration numbers verboseIter=TRUE

'''


# SNP 중요도를 저장할 행렬
importance_scores <- matrix(NA, nrow = 10, ncol = ncol(X) - 1)  

'''
    Importance score was calculated using varImp functions.

    e.g)
        model = train(formula, data=data, method=svr, ~~~)
        importance = varImp(model, scale=~~)
'''




# 10번 반복하여 경험적 임계값 계산
for (i in 1:10) {
  
  set.seed(i) # 매 반복마다 다른 시드를 설정하여, 교차 검증 분할이 다르게 이루어짐
  
  # SVM 모델 훈련
  svm_model <- train(y ~ ., data = X, method = "svmRadial", trControl = train_control,  
                     preProcess = c("center", "scale"))
  
  # 변수 중요도 추출
  var_imp <- varImp(svm_model, scale = FALSE)  # 중요도 추출
  
  # 중요도 데이터를 데이터프레임으로 변환
  importance_scores_df <- as.data.frame(var_imp$importance)
  
  # 중요도 값 저장
  importance_scores[i, ] <- importance_scores_df$Overall
}

# 0-100 스케일로 중요도 값 변환
importance_scaled <- apply(importance_scores, 2, function(x) rescale(x, to = c(0, 100)))

# 모든 반복을 하나의 데이터프레임으로 결합
all_importance_scores <- data.frame(importance_scaled)

# 각 SNP의 중요도 평균 계산
mean_importance_scores <- rowMeans(all_importance_scores)

# 유의미한 임계값 설정 (상위 5%의 중요도 값을 임계값으로 설정)
threshold <- quantile(mean_importance_scores, 0.95)
print(paste("Empirical significance threshold (95th percentile):", threshold))

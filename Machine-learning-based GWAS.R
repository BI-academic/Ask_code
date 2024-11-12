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

# SNP 중요도를 저장할 행렬
importance_scores <- matrix(NA, nrow = 10, ncol = ncol(X) - 1)  

# Empirical distribution 기반 threshold 계산
estimate_threshold <- function(X, model, num_iter, alpha) {
  
  # 교차검증 설정 (5-fold)
  train_control <- trainControl(method = "cv", number = 5)
  
  # 변수 중요도를 저장할 벡터 (각 반복에서 가장 높은 중요도 추출)
  importance_scores <- vector("list", num_iter)
  
  # 반복문을 통한 모델 학습
  for (i in 1:num_iter) {
    set.seed(i)  # 반복마다 시드를 다르게 설정
    train_idx <- createDataPartition(y = X$y, p = 0.8, list = FALSE)
    X_train <- X[train_idx, ]
    
    # 모델 학습 (RF 또는 SVR)
    if (model == "RF") {
      estimate <- train(y ~ ., data = X_train, method = "rf", trControl = train_control,
                        preProcess = c("center", "scale"))
    } else if (model == "SVR") {
      estimate <- train(y ~ ., data = X_train, method = "svmRadial", trControl = train_control,
                        preProcess = c("center", "scale"))
    }
    
    # 변수 중요도 추출
    var_imp <- varImp(estimate, scale = FALSE)
    
    # 각 변수의 중요도 중에서 가장 큰 값 추출하여 저장
    highest_importance_score <- max(var_imp$importance[, 1])
    
    # 중요도 점수를 리스트에 저장 (반복마다 최고 점수)
    importance_scores[[i]] <- highest_importance_score
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

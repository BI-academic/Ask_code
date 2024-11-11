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
registerDoMC(cores = 5)

# 데이터 불러오기
setwd("~/Google Drive/내 드라이브/Park-Jeong_Woon/개인_공부/GWAS")
X <- read.csv("wheat.X", sep = " ", header = FALSE)
y <- read.csv("wheat.Y", sep = " ", header = FALSE)

# 데이터 전처리
colnames(X) <- paste0("SNP", seq(1, ncol(X)))
colnames(y) <- paste0("Env", seq(1, ncol(y)))

# X에 첫 번째 환경 변수 추가
X <- cbind(X, y$Env1)
colnames(X)[ncol(X)] <- "y"

# Empirical distribution 기반 threshold 계산
estimate_threshold <- function(X, model, num_iter, alpha) {
  
  # 교차검증 설정 (5-fold)
  train_control <- trainControl(method = "cv", number = 5)
  
  # 변수 중요도를 저장할 행렬 (각 변수에 대한 중요도)
  importance_scores <- matrix(NA, nrow = num_iter, ncol = ncol(X) - 1)
  
  # 반복문을 통한 모델 학습
  for (i in 1:num_iter) {
    set.seed(i)  # 반복마다 시드를 다르게 설정
    
    # 모델 학습 (RF 또는 SVR)
    if (model == "RF") {
      estimate <- train(y ~ ., data = X, method = "rf", trControl = train_control,
                        preProcess = c("center", "scale"))
    } else if (model == "SVR") {
      estimate <- train(y ~ ., data = X, method = "svmRadial", trControl = train_control,
                        preProcess = c("center", "scale"))
    }
    
    # 변수 중요도 추출
    var_imp <- varImp(estimate, scale = FALSE)
    
    # 각 변수의 중요도 중에서 가장 큰 값 추출하여 저장
    highest_importance_score <- max(var_imp$importance[, 1])
    
    # 중요도 값을 모두 저장 (각 SNP에 대해 중요도 저장)
    importance_scores[i, ] <- highest_importance_score
  }
  
  # 중요도 점수를 0에서 100 사이로 스케일링
  importance_scaled <- apply(importance_scores, 2, function(x) rescale(x, to = c(0, 100)))
  
  # 모든 반복에서 나온 중요도 점수를 결합
  all_importance_scores <- data.frame(importance_scaled)
  
  # 유의미한 임계값 계산 (alpha 비율을 기반으로)
  threshold <- apply(all_importance_scores, 2, function(x) quantile(x, 1 - alpha))
  
  return(threshold)
}

# RF 모델을 사용하여 threshold 계산
threshold <- estimate_threshold(X, model = "RF", num_iter = 10, alpha = 0.05)
print(paste("Empirical significance threshold :", threshold))

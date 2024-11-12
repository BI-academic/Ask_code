## 머신러닝 기반 GWAS 분석
# 참고논문 :
# 1. https://doi.org/10.3390/ijms23105538
# 2. https://doi.org/10.3390/plants12142659
# 데이터 출처 : https://github.com/miguelperezenciso/DLpipeline/tree/master/DATA

# 필요한 라이브러리 로드
library(tidyverse)
library(caret)
library(scales)
library(doParallel)

# 병렬 처리 설정 (8개 코어)
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 데이터 불러오기
setwd("~/Google Drive/내 드라이브/Park-Jeong_Woon/개인_공부/GWAS")
X <- read.csv("wheat.X", sep = " ", header = FALSE)
y <- read.csv("wheat.Y", sep = " ", header = FALSE)

colnames(X) <- paste0("SNP", seq(1, ncol(X)))
colnames(y) <- paste0("Env", seq(1, ncol(y)))

X <- cbind(X, y$Env1)
colnames(X)[ncol(X)] <- "y"

estimate_threshold <- function(X, model, num_iter, alpha) {
  
  ged <- c()
  
  # trainControl에서 allowParallel = FALSE로 설정 -> caret package의 자체 병렬화 방지.
  train_control <- trainControl(method = "cv", number = 5, allowParallel = FALSE)
  
  ged <- foreach(i = 1:num_iter, .combine = c, .packages = c("caret")) %dopar% {
    set.seed(i)
    if (model == "RF") {
      estimate <- train(y ~ ., data = X, method = "rf", trControl = train_control,
                        preProcess = c("center", "scale"))
    } else if (model == "SVR") {
      estimate <- train(y ~ ., data = X, method = "svmRadial", trControl = train_control,
                        preProcess = c("center", "scale"))
    }
    max(varImp(estimate, scale = FALSE)$importance[, 1])
  }
  
  global_empirical_threshold <- quantile(ged, 1 - alpha)
  return(global_empirical_threshold)
}

# RF 모델을 사용하여 threshold 계산
threshold <- estimate_threshold(X, model = "RF", num_iter = 10, alpha = 0.05)
print(paste("Empirical significance threshold:", threshold))

# 클러스터 해제
stopCluster(cl)
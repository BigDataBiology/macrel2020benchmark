library(protr)

## PseAAC features per residue
myprops <- data.frame(AccNo = c("HydroB","HydroF","MW","pk1","pk2","pI25"),
                      A = c(0.62, -0.5, 15.0, 2.35, 9.87, 6.11),
                      R = c(-2.53, 3.0, 101.0, 2.18, 9.09, 10.76), 
                      N = c(-0.78, 0.2, 58.0, 2.18, 9.09, 10.76), 
                      D = c(-0.90, 3.0, 59.0, 1.88, 9.60, 2.98), 
                      C = c(0.29, -1.0, 47.0, 1.71, 10.78, 5.02), 
                      E = c(-0.74, 3.0, 73.0, 2.19, 9.67, 3.08), 
                      Q = c(-0.85, 0.2, 72.0, 2.17, 9.13, 5.65), 
                      G = c(0.48, 0.0, 1.0, 2.34, 9.60, 6.06), 
                      H = c(-0.40, -0.5, 82.0, 1.78, 8.97, 7.64), 
                      I = c(1.38, -1.8, 57.0, 2.32, 9.76, 6.04), 
                      L = c(1.06, -1.8, 57.0, 2.36, 9.60, 6.04), 
                      K = c(-1.50, 3.0, 73.0, 2.20, 8.90, 9.47), 
                      M = c(0.64, -1.3, 75.0, 2.28, 9.21, 5.74), 
                      F = c(1.19, -2.5, 91.0, 2.58, 9.24, 5.91),
                      P = c(0.12, 0.0, 42.0, 1.99, 10.60, 6.30), 
                      S = c(-0.18, 0.3, 31.0, 2.21, 9.15, 5.68), 
                      T = c(-0.05, -0.4, 45.0, 2.15, 9.12, 5.60), 
                      W = c(0.81, -3.4, 130.0, 2.38, 9.39, 5.88), 
                      Y = c(0.26, -2.3, 107.0, 2.20, 9.11, 5.63), 
                      V = c(1.08, -1.5, 43.0, 2.29, 9.74, 6.02))

### 1:1
x <- readFASTA("sets/1_1_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_1_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:5
x <- readFASTA("sets/1_5_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_5_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:10
x <- readFASTA("sets/1_10_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_10_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:20
x <- readFASTA("sets/1_20_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_20_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:30
x <- readFASTA("sets/1_30_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_30_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:40
x <- readFASTA("sets/1_40_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_40_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### 1:50
x <- readFASTA("sets/1_50_trainset", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "1_50_trainset.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

### Test set
x <- readFASTA("sets/test_set_1_1", seqonly = TRUE)
feat = data.frame(matrix(vector(), 24, length(x)))
for (i in (1:length(x))){
  feat[i] <- extractPAAC(x[[i]], customprops = myprops, props = c("HydroB", "MW", "pk1","pk2","pI25"), lambda = 4, w = 0.1)
}
feat <- t(feat)
write.table(feat, "test_set_1_1.tmp", sep="\t", col.names=FALSE, row.names=FALSE)

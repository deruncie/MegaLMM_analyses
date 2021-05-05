library(data.table)
library(rrBLUP)
library(tidyr)
library(rrBLUP)
library(BGLR)
source('Estimate_gcor_prediction.R')
results_dir = 'Results'

trial = as.numeric(commandArgs(t=T)[1])
if(is.na(trial)) trial = 1
foldid = as.numeric(commandArgs(t=T)[2])
if(is.na(foldid)) foldid = 1

folder = sprintf('Trial_%02d',trial)
setwd(folder)

# load data
BLUEs = fread('BLUES.csv',data.table=F)

K = fread('K.csv',data.table=F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

D = fread('D.csv',data.table=F,h=T)
rownames(D) = D[,1]
D = as.matrix(D[,-1])

X = fread('X.csv',data.table=F)
rownames(X) = X[,1]
X = as.matrix(X[,-1])

HTP = fread('HTP.csv',data.table=F)
rownames(HTP) = HTP[,1]
HTP = as.matrix(HTP[,-1])

# select training
set.seed(foldid)
BLUEs$partition = sample(c('Training','Validation'),size=nrow(BLUEs),replace=T)
BLUEs$yNA = BLUEs$Grain_Yield_BLUE
nas = BLUEs$partition == 'Validation'
BLUEs$yNA[nas] = NA

## rrBLUP - K
res_K = mixed.solve(BLUEs$yNA,K = K)
h2_BLUE = res_K$Vu/(res_K$Vu+res_K$Ve)
res_K = res_K$u[nas]

## rrBLUP - RKHS
res_RKHS = kin.blup(BLUEs,geno = 'GID',pheno = 'yNA',GAUSS = T,K = D)$pred
res_RKHS = res_RKHS[nas]

##BGLR - BL - default priors
try(dir.create('BGLR_dir'))
bglr_BL = BGLR(y = BLUEs$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F,saveAt = sprintf('BGLR_dir/foldid_%d',foldid))
bglr_BLs = X[nas,] %*% bglr_BL$ETA[[1]]$b

Knn = K[nas,nas]
sKnn = svd(Knn)

results = rbind(
  data.frame(Method = 'GBLUP_rrBLUP',
             pearson = cor(BLUEs$BLUP[nas],res_K)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = res_K),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'RKHS',
             pearson = cor(BLUEs$BLUP[nas],res_RKHS)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = res_RKHS),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'BL',
             pearson = cor(BLUEs$BLUP[nas],bglr_BLs)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = bglr_BLs),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
)
results$fold = foldid
write.csv(results,file = sprintf('%s/results_univariate_fold_%d.csv',results_dir,foldid))

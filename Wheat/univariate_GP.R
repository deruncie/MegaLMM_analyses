library(data.table)
library(rrBLUP)
library(tidyr)
library(rrBLUP)
library(BGLR)

foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)

source('data_prep_Krause.R')
source('../Simulation/Estimate_gcor_prediction.R')
try(dir.create('BGLR_dir'))

## rrBLUP - K
res_K = mixed.solve(data$yNA,K = K)$u
res_K = res_K[nas]

## rrBLUP - RKHS
D = fread('Data/D.csv',data.table=F,h=T)
rownames(D) = D[,1]
D = as.matrix(D[,-1])
res_RKHS = kin.blup(data,geno = 'GID',pheno = 'yNA',GAUSS = T,K = D)$pred
res_RKHS = res_RKHS[nas]

##BGLR - BL - default priors
bglr_BL = BGLR(y = data$yNA,ETA = list(list(X = geno2,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F,saveAt = sprintf('BGLR_dir/foldid_%d',foldid))
bglr_BLs = geno2[nas,] %*% bglr_BL$ETA[[1]]$b



results = rbind(
  data.frame(Method = 'GBLUP_rrBLUP',
             pearson = cor(data$BLUP[nas],res_K)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = res_K),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'RKHS',
             pearson = cor(data$BLUP[nas],res_RKHS)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = res_RKHS),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'BL',
             pearson = cor(data$BLUP[nas],bglr_BLs)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = bglr_BLs),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
)
results$fold = foldid
write.csv(results,file = sprintf('%s/results_univariate_fold_%d.csv',results_dir,foldid))
